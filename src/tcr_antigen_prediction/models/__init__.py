"""Models for predicting intractions between TCRs and pMHCs."""

import torch
from torch import nn


class TCRContactMapPredictor(nn.Module):
    """Model for predicting interactions between TCRs and pMHCs using structurally constrained sequences.

    Args:
        cdr_peptide_contact_maps: tensors with proportion of interacting residues between each CDR and the peptide
        cdr_mhc_contact_maps: tensors with proportion of interacting residues between each CDR and the mhc molecule
            pseudo sequence
        cdr_1_length: maximum allowed length of CDR 1s
        cdr_2_length: maximum allowed length of CDR 2s
        cdr_3_length: maximum allowed length of CDR 3s
        peptide_length: maximum allowed length of peptides
        mhc_length: length of MHC pseudo sequence used
        drop_out_rate: rate of dropout in the encoding layer during model training
        input_depth: size of input embedding space for each residue (20 for one-hot-encoding)
        learn_contact_maps: boolean to unfreeze the contact maps proportions during model training

    Attributes:
        cdr_peptide_contact_maps: tensors with proportion of interacting residues between each CDR and the peptide
        cdr_mhc_contact_maps: tensors with proportion of interacting residues between each CDR and the mhc molecule
                pseudo sequence
        cdr_embedding_layers: layers for embedding each of the 6 CDR loops
        peptide_embedding_layer: layer for embedding peptides
        mhc_embedding_layer: layer for embedding MHC pseudo sequences
        combined_fc: layer to encode relationships between TCR and the peptide or MHC
        drop_out: dropout to apply to the combined_fc layer during training
        output_fc: output layer to map latent space to binding prediction

    """

    def __init__(
        self,
        cdr_peptide_contact_maps: tuple[torch.Tensor] | None = None,
        cdr_mhc_contact_maps: tuple[torch.Tensor] | None = None,
        *,
        cdr_1_length: int = 8,
        cdr_2_length: int = 8,
        cdr_3_length: int = 24,
        peptide_length: int = 12,
        mhc_length: int = 26,
        drop_out_rate: int = 0.2,
        input_depth: int = 20,
        learn_contact_maps: bool = False,
    ) -> None:
        super().__init__()

        if cdr_peptide_contact_maps is None:
            cdr_peptide_contact_maps = tuple(
                [
                    torch.ones((cdr_length, peptide_length), dtype=torch.float32) / (cdr_length * peptide_length)
                    for _ in range(2)
                    for cdr_length in (cdr_1_length, cdr_2_length, cdr_3_length)
                ]
            )

        if cdr_mhc_contact_maps is None:
            cdr_mhc_contact_maps = tuple(
                [
                    torch.ones((cdr_length, mhc_length), dtype=torch.float32) / (cdr_length * mhc_length)
                    for _ in range(2)
                    for cdr_length in (cdr_1_length, cdr_2_length, cdr_3_length)
                ]
            )

        self.cdr_peptide_contact_maps = nn.ParameterList(
            [nn.Parameter(contact_map, requires_grad=learn_contact_maps) for contact_map in cdr_peptide_contact_maps],
        )
        self.cdr_mhc_contact_maps = nn.ParameterList(
            [nn.Parameter(contact_map, requires_grad=learn_contact_maps) for contact_map in cdr_mhc_contact_maps],
        )

        self.cdr_embedding_layers = nn.ModuleList(
            [nn.Linear(size * input_depth, size) for size in [cdr_1_length, cdr_2_length, cdr_3_length] * 2],
        )

        self.peptide_embedding_layer = nn.Linear(peptide_length * input_depth, peptide_length)
        self.mhc_embedding_layer = nn.Linear(mhc_length * input_depth, mhc_length)

        self.combined_fc = nn.Linear(
            2 * (cdr_1_length * peptide_length + cdr_2_length * peptide_length + cdr_3_length * peptide_length)
            + 2 * (cdr_1_length * mhc_length + cdr_2_length * mhc_length + cdr_3_length * mhc_length),
            512,
        )
        self.drop_out = nn.Dropout(drop_out_rate)
        self.output_fc = nn.Linear(512, 1)

    def forward(
        self,
        cdr1_alpha: torch.Tensor,  # batch_size x 8 x 20
        cdr2_alpha: torch.Tensor,  # batch_size x 8 x 20
        cdr3_alpha: torch.Tensor,  # batch_size x 24 x 20
        cdr1_beta: torch.Tensor,  # batch_size x 8 x 20
        cdr2_beta: torch.Tensor,  # batch_size x 8 x 20
        cdr3_beta: torch.Tensor,  # batch_size x 24 x 20
        peptide: torch.Tensor,  # batch_size x 12 x 20
        mhc_pseudo: torch.Tensor,  # batch_size x 26 x 20
    ) -> torch.Tensor:  # batch_size
        """Forward pass of neural network model.

        Args:
            cdr1_alpha: encoded batch of cdr1 alpha sequences (batch_size x cdr_1_length)
            cdr2_alpha: encoded batch of cdr2 alpha sequences (batch_size x cdr_2_length)
            cdr3_alpha: encoded batch of cdr3 alpha sequences (batch_size x cdr_3_length)
            cdr1_beta: encoded batch of cdr1 beta sequences (batch_size x cdr_1_length)
            cdr2_beta: encoded batch of cdr2 beta sequences (batch_size x cdr_2_length)
            cdr3_beta: encoded batch of cdr3 beta sequences (batch_size x cdr_3_length)
            peptide: encoded batch of peptide sequences (batch_size x peptide_length)
            mhc_pseudo: encoded batch of MHC pseudo sequences (batch_size x mhc_length)

        Returns:
            tensor (batch_size x 1) with the binding predictions (between 0 and 1) for each sequence in the batch

        """
        peptide_emb = torch.flatten(peptide, start_dim=1)  # batch_size x (12 * 20)
        peptide_emb = self.peptide_embedding_layer(peptide_emb)  # batch_size x 12

        mhc_emb = torch.flatten(mhc_pseudo, start_dim=1)  # batch_size x (26 * 20)
        mhc_emb = self.mhc_embedding_layer(mhc_emb)  # batch_size x 26

        peptide_interactions = []
        mhc_interactions = []
        for i, cdr in enumerate((cdr1_alpha, cdr2_alpha, cdr3_alpha, cdr1_beta, cdr2_beta, cdr3_beta)):
            cdr_emb = torch.flatten(cdr, start_dim=1)  # batch_size x (cdr_length * 20)
            cdr_emb = self.cdr_embedding_layers[i](cdr_emb)  # batch_size x cdr_length

            # batch_size x cdr_length x 12
            cdr_peptide = (cdr_emb.unsqueeze(-1) * peptide_emb.unsqueeze(-2)) * self.cdr_peptide_contact_maps[i]
            cdr_peptide = torch.flatten(cdr_peptide, start_dim=1)  # batch_size x (cdr_length * 12)
            peptide_interactions.append(cdr_peptide)

            # batch_size x cdr_length x 26
            cdr_mhc = (cdr_emb.unsqueeze(-1) * mhc_emb.unsqueeze(-2)) * self.cdr_mhc_contact_maps[i]
            cdr_mhc = torch.flatten(cdr_mhc, start_dim=1)  # batch_size x (cdr_length * 26)
            mhc_interactions.append(cdr_mhc)

        x = torch.cat((*peptide_interactions, *mhc_interactions), dim=-1)  # batch_size x 3040

        x = self.combined_fc(x)  # batch_size x 512
        x = nn.functional.relu(x)  # batch_size x 512
        x = self.drop_out(x)  # batch_size x 512

        x = self.output_fc(x)  # batch_size x 1
        prediction = nn.functional.sigmoid(x)  # batch_size x 1

        return prediction
