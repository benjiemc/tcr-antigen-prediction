"""Run NetTCR predictions the data provided."""

import glob
import os
import re

import h5py
import numpy as np
import pandas as pd
import torch
from nettcr.models import NetTCR
from snakemake.script import snakemake

with h5py.File(snakemake.input.sequences) as fh:
    cdr1_alphas = fh['cdr1_alpha'][:]
    cdr2_alphas = fh['cdr2_alpha'][:]
    cdr3_alphas = fh['cdr3_alpha'][:]
    cdr1_betas = fh['cdr1_beta'][:]
    cdr2_betas = fh['cdr2_beta'][:]
    cdr3_betas = fh['cdr3_beta'][:]
    peptides = fh['peptide'][:]
    mhcs = fh['mhc_pseudo'][:]

indices = np.arange(cdr1_alphas.shape[0])

predictions = []

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print('Using device:', device)

for sub_model in glob.glob(os.path.join(snakemake.input.model, 'model_*.pt')):
    print('Loading model:', sub_model)
    model = NetTCR()
    model.load_state_dict(
        torch.load(
            sub_model,
            map_location=device,
            weights_only=True,
        )
    )
    model.to(device)

    print('Making predictions')
    model.eval()

    prediction = model(
        torch.tensor(cdr1_alphas, dtype=torch.float32),
        torch.tensor(cdr2_alphas, dtype=torch.float32),
        torch.tensor(cdr3_alphas, dtype=torch.float32),
        torch.tensor(cdr1_betas, dtype=torch.float32),
        torch.tensor(cdr2_betas, dtype=torch.float32),
        torch.tensor(cdr3_betas, dtype=torch.float32),
        torch.tensor(peptides, dtype=torch.float32),
    )

    sub_model_number = re.match(r'model_(\d+).pt', os.path.basename(sub_model)).group(1)

    predictions.append(
        pd.DataFrame(
            {
                'index': indices,
                'sub_model_number': [sub_model_number for _ in range(len(indices))],
                'prediction': prediction.squeeze(-1).detach().cpu().numpy(),
            }
        )
    )

predictions = pd.concat(predictions)
predictions = predictions.sort_values(['index', 'sub_model_number'])

predictions.to_csv(snakemake.output[0], index=False)
