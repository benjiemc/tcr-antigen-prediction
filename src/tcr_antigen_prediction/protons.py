from subprocess import Popen, PIPE


def reprotonate(path: str, out_path: str):
    '''Remove hydrogens (if any) and re-protonate a structure.'''
    # Remove hydrogens
    args = ["reduce", "-Trim", path]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p2.communicate()

    with open(out_path, "w") as outfile:
        outfile.write(stdout.decode('utf-8').rstrip())

    # Re-add hydrogens
    args = ["reduce", "-HIS", out_path]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p2.communicate()

    with open(out_path, "w") as outfile:
        outfile.write(stdout.decode('utf-8'))

