import argparse
from pathlib import Path
import pandas as pd
import numpy as np

import utils

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--name", type=str, default=None)
    args = parser.parse_args()

    name = args.name 

    counts_file = Path(__file__).resolve().parents[1].joinpath("data", "proteinDomain", "{}_samples_count.csv".format(name))
    lr_file = Path(__file__).resolve().parents[1].joinpath("data", "proteinDomain", "{}_lr.csv".format(name))
    outfile = Path(__file__).resolve().parents[1].joinpath("data", "proteinDomain", "{}_count_lr.csv".format(name))

    counts = pd.read_csv(counts_file, usecols=np.arange(1,7,1))
    lr = pd.read_csv(lr_file)

    variant_types = ["benign_count", "damaging_count", "mis_count", "mpc>=2_count", "ptv_count", "syn_count"]

    counts.loc["count"] = np.zeros(6)

    for vt in variant_types:
        col = counts.filter(like=vt, axis=1)
        counts.loc["count", col.columns[0]] = int(col.sum()[0])

    # add counts to lr df
    lr["count"] = counts.loc["count"].tolist()

    # reorder columns so count is 2nd
    cols = list(lr.columns)
    cols.insert(1, cols.pop(cols.index("count")))
    lr = lr[cols]

    lr.to_csv(outfile, index=False)

if __name__ == "__main__":
    main()