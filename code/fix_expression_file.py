from pathlib import Path
import argparse
import pandas as pd

from utils import find_subset

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, default=None)
    parser.add_argument("-r", "--region", type=str, default=None)
    args = parser.parse_args()

    if args.input is not None:
        datafile = Path(args.input)
        gnomad = pd.read_csv(datafile)
    else:
        datafile = (
            Path(__file__)
            .resolve()
            .parents[1]
        .joinpath("data", "rawData", "2020_07_22_ALS_rare_variants_0.01_with_gnomAD_NFE_MPCvalues.txt")
    )
        gnomad = pd.read_csv(datafile, sep="\t")

    if args.region is not None:
        region = args.region
    else:
        region = ""
    
    infile = Path("~/Downloads/merged_with_means.csv")
    outfile = Path("../data/proteinDomain/{}variants_mean_exp.csv".format(region))

    df = pd.read_csv(infile, usecols=["variant", "mean_brain", "mean_all"])

    exps = ["high", "medium", "low"]
    names = ["brain", "all"]

    for name in names:
        for exp in exps:
            filename = Path("../data/proteinDomain/{}{}_{}.csv".format(region, name, exp))

            df1 = find_subset(df, "mean_{}".format(name), exp)
            var_names = df1.variant
            df1 = gnomad.loc[gnomad['variant'].isin(var_names)]

            df1.to_csv(filename, index=False)
            

    df.to_csv(outfile, index=False)

if __name__ == "__main__":
    main()