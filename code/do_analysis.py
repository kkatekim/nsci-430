import argparse
from pathlib import Path
import pandas as pd

import utils

def get_sample_counts(df, samples, name, save=False, verbose=False):
    rows = []

    for sample in samples:
        df1 = df[df["sample"] == sample]
        mis = df1[df1["consequence"] == "missense"]

        new_row = {}
        new_row["sample"] = sample
        new_row["{}_syn_count".format(name)] = len(df1[df1["consequence"] == "synonymous"])
        new_row["{}_mis_count".format(name)] = len(mis)
        new_row["{}_mis_mpc>=2_count".format(name)] = len(mis[mis["MPC"] >= 2])
        new_row["{}_mis_damaging_count".format(name)] = len(mis[ (mis["sift_prediction"] == "deleterious") & (mis["polyphen2_prediction"] == "probably_damaging") ])
        new_row["{}_mis_benign_count".format(name)] = len(mis[ (mis["sift_prediction"] == "tolerated") & (mis["polyphen2_prediction"] == "benign") ])
        new_row["{}_ptv_count".format(name)] = len(df1[df1["consequence"] == "lof"])
        new_row["{}_sum".format(name)] = new_row["{}_syn_count".format(name)] + new_row["{}_mis_count".format(name)] + new_row["{}_ptv_count".format(name)]
        rows.append(new_row)
        if verbose:
            #print(len(mis), new_row["RV_noFilter_mis_damaging_count"] + new_row["RV_noFilter_mis_benign_count"])
            print(new_row)
    df1 = pd.DataFrame(rows)
    if save:
        df1.to_csv(Path(__file__).resolve().parents[1].joinpath("data", "rawData", "{}_samples_count.csv".format(name)), index=False)
    return df1

def merge_info(df, name, save=False, verbose=False):
    info_path = Path(__file__).resolve().parents[1].joinpath("data", "rawData", "2021_06_17_ALS_log_reg_Untouched.txt")
    info_df = pd.read_csv(info_path, sep="\t")
    merged = info_df.merge(df, on="sample")
    if save:
        merged.to_csv(Path(__file__).resolve().parents[1].joinpath("data", "proteinDomain", "{}_log_reg_info.csv".format(name)), index=False)
    if verbose:
        print(merged)
    return merged
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, default=None)
    parser.add_argument("-o", "--output", type=str, default=None)
    parser.add_argument("-n", "--name", type=str, default=None)
    args = parser.parse_args()

    if args.input is not None:
        infile = Path(args.input)
        df = pd.read_csv(infile)
    else:
        infile = (
            Path(__file__)
            .resolve()
            .parents[1]
        .joinpath("data", "rawData", "2020_07_22_ALS_rare_variants_0.01_with_gnomAD_NFE_MPCvalues.txt")
    )
        df = pd.read_csv(infile, sep="\t")

    if args.name is not None:
        name = args.name
    else:
        name = "RV_noFilter"
        
    if args.output is not None:
        outfile = args.output
    else:
        outfile = "../data/proteinDomain/{}_samples_count.csv".format(name)

    

    samples = list(dict.fromkeys(df["sample"].tolist()))

    s_df = get_sample_counts(df, samples, name)
    s_df.to_csv(outfile, index=False)
    
    #s_df = pd.read_csv("../data/proteinDomain/exome_samples_count.csv")

    #merged_df = merge_info(s_df, name, save=True)
    #df = find_subset(df, "MPC", 10, "<=")

if __name__ == "__main__":
    main()