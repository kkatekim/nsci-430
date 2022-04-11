import pandas as pd
import argparse
from pathlib import Path
import hail as hl

from utils import fishers_test, get_case_control_per_variant, find_subset

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, default=None)
parser.add_argument("-t", "--type", type=str, default=None)
args = parser.parse_args()

if args.input is not None:
    docs_file = Path(args.input)
else:
    docs_file = (
        Path(__file__)
        .resolve()
        .parents[1]
    .joinpath("data", "rawData", "total_als_file.csv")
)

if args.type is not None:
    level = str(args.type)
else:
    level = "high"

df = pd.read_csv(docs_file, usecols=["variant", "MPC", "consequence", "gene", "case_control", "sample", "mean_brain", "mean_all"], index_col="variant")
#df = df[df["gene"] != "TTN"]
#df = find_subset(df, "MPC", 10, "<=")
print("Loaded df")

for exp in ["mean_brain", "mean_all"]:
    exp_df = find_subset(df, exp, level)

    exp_df = get_case_control_per_variant(exp_df, mpc=True).reset_index()
    exp_df.to_csv("../data/proteinDomain/exp_{}_case_control_all.csv".format(level), index=False)
    #exp_df.to_csv("../data/proteinDomain/exp_{}_case_control_filtered.csv".format(level), index=False)

    print("Got subset for {}".format(exp))
    print(exp_df.head())

    fishers_list = []

    for variant in ["synonymous", "missense", "PTVs"]:
        print("Calculating fishers for {}".format(variant))
        fishers = find_subset(exp_df, "case_{}".format(variant), 3864, "<=")
        fishers = find_subset(exp_df, "control_{}".format(variant), 7839, "<=")
        fishers = fishers_test(fishers, variant, 3864, 7839)
        fishers_list.append(fishers)

    exp_df = exp_df.join(pd.concat(fishers_list, axis=1), how="right")

    exp_df.to_csv("../data/proteinDomain/exp_{}_{}_fishers.csv".format(exp, level), index=False)