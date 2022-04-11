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
    name = str(args.input)
else:
    name = "RV_brain_medium"

if args.type is not None:
    mutation = str(args.type)
else:
    mutation = "missense_mpc_2"

'''
docs_file = (
        Path(__file__)
        .resolve()
        .parents[1]
    .joinpath("data", "proteinDomain", "{}.csv".format(name))
)
'''
docs_file = (Path(__file__)
            .resolve()
            .parents[1]
            .joinpath("data", "proteinDomain", "variants_in_protein_domain.csv"))

df = pd.read_csv(docs_file, usecols=["variant", "MPC", "consequence", "gene", "case_control", "sample"], index_col="variant")
print("Loaded df")

if mutation == "missense_mpc_2":
    exp_df = get_case_control_per_variant(df, mpc=True).reset_index()
else:
    exp_df = get_case_control_per_variant(df).reset_index()

#exp_df.to_csv("../data/proteinDomain/{}_{}_case_control.csv".format(name,mutation), index=False)
#exp_df.to_csv("../data/proteinDomain/exp_{}_case_control_filtered.csv".format(level), index=False)

print("Calculating fishers for {}".format(mutation))
fishers = find_subset(exp_df, "case_{}".format(mutation), 3864, "<=")
fishers = find_subset(exp_df, "control_{}".format(mutation), 7839, "<=")
fishers = fishers_test(fishers, mutation, 3864, 7839)


final = exp_df.join(fishers, how="right")

final.to_csv("../data/proteinDomain/exp_{}_{}_fishers_case_control.csv".format(name, mutation), index=False)