import pandas as pd

def main():

    l = ["noFilter", "proteinDomain", "allHigh", "allMedium", "allLow", "brainHigh", "brainMedium", "brainLow", \
    "pd_allHigh", "pd_allMedium", "pd_allLow", "pd_brainHigh", "pd_brainMedium", "pd_brainLow"]

    dfs = []

    for t in l:
        name = "RV_" + t
        dfs.append(pd.read_csv("../data/proteinDomain/{}_count_lr.csv".format(name)))

    df = pd.concat(dfs)
    df.to_csv("../data/proteinDomain/all_lr.csv", index=False)

if __name__ == "__main__":
    main()