import pandas as pd
from pathlib import Path
import argparse

import utils


def split_data():

    file = Path("../data/rawData/2020_07_22_ALS_rare_variants_0.01_with_gnomAD_NFE_MPCvalues.txt")
    df = pd.read_csv(file, delimiter="\t")
    df = df.drop_duplicates(subset="variant")

    keep = ["variant", "chr", "pos", "gene"]
    columns = [x for x in df.columns.tolist() if x not in keep]
    df = df.drop(columns, axis=1)

    start = 0
    for i in range(250000, 1243022, 250000):

        df1 = df.iloc[start:i].reindex(columns=df.columns.tolist())
        df1.to_csv("../data/rawData/subset_{}.csv".format(int(i/250000)), index=False)
        start = i+1

    df1 = df.iloc[start:1243022].reindex(columns=df.columns.tolist())
    df1.to_csv("../data/rawData/subset_5.csv")


def run_data(variants, protein, out_file):
    ''' takes as input 2 dfs (variants is a list of all gene variants and protein is a list of all protein domains)'''

    #TODO: switch everything to iteritems from df

    chr_v = variants.chr.tolist()
    pos_v = variants.pos.tolist()

    chr_short = [utils.get_chr_num(x) for x in protein["chrom"].tolist()]
    protein["chr_short"] = pd.Series(chr_short)

    crnt_chr = str(chr_v[0])
    df = utils.find_subset(protein, "chr_short", crnt_chr)

    start = df.chr_start.tolist()
    end = df.chr_end.tolist()
    gene = df.gene.tolist()

    variants_in_protein_domain = pd.DataFrame()

    col = []

    for i in range(len(chr_v)):

        if str(chr_v[i]) != crnt_chr:
            crnt_chr = str(chr_v[i])
            df = utils.find_subset(protein, "chr_short", crnt_chr)
            start = df.chr_start.tolist()
            end = df.chr_end.tolist()
            gene = df.gene.tolist()
        
        print(i, chr_v[i], crnt_chr)
        for j in range(df.shape[0]):
            if pos_v[i] >= start[j] and pos_v[i] <= end[j]:
                variants_in_protein_domain = variants_in_protein_domain.append(variants.iloc[i])
                col.append(
                    gene[j] +
                    ":" +
                    crnt_chr +
                    ":" + 
                    str(int(start[j])) + 
                    "-" + 
                    str(int(end[j])) +
                    ":" + 
                    str(int(pos_v[i]))
                )

    variants_in_protein_domain = variants_in_protein_domain.reindex(columns=variants.columns.tolist())
    variants_in_protein_domain["gene:chr:domain:pos"] = col
    variants_in_protein_domain.to_csv(out_file, index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--docs", type=str, default=None)
    parser.add_argument("-o", "--output", type=str, default=None)
    args = parser.parse_args()

    if args.docs is not None:
        file = Path(args.docs)
        variants = pd.read_csv(file)
    else:
        file = Path("../data/rawData/2020_07_22_ALS_rare_variants_0.01_with_gnomAD_NFE_MPCvalues.txt")
        variants = pd.read_csv(file, delimiter="\t")

    if args.output is not None:
        out_file = Path(args.output)
    else:
        out_file = Path("../data/variants_in_domain.csv")

    file_p = Path("../data/datasets/prot2hg_1938_112019.csv")
    protein = pd.read_csv(file_p, delimiter=";")

    #split_data()
    run_data(variants, protein, out_file)
    # merge_data()

if __name__ == "__main__":
    main()
    