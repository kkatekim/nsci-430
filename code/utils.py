import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import hail as hl

'''Module with helper functions used in both projects (and for mt + tRNA).'''

def find_subset(df, column_name, factor, condition="="):
    '''Returns df subsetted by factor in specified column (==, !=, <, >, <=, >=).'''

    if condition == "=":
        return df[df[column_name] == factor]
    
    elif condition == "!=":
        return df[df[column_name] != factor]
    
    elif condition == "<":
        return df[df[column_name] < factor]

    elif condition == ">":
        return df[df[column_name] > factor]

    elif condition == "<=":
        return df[df[column_name] <= factor]

    elif condition == ">=":
        return df[df[column_name] >= factor]


def get_geneIDs(file, gene_type, column_name, sheet=None):
    '''Extracts genes names as series, writes them into csv, and returns it'''
    
    #out_file = Path("../data/{}_gene_IDs.csv".format(gene_type))

    if sheet is None:
        df = pd.read_excel(file)
    else:
        df = pd.read_excel(file, sheet_name=sheet)

    gene_IDs = df[column_name].dropna().reset_index(drop=True)
    #gene_IDs.to_csv(out_file, index=False)

    return gene_IDs


def get_genes_of_interest(dataset_df, gene_df):
    '''Finds genes of interest, writes them into csv, and returns df'''

    gene_file = Path(gene_file)
    out_file = Path("../data/{}_gene_variants.csv".format(gene_type))

    gene = get_geneIDs(gene_file, gene_type, column_name, sheet)
    gene.rename("gene", inplace=True)

    gene = df.merge(gene, how="inner", on="gene")
    gene.to_csv(out_file, index=False)

    return gene


def calculate_fishers(case_carrier, control_carrier, case_noncarrier, control_noncarrier):
    '''Runs fisher's test one one gene and returns a list with pvalue, OR, and CI.'''
    result = hl.eval(hl.fisher_exact_test(case_carrier, case_noncarrier, control_carrier, control_noncarrier))
    return [result["p_value"], result["odds_ratio"], result["ci_95_lower"], result["ci_95_upper"]]


def fishers_test(df, variant, n_cases, n_controls):
    '''Runs fisher's test on specified variant and returns df with pval, OR, and CI added.'''
    case = "case_" + variant
    control = "control_" + variant

    fishers_df = pd.DataFrame(index=df.index)
    fishers_df["case_carrier"] = df[case].values
    fishers_df["control_carrier"] = df[control].values
    fishers_df["case_noncarrier"] = n_cases - df[case].values
    fishers_df["control_noncarrier"] = n_controls - df[control].values 
    fishers_df = fishers_df.astype(np.int32)

    # col names
    col_p = "pval_" + variant
    col_or = "OR_" + variant
    col_lowci = "lowci_" + variant
    col_highci = "highci_" + variant
    
    fishers_df[[col_p, col_or, col_lowci, col_highci]] = fishers_df.apply(lambda x: calculate_fishers(
                            x.case_carrier, x.control_carrier, x.case_noncarrier, x.control_noncarrier),
                            axis=1, result_type="expand")

    return fishers_df.drop(labels=["case_carrier", "control_carrier", "case_noncarrier", "control_noncarrier"], axis=1)


def generate_table(df, group_list, index, columns):
    '''Groups data together based on columns of interest.'''
    # ex. of how to use: print(generate_table(df, ["gene", "case_control"], "gene", "case_control"))

    new_df = pd.DataFrame({"count": df.groupby(group_list).size()}).reset_index()
    return new_df.pivot(index=index, columns=columns).droplevel(0, axis=1)


def rename_and_reorder_columns(df, col_names, col_order):
    '''If column exists in df1, adds to df2'''
    
    df = df.rename(columns={k: v for k, v in col_names.items() if k in df})

    for col in col_order:
        if col not in df.columns:
            df[col] = [0] * len(df.index)

    return df[col_order]

def get_case_control_per_variant(df, mpc=False, file_name=None):
    '''Returns df with number of cases and controls per gene (and number of samples for mutation type).'''

    # find case count for each variant
    case_control = pd.DataFrame({"count": df.groupby(["gene", "consequence", "case_control"]).size()}).reset_index()
   
    # merge consequence and case_control column and add to the df
    case_control["consequence_case_control"] = case_control["case_control"].str.lower() + "_" + case_control["consequence"]

    # remove unneeded labels, switch consequence_case_control to column header
    # then drop multilevel index (ie count and consequence_case_control)
    # rename columns as needed and reorder them
    case_control = (case_control.drop(labels=["consequence", "case_control"], axis=1)
                                .pivot(index="gene", columns="consequence_case_control")
                                .droplevel(0, axis=1)     
                    )

    if not mpc:

        case_control = rename_and_reorder_columns(case_control,
                        {"case_lof": "case_PTVs", "control_lof": "control_PTVs"},
                        ["case_synonymous", "control_synonymous", "case_missense", "control_missense", "case_PTVs", "control_PTVs"])
        
        # get counts per mutations
        #mutation_count = pd.DataFrame({"count": df.groupby(["gene", "consequence"]).size()}).reset_index()
        mutation_count = generate_table(df, ["gene", "consequence"], "gene", "consequence")
        mutation_count = rename_and_reorder_columns(mutation_count, 
                        {"missense": "sample_missense", "lof": "sample_PTVs", "synonymous": "sample_synonymous"}, 
                        ["sample_synonymous", "sample_missense", "sample_PTVs"])

    else:

        case_control = rename_and_reorder_columns(case_control, 
                        {"case_missense": "case_missense_mpc_2", "control_missense": "control_missense_mpc_2"}, 
                        ["case_missense_mpc_2", "control_missense_mpc_2"])
        
        mutation_count = generate_table(df, ["gene", "consequence"], "gene", "consequence")
        mutation_count = rename_and_reorder_columns(mutation_count, 
                        {"missense": "sample_missense_mpc_2"},
                        ["sample_missense_mpc_2"])
 
     # merge, fill NaN as 0, and convert everything to type int
    complete_df = pd.merge(case_control, mutation_count, on="gene").fillna(0).astype(int)

    # if given a gene_type, write output to file with gene_type as name
    if file_name is not None:
        out_file = Path("../data/{}_case_control_variants_per_gene.csv".format(file_name))
        complete_df.to_csv(out_file)

    return complete_df 


def get_chr_num(chrom):
    '''Returns chrom number extracted from hg19 chrom string.''' 
    if "_" in chrom:
        return chrom[3:chrom.index("_")]
    return chrom[3:]


def plot_qq(df, name, variant):
    '''plots QQ graph of expected vs observed using plt scatter plot. doesn't include CI rn'''

    # finds relevant column names
    col_p = "pval_" + variant
    col_low = "lowci_" + variant
    col_high = "highci_" + variant
    
    # removes all NaN and sort
    cleaned = df.dropna(subset=[col_p, col_low, col_high], how="all").sort_values(col_p)

    # -log to get observed values
    p_list = cleaned[col_p].to_numpy()
    obs = -1 * np.log10(p_list)
    exp = -1 * np.log10(np.arange(1, p_list.shape[0]+1) / p_list.shape[0])
    
    fig, ax = plt.subplots()
    ax.scatter(exp, obs)
    xpoints = ypoints = ax.get_xlim()
    ax.plot(xpoints, ypoints, color='black', scalex=False, scaley=False)
    plt.title("{} {} Q-Q Plot".format(name, variant))
    plt.xlabel("Expected value")
    plt.ylabel("Observed value")
    plt.savefig("../data/proteinDomain/figures/{}_{}_QQ.png".format(name, variant))
    plt.show()

    return cleaned


def find_significant_genes(df, variant, alpha=0.05, file_name=None):
    '''find genes less than p value based on mutation variant'''
    
    # multiple test correction, divide alpha by # of genes
    adjusted_p = alpha/df.shape[0]

    col = "pval_" + variant
    df.sort_values(col, ignore_index=True, inplace=True)

    significant_genes = find_subset(df, col, adjusted_p, "<=")
    significant_genes = find_subset(significant_genes, "pval_synonymous", adjusted_p, ">")

    if file_name is not None:
        out_file = Path("../data/datasets/{}_{}_significant_genes.csv".format(file_name, variant))
        significant_genes.to_csv(out_file, index=False)
    
    return significant_genes


def find_sample_variants_for_genes(df, filename, name, variant):

    all_genes = pd.read_csv(filename, usecols=["gene", "variant", "consequence", "case_control"])

    df = df[["gene", "pval_{}".format(variant), "OR_{}".format(variant), "lowci_{}".format(variant), "highci_{}".format(variant)]]
    signficant_gene_variants = pd.merge(all_genes, df, how="inner", on="gene")
    signficant_gene_variants = find_subset(signficant_gene_variants, "consequence", variant, "=")

    signficant_gene_variants.to_csv("../data/{}_{}_significant_variants.csv".format(name,variant), index=False)

    return signficant_gene_variants


def create_variant_count_by_case_control(df, name, mutation):

    df_fishers = df.drop(["gene", "case_control", "consequence"], axis=1)
    grouped = pd.DataFrame({"variant_count": df.groupby(["gene", "variant"]).size()}).reset_index()

    df_list = []

    for gene in set(df["gene"].tolist()):

        subset = find_subset(df, "gene", gene, "=")

        for variant in set(subset["variant"].tolist()):
            variant_subset = find_subset(subset, "variant", variant, "=")

            variant_list={}
            variant_list["variant"] = variant
            case = find_subset(variant_subset, "case_control", "Case", "=")
            control = find_subset(variant_subset, "case_control", "Control", "=")
            variant_list["case"] = case.shape[0]
            variant_list["control"] = control.shape[0]

            df_list.append(variant_list)

        #case_control = pd.DataFrame({"case_control_count": df.groupby(["variant", "case_control"]).size()}).reset_index()

    all_case_control = pd.DataFrame(df_list)

    counts = pd.merge(grouped, all_case_control, on="variant")
    df = pd.merge(counts, df_fishers, on="variant").drop_duplicates("variant")
    df.to_csv("../data/{}_{}_signficant_genes_count.csv".format(name, mutation), index=False)