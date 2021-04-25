import wget
import re
import os
import pandas
import tqdm

DESIRED_CHINESE_SPRING_GENES = pandas.read_csv("239_my_only_genes_all_yes_no_dif_TE_df.txt", sep="\r", header=None,
                                               squeeze=True)
# DESIRED_CHINESE_SPRING_GENES = DESIRED_CHINESE_SPRING_GENES.str.replace("02G", "01G")

TRANSLATION_FILE = pandas.read_csv("geneid_2_chinese.sourceid.txt", sep="\t",
                                   names=["10wheatGeneID", "chineseSpringGeneID"], header=0)
TRANSLATION_FILE["chineseSpringGeneID"] = TRANSLATION_FILE["chineseSpringGeneID"].str.rstrip(".1")


def get_full_file_names_in_directory(directory_path, file_name_regex_pattern):
    downloaded_file_path = wget.download(directory_path)
    file_content = open(downloaded_file_path).read()

    matches = re.findall(file_name_regex_pattern, file_content, re.MULTILINE)
    striped_matches = [directory_path + match.replace("\"", "") for match in matches]
    os.remove(downloaded_file_path)
    return striped_matches


def download_genes():
    regex = r"\"projectedGenes.*\.gff\""
    link = r'https://webblast.ipk-gatersleben.de/downloads/wheat/gene_projection/'
    files_paths = get_full_file_names_in_directory(link, regex)
    for file_path in files_paths:
        wget.download(file_path)


def download_TEs():
    regex = r"\"PGSB.*\.gff\.gz\""
    link = r'https://webblast.ipk-gatersleben.de/downloads/wheat/TE_annotation/Wheat_11__Transposons__via_Triticeae-TE-library/'
    files_paths = get_full_file_names_in_directory(link, regex)
    for file_path in files_paths:
        wget.download(file_path)


def open_gene(gene_file_name):
    df = pandas.read_csv(gene_file_name, "\t",
                         header=None,
                         names=["Chr", "Ver", "Type", "startIndex", "endIndex", "Irrelevant1", "Irrelevant2",
                                "Irrelevant3",
                                "10wheatGeneID"])

    df = df.drop(columns=["Ver", "Irrelevant1", "Irrelevant2", "Irrelevant3"])
    only_genes = df[df["Type"] == "gene"]
    only_genes["10wheatGeneID"] = only_genes["10wheatGeneID"].str.lstrip("ID=")

    only_genes = only_genes.join(TRANSLATION_FILE.set_index("10wheatGeneID"), "10wheatGeneID")

    # TODO: Fix line below with inbar.
    # only_desired_genes = only_genes[only_genes["chineseSpringGeneID"].isin(DESIRED_CHINESE_SPRING_GENES)]

    only_genes["startIndex"], only_genes["endIndex"] = only_genes[["startIndex", "endIndex"]].min(axis=1), \
                                                       only_genes[["startIndex", "endIndex"]].max(axis=1)

    return only_genes


def open_TE(te_file_name):
    NUMBER_OF_LINES_THAT_REPRESENTS_GFF_HEADER = 13
    return pandas.read_csv(te_file_name, "\t", skiprows=NUMBER_OF_LINES_THAT_REPRESENTS_GFF_HEADER, header=None,
                           names=["Chr", "Ver", "Type", "startIndex", "endIndex", "Irrelevant1", "Irrelevant2",
                                  "Irrelevant3", "TEName"]).drop(
        columns=["Irrelevant1", "Irrelevant2", "Irrelevant3", "Ver", "Type"])


def has_overlap(gene_start, gene_end, tra_start, tra_end):
    if gene_start > tra_end:
        return False
    if gene_end < tra_start:
        return False

    return True


a = open_TE("PGSB_Transposon_annotation-v1__Triticum_aestivum_ArinaLrFor_v3.0.gff")

b = open_gene(r"D:\Code\WheatGenome\projectedGenes__Triticum_aestivum_ArinaLrFor_v3.0.gff")

results = pandas.DataFrame({"chineseSpringGeneID": [], "overlap": []})

for _, gene in tqdm.tqdm(b.iterrows()):
    chromosome_matched_transposons = a[a["Chr"] == gene["Chr"]]
    start, end = gene["startIndex"], gene["endIndex"]
    found_overlapping_transposon = False
    for _, transposon in chromosome_matched_transposons.iterrows():
        tra_start, tra_end = transposon["startIndex"], transposon["endIndex"]
        if has_overlap(start, end, tra_start, tra_end):
            found_overlapping_transposon = True
            break
    results = results.append({"chineseSpringGeneID": gene["chineseSpringGeneID"],
                              "overlap": found_overlapping_transposon}, ignore_index=True)
    results.to_csv("results.csv")


# download_genes()
# download_TEs()

### DESIRED FLOW
# TODO 1. REading all dataframes of genes
# TODO 2. FIltering unwanted genes from each dataframe
# TODO 3. Combining all dataframes into one dataframe and adding a column indicating the source of each genome
#   * Naming the source dataframe using the regex: "_aestivum_(.*)_.*\.gff"

# ---Continuing with TE in GENE Flow---
