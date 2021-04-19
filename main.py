import wget
import re
import os
import pandas


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


# download_genes()
# download_TEs()
df = pandas.read_csv("D:\Code\WheatGenome\projectedGenes__Triticum_aestivum_ArinaLrFor_v3.0.gff", "\t", header=None,
                     names=["Chr", "Ver", "Type", "startIndex", "endIndex", "Irrelevant1", "Irrelevant2", "Irrelevant3",
                            "geneID"])
df = df.drop(columns=["Ver", "Irrelevant1", "Irrelevant2", "Irrelevant3"])
only_genes = df[df["Type"] == "gene"]
only_genes["startIndex"], only_genes["endIndex"] = only_genes[["startIndex", "endIndex"]].min(axis=1), \
                                                   only_genes[["startIndex", "endIndex"]].max(axis=1)
only_genes

