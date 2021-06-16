import os
import random
import re
from urllib.error import HTTPError
import pandas
import requests
import tqdm
import wget

DESIRED_CHINESE_SPRING_GENES = pandas.read_csv("239_my_only_genes_all_yes_no_dif_TE_df.txt", sep="\r", header=None,
                                               squeeze=True)
# TODO: Fix line below with Inbar.
# DESIRED_CHINESE_SPRING_GENES = DESIRED_CHINESE_SPRING_GENES.str.replace("02G", "01G")

TRANSLATION_FILE = pandas.read_csv("geneid_2_chinese.sourceid.txt", sep="\t",
                                   names=["10wheatGeneID", "chineseSpringGeneID"], header=0)
TRANSLATION_FILE["chineseSpringGeneID"] = TRANSLATION_FILE["chineseSpringGeneID"].str.rstrip(".1")


def get_regex_from_web_pages(web_page_path, regex):
    downloaded_file_path = wget.download(web_page_path)
    file_content = open(downloaded_file_path).read()
    os.remove(downloaded_file_path)
    return re.findall(regex, file_content, re.MULTILINE)


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

    # TODO: Fix line below with Inbar.
    # only_desired_genes = only_genes[only_genes["chineseSpringGeneID"].isin(DESIRED_CHINESE_SPRING_GENES)]

    # TODO change "only_genes" to "only_desired_genes"
    only_genes["startIndex"], only_genes["endIndex"] = only_genes[["startIndex", "endIndex"]].min(axis=1), \
                                                       only_genes[["startIndex", "endIndex"]].max(axis=1)

    # TODO change "only_genes" to "only_desired_genes"
    return only_genes


def open_TE(te_file_name):
    NUMBER_OF_LINES_THAT_REPRESENTS_GFF_HEADER = 13
    df = pandas.read_csv(te_file_name, "\t", skiprows=NUMBER_OF_LINES_THAT_REPRESENTS_GFF_HEADER, header=None,
                         names=["Chr", "Ver", "Type", "startIndex", "endIndex", "Irrelevant1", "Irrelevant2",
                                "Irrelevant3", "TEName"], compression=None)
    df.drop(columns=["Irrelevant1", "Irrelevant2", "Irrelevant3", "Ver", "Type", "TEName"])
    df = df[df["Chr"] != "chrUn"]

    df["startIndex"], df["endIndex"] = df[["startIndex", "endIndex"]].min(axis=1), \
                                       df[["startIndex", "endIndex"]].max(axis=1)
    return df


def has_overlap(gene_start, gene_end, tra_start, tra_end):
    if gene_start > tra_end:
        return False
    if gene_end < tra_start:
        return False

    return True


def get_te_gene_overlap(te_df, gene_df):
    results = pandas.DataFrame({"chineseSpringGeneID": [], "overlap": []})
    # for _, gene in tqdm.tqdm(gene_df.iterrows()):
    for _, gene in tqdm.tqdm(gene_df.iloc[0:5,:].iterrows()):
        chromosome_matched_transposons = te_df[te_df["Chr"] == gene["Chr"]]
        start, end = gene["startIndex"], gene["endIndex"]
        found_overlapping_transposon = False
        for _, transposon in chromosome_matched_transposons.iterrows():
            tra_start, tra_end = transposon["startIndex"], transposon["endIndex"]
            if has_overlap(start, end, tra_start, tra_end):
                found_overlapping_transposon = True
                break
        results = results.append({"chineseSpringGeneID": gene["chineseSpringGeneID"],
                                  "overlap": found_overlapping_transposon}, ignore_index=True)
    return results


def fetch_accessions_names():
    TEs_directory = r"https://webblast.ipk-gatersleben.de/downloads/wheat/TE_annotation/Wheat_11__Transposons__via_Triticeae-TE-library/"
    accession_name_and_version_regex = r'Triticum_aestivum_(.*)\.gff\.gz"'
    return get_regex_from_web_pages(TEs_directory, accession_name_and_version_regex)


def download_tdqm(path, out):
    # Streaming, so we can iterate over the response.
    response = requests.get(path, stream=True)
    total_size_in_bytes = int(response.headers.get('content-length', 0))
    block_size = 1024  # 1 Kibibyte
    progress_bar = tqdm.tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
    with open(out, 'wb') as file:
        for data in response.iter_content(block_size):
            progress_bar.update(len(data))
            file.write(data)
    progress_bar.close()
    return out


def download_TE(accession_name):
    TEs_directory = r"https://webblast.ipk-gatersleben.de/downloads/wheat/TE_annotation/Wheat_11__Transposons__via_Triticeae-TE-library/"
    file_url = "{}PGSB_Transposon_annotation-v1__Triticum_aestivum_{}.gff.gz".format(TEs_directory, accession_name)
    return download_tdqm(file_url, "{}_TE.gff".format(accession_name))


def download_gene(accession_name):
    TEs_directory = r"https://webblast.ipk-gatersleben.de/downloads/wheat/gene_projection/"
    out_file_name = "{}_gene.gff"
    try:
        file_url = "{}projectedGenes__Triticum_aestivum_{}.gff".format(TEs_directory, accession_name)

        return download_tdqm(file_url, out_file_name.format(accession_name))
    except HTTPError:
        file_url = "{}projectedGenes__Triticum_aestivum_{}_scaffolds.fasta.gz.gff".format(TEs_directory, accession_name)
        return download_tdqm(file_url, out_file_name.format(accession_name))


accessions = fetch_accessions_names()
all_accessions_results = pandas.DataFrame({"chineseSpringGeneID": [], "overlap": [], "accession": []})
# TODO: Remove the shuffle below
random.shuffle(accessions)
for accession in accessions:
    results = get_te_gene_overlap(open_TE(download_TE(accession)), open_gene(download_gene(accession)))
    results["accession"] = accession
    results.to_csv("{}.results.csv".format(accession))
    all_accessions_results = all_accessions_results.append(results)

all_accessions_results.to_csv("results.csv")
