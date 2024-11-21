#!/usr/bin/env python
"""A Python script designed to calculate frequencies, FishersExact and Variant Effect Predictions.
"""


import os
from os.path import join

from pandas import read_csv
from snakemake import shell

__author__ = "Graeme Ford"
__credits__ = [
    "Graeme Ford",
    "Prof. Michael S. Pepper",
    "Prof. Fourie Joubert",
    "Antionette Colic",
    "Fatima Barmania",
    "Sarah Turner",
    "Megan Ryder",
]
__version__ = "1.0.0"
__maintainer__ = "Graeme Ford"
__email__ = "graeme.ford@tuks.co.za"
__status__ = "Development"


# Define constant:
config = snakemake.config
wildcards = snakemake.wildcards
params = snakemake.params
datasets = read_csv(join("input", "datasets.csv"))


# Define functions:
def directoryExists(path: str):
    """Test weather or not a directory exists. If not, create it.

    Args:
        path (str): file path of the directory to test.
    """
    if not os.path.exists(path):
        os.makedirs(path)


print("Determining Liftover requirements now...")
listed_refs: list = datasets.loc[
    datasets["dataset_name"] == wildcards.sample, "reference_genome"
].item()
if listed_refs != "GRCh38":
    shell(
        r"echo 'Liftover required. All datasets have been mapped to {}'".format(
            listed_refs
        )
    )
    if listed_refs == "GRCh37" or listed_refs == "Hg19":
        shell(r"echo 'Lifting from GRCh37 to GRCh38.'"),
        directoryExists("results/LIFTOVER")
        shell(
            r"plink2 --vcf results/PREP/{wildcards.sample}.vcf.gz --set-all-var-ids @:#\$r-\$a --allow-extra-chr --new-id-max-allele-len 40 truncate --chr 1-22 --out results/LIFTOVER/{wildcards.sample}_PREP --export vcf-4.2 bgz --output-chr chr26"
        ),
        shell(
            r"sleep 60; tabix -p vcf results/LIFTOVER/{wildcards.sample}_PREP.vcf.gz"
        ),
        shell(
            r"java -Xmx{params.mem} -jar tools/picard.jar LiftoverVcf I=results/LIFTOVER/{wildcards.sample}_PREP.vcf.gz O=results/LIFTOVER/{wildcards.sample}_LIFTOVER.vcf.gz C={params.chainFile} REJECT=results/LIFTOVER/{wildcards.sample}_REJECTED.vcf.gz R={params.ref}"
        ),
        shell(
            r"bcftools sort -m 1G -T results/LIFTOVER -O z -o results/LIFTOVER/{wildcards.sample}.vcf.gz results/LIFTOVER/{wildcards.sample}_LIFTOVER.vcf.gz"
        ),
    # TODO: Add conditionals for other human reference genome builds
    else:
        print(
            "No liftover required. Dataset {} is already mapped to GRCh38.".format(
                wildcards.sample
            )
        ),
        shell(r"touch results/LIFTOVER/{wildcards.sample}_EXCLUDE.dat"),
        shell(
            r"plink --map input/{wildcards.sample}.map --ped input/{wildcards.sample}.ped --allow-extra-chr --chr 1-22 --recode vcf --keep-allele-order --exclude {params.exclusionList} --out results/LIFTOVER/{wildcards.sample}"
        ),
    # shell("bgzip results/LIFTOVER/{wildcards.sample}.vcf"),
    shell(r"sleep 1m; tabix -f -p vcf results/LIFTOVER/{wildcards.sample}.vcf.gz"),
    shell(
        r"echo 'results/LIFTOVER/{wildcards.sample}.vcf.gz' >> results/LIFTOVER/merge.list"
    )
else:
    shell(
        r"plink2 --vcf results/PREP/{wildcards.sample}.vcf.gz --set-all-var-ids @:#\$r-\$a --allow-extra-chr --new-id-max-allele-len 400 truncate --chr 1-22 --out results/LIFTOVER/{wildcards.sample} --export vcf-4.2 bgz --output-chr chr26"
    ),
    shell(r"sleep 60; tabix -p vcf results/LIFTOVER/{wildcards.sample}.vcf.gz"),
    shell(
        r"echo 'results/LIFTOVER/{wildcards.sample}.vcf.gz' >> results/LIFTOVER/merge.list"
    )
