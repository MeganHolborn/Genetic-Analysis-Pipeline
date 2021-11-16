import pandas as pd
from os.path import join

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

cluster = pd.read_excel(join("config", snakemake.config["cluster"]["file"]))
cluster["FID"] = cluster["ID"]
cluster[["ID", "FID", snakemake.wildcards.cluster]].to_csv(
    join("results", "REFERENCE", "cluster_{}.txt".format(snakemake.wildcards.cluster)),
    sep="\t",
    index=False,
)
