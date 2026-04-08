import scanpy as sc
import infercnvpy as cnv
import pandas as pd

adata = sc.read_10x_h5("/home/sduan/projects/immune-atlas/data/processed/scrup/harmonized.counts.h5")

metadata = pd.read_csv("/home/sduan/projects/immune-atlas/data/processed/scrup/harmonized.cell_metadata.integrated.csv")

metadata = metadata.set_index("barcode").loc[adata.obs_names]
adata.obs = metadata

cnv.io.genomic_position_from_gtf(
    gtf_file="/home/sduan/data/scrup/reference/GRCh38/v48/sources/filtered.gtf",
    adata=adata
)

reference_types = [
        "B cells",
        "Innate lymphoid cells"
]

adata_sub = adata[adata.obs["cell_type_level1"].isin(reference_types + ["Epithelial cells"])].copy()

cnv.tl.infercnv(
    adata_sub,
    reference_key="cell_type_level1",
    reference_cat=reference_types,
    window_size=250,
)

import os

# Define output directory and ensure it exists.
outdir = "/home/sduan/projects/immune-atlas/data/processed/scrup/cnv"
os.makedirs(outdir, exist_ok=True)

# Define output file path.
outfile = os.path.join(outdir, "harmonized_infercnv_epithelial_vs_reference.h5ad")

# Write the AnnData object with CNV results.
adata_sub.write(outfile)

print(f"InferCNV results have been saved to: {outfile}")


import pandas as pd

# Convert CNV matrix to a DataFrame.
cnv_df = pd.DataFrame(
    adata_sub.obsm["X_cnv"],
    index=adata_sub.obs_names
)

# Save CNV matrix.
cnv_matrix_file = os.path.join(outdir, "harmonized_infercnv_X_cnv.tsv.gz")
cnv_df.to_csv(cnv_matrix_file, sep="\t", compression="gzip")

print(f"CNV matrix has been saved to: {cnv_matrix_file}")



