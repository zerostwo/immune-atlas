# Uploading the immune-atlas source dataset to Zenodo

This document describes how to use `scripts/r/upload_data.R` to publish the immune-atlas source dataset to Zenodo with `zen4R`.

The upload script defaults to a dry run. It will only create or update a Zenodo deposition when explicitly run with `--execute`.

## Files to upload

The current source dataset consists of:

- `data/raw/reference/reference.raw.qs`
- `data/raw/reference/bpcells.tar.gz`

For compatibility, the upload script will look for either:

- `data/raw/reference/raw.qs`
- `data/raw/reference/reference.raw.qs`

In the current repository, the actual file present is:

- `data/raw/reference/reference.raw.qs`

## Data provenance

The source/reference dataset is generated from:

- `scripts/r/02.building_and_annotating_the_atlas_core/01_build_reference_source_data.R`

That script writes the serialized reference object to:

- `data/raw/reference/reference.raw.qs`

and writes the BPCells-backed matrices to:

- `data/raw/reference/bpcells/`

The BPCells directory is then archived as:

- `data/raw/reference/bpcells.tar.gz`

The upload script does not rebuild the data. It only validates metadata and uploads the files to Zenodo.

## Dataset summary

The current `reference.raw.qs` object contains:

- `4,255,241` cells
- `104` donors
- `772` samples
- `10` source studies
- `32` tissue and organ categories
- `54,872` RNA features
- `269` ADT features across `1,150,344` cells
- `5` assay chemistries: `10x3'v2`, `10x3'v3`, `10x3'v3.1`, `10x5'v1`, `10x5'v2`

Included tissue and organ categories:

- `Bile duct, Bladder, Blood, Bone marrow, Ear, Esophagus, Eye, Fat, Heart, Kidney, Large intestine, Liver, Lung, Lymph node, Mammary, Muscle, Omentum, Ovary, Pancreas, Prostate, Salivary gland, Skin, Small intestine, Spleen, Stomach, Testis, Thymus, Tongue, Tonsil, Trachea, Uterus, Vasculature`

Study-level cell counts:

- `TabulaSapiens2022`: `2,280,879`
- `StevenBWells2025`: `1,065,585`
- `EMadissoon2020`: `284,573`
- `CDominguezConde2022`: `163,992`
- `ShuaiHe2020`: `111,967`
- `AndrewDHildreth2021`: `110,967`
- `NataliaJaeger2024`: `107,679`
- `DanielPCaron2025`: `95,581`
- `ZhenlongLi2024`: `31,572`
- `MarylineFalquet2023`: `2,446`

Study-level donor, sample, and tissue coverage:

- `TabulaSapiens2022`: `22` donors, `318` samples, `28` tissues
- `StevenBWells2025`: `20` donors, `252` samples, `10` tissues
- `EMadissoon2020`: `12` donors, `58` samples, `3` tissues
- `CDominguezConde2022`: `6` donors, `58` samples, `11` tissues
- `ShuaiHe2020`: `1` donor, `15` samples, `15` tissues
- `AndrewDHildreth2021`: `12` donors, `20` samples, `1` tissue
- `NataliaJaeger2024`: `19` donors, `19` samples, `4` tissues
- `DanielPCaron2025`: `3` donors, `23` samples, `6` tissues
- `ZhenlongLi2024`: `6` donors, `6` samples, `1` tissue
- `MarylineFalquet2023`: `3` donors, `3` samples, `1` tissue

Assays in the object:

- `RNA`: `54,872` features across `4,255,241` cells
- `ADT`: `269` features across `1,161,166` cells in the assay object

Metadata columns available in `reference.raw.qs` include:

- study and sample provenance
- donor-level covariates
- tissue annotations
- assay chemistry
- dataset role/tier/set annotations
- RNA and ADT QC metrics

## Source studies

The current source dataset is derived from 10 published studies:

1. Wells et al. `Multimodal profiling reveals tissue-directed signatures of human immune cells altered with age`
2. Tabula Sapiens Consortium. `Tabula sapiens reveals transcription factor expression, senescence effects, and sex-specific features in cell types from 28 human organs and tissues`
3. Caron et al. `Multimodal hierarchical classification of CITE-seq data delineates immune cell states across lineages and tissues`
4. Jaeger et al. `Diversity of group 1 innate lymphoid cells in human tissues`
5. Dominguez Conde et al. `Cross-tissue immune cell analysis reveals tissue-specific features in humans`
6. Hildreth et al. `Single-cell sequencing of human white adipose tissue identifies new cell states in health and obesity`
7. He et al. `Single-cell transcriptome profiling of an adult human cell atlas of 15 major organs`
8. Madissoon et al. `scRNA-seq assessment of the human lung, spleen, and esophagus tissue stability after cold preservation`
9. Li et al. `Therapeutic application of human type 2 innate lymphoid cells via induction of granzyme B-mediated tumor cell death`
10. Falquet et al. `Dynamic single-cell regulomes characterize human peripheral blood innate lymphoid cell subpopulations`

## What the upload script does

`scripts/r/upload_data.R` performs the following steps:

1. Reads configuration from environment variables instead of hard-coding credentials.
2. Validates that the upload files exist.
3. Validates that at least one creator is provided.
4. Creates a new Zenodo draft deposition or updates an existing one.
5. Uploads `reference.raw.qs` and `bpcells.tar.gz`.
6. Leaves the record as a draft by default.
7. Only publishes when both `ZENODO_PUBLISH=true` and `--execute` are used.

## R package dependencies

The upload workflow depends on:

- `zen4R`
- `jsonlite`
- `glue`

Install them if needed:

```r
install.packages(c("zen4R", "jsonlite", "glue"))
```

## Configuration

Use the example environment file:

- `scripts/r/upload_data.env.example`

For example:

```bash
cp scripts/r/upload_data.env.example .env.zenodo
```

Then edit `.env.zenodo` with your actual values.

### Required fields

These fields must be set:

- `ZENODO_TOKEN` or `ZENODO_SANDBOX_TOKEN`
- `ZENODO_CREATORS_JSON`

Example:

```bash
export ZENODO_CREATORS_JSON='[
  {
    "firstname": "Songqi",
    "lastname": "Duan",
    "affiliations": ["City of Hope"],
    "orcid": "0000-0002-0822-5883"
  }
]'
```

### Common fields

- `ZENODO_USE_SANDBOX=false`
- `ZENODO_DEPOSITION_ID=`
- `ZENODO_PUBLISH=false`
- `ZENODO_OVERWRITE_FILES=true`
- `ZENODO_TITLE`
- `ZENODO_VERSION`
- `ZENODO_PUBLICATION_DATE`
- `ZENODO_LICENSE`
- `ZENODO_RELATED_IDENTIFIERS_JSON`
- `ZENODO_KEYWORDS_JSON`
- `ZENODO_DESCRIPTION`
- `ZENODO_NOTES`

## Recommended workflow

### 1. Load environment variables

```bash
set -a
source .env.zenodo
set +a
```

### 2. Run a dry run

This validates files and metadata and prints the upload plan:

```bash
Rscript scripts/r/upload_data.R
```

### 3. Create or update the draft deposition

```bash
Rscript scripts/r/upload_data.R --execute
```

### 4. Publish after final verification

```bash
export ZENODO_PUBLISH=true
Rscript scripts/r/upload_data.R --execute
```

Zenodo publication is irreversible in the same way a draft is editable, so the draft should be checked carefully before publication.

## Typical usage scenarios

### First upload to production Zenodo

```bash
export ZENODO_USE_SANDBOX=false
export ZENODO_TOKEN='your-production-token'
export ZENODO_CREATORS_JSON='[{"firstname":"Songqi","lastname":"Duan","affiliations":["City of Hope"],"orcid":"0000-0002-0822-5883"}]'

Rscript scripts/r/upload_data.R
Rscript scripts/r/upload_data.R --execute
```

### First upload to Zenodo sandbox

```bash
export ZENODO_USE_SANDBOX=true
export ZENODO_SANDBOX_TOKEN='your-sandbox-token'
export ZENODO_CREATORS_JSON='[{"firstname":"Songqi","lastname":"Duan","affiliations":["City of Hope"],"orcid":"0000-0002-0822-5883"}]'

Rscript scripts/r/upload_data.R
Rscript scripts/r/upload_data.R --execute
```

### Update an existing draft

```bash
export ZENODO_DEPOSITION_ID='your-existing-draft-id'
Rscript scripts/r/upload_data.R --execute
```

## Related identifiers

The Zenodo record should link back to the project repository:

```bash
export ZENODO_RELATED_IDENTIFIERS_JSON='[
  {
    "identifier": "https://github.com/zerostwo/immune-atlas",
    "scheme": "url",
    "relation_type": "isdocumentedby",
    "resource_type": "software"
  }
]'
```

The record can also include the source studies as related identifiers when DOI values are provided.

## Download and local BPCells rebinding

`reference.raw.qs` contains BPCells-backed layers. After users download the dataset, they must extract `bpcells.tar.gz` and rebind the BPCells directories to the Seurat object locally.

### Suggested download layout

```text
immune-atlas-source/
├── reference.raw.qs
└── bpcells.tar.gz
```

Extract the archive:

```bash
cd immune-atlas-source
tar -xzf bpcells.tar.gz
```

The directory should then look like:

```text
immune-atlas-source/
├── reference.raw.qs
├── bpcells.tar.gz
└── bpcells/
    ├── counts/
    ├── decontaminated_counts/
    └── adt/
```

### Minimal R example

```r
library(Shennong)
library(SeuratObject)
library(BPCells)

ref <- sn_read("reference.raw.qs")

LayerData(ref, assay = "RNA", layer = "counts") <-
  BPCells::open_matrix_dir("bpcells/counts")

LayerData(ref, assay = "RNA", layer = "decontaminated_counts") <-
  BPCells::open_matrix_dir("bpcells/decontaminated_counts")

LayerData(ref, assay = "ADT", layer = "counts") <-
  BPCells::open_matrix_dir("bpcells/adt")

ref
table(ref$study)
table(ref$tissue_level1)
```

### Save a reusable local copy

```r
sn_write(ref, "reference.rebound.qs")
```

### Command-line helper script

The repository also includes:

- `scripts/r/rebind_reference_bpcells.R`

Usage:

```bash
Rscript scripts/r/rebind_reference_bpcells.R \
  reference.raw.qs \
  bpcells \
  reference.rebound.qs
```

Arguments:

- argument 1: path to `reference.raw.qs`
- argument 2: path to the extracted `bpcells/` directory
- argument 3: optional output filename for the rebound object

## Recommended short Zenodo usage note

The following note can be included in the Zenodo description or notes:

```text
After download, extract bpcells.tar.gz and rebind the BPCells-backed layers to the Seurat object in reference.raw.qs. Example in R:

library(Shennong)
library(SeuratObject)
library(BPCells)

ref <- sn_read("reference.raw.qs")
LayerData(ref, assay = "RNA", layer = "counts") <- BPCells::open_matrix_dir("bpcells/counts")
LayerData(ref, assay = "RNA", layer = "decontaminated_counts") <- BPCells::open_matrix_dir("bpcells/decontaminated_counts")
LayerData(ref, assay = "ADT", layer = "counts") <- BPCells::open_matrix_dir("bpcells/adt")
```

## Notes

- Do not commit a real Zenodo token to the repository.
- `bpcells.tar.gz` is approximately `11 GB`, so upload time will be substantial.
- For testing, validate metadata and file behavior in Zenodo sandbox first.
- Once a production record is published, future updates should usually be released as a new version rather than overwriting the original deposition.
