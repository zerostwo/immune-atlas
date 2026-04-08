library(tidyverse)
library(Shennong)
library(glue)

# TapsiKumar2023 ----------------------------------------------------------
# breast
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195665
outdir <- sn_set_path("/home/sduan/data/TapsiKumar2023/fastqs")

# MonikaLitvinukova2020 ---------------------------------------------------
# heart
# https://www.ebi.ac.uk/ena/browser/view/ERP123138
outdir <- sn_set_path("/home/sduan/data/MonikaLitvinukova2020/fastqs")
metadata <- sn_read("/home/sduan/data/MonikaLitvinukova2020/fastqs/ERP123138.tsv")
metadata$submitted_ftp
a1 <- str_split_fixed(metadata$submitted_ftp,";", 3)[,1]
a2 <- str_split_fixed(metadata$submitted_ftp,";", 3)[,2]
a3 <- str_split_fixed(metadata$submitted_ftp,";", 3)[,3]
urls <- c(a1,a2,a3)
glue("axel -n10 {a1}") |> 
  cat(file = glue("{outdir}/download_i.sh"), sep = "\n")
glue("axel -n10 {a2}") |> 
  cat(file = glue("{outdir}/download_r1.sh"), sep = "\n")
glue("axel -n10 {a3}") |> 
  cat(file = glue("{outdir}/download_r2.sh"), sep = "\n")

df <- tibble(
  i1_url = a1,
  r1_url = a2,
  r2_url = a3,
  i1 = glue("{outdir}/{basename(a1)}"),
  r1 = glue("{outdir}/{basename(a2)}"),
  r2 = glue("{outdir}/{basename(a3)}")
)

df <- tibble(
  url = c(a1,a2,a3),
  file = glue("{outdir}/{basename(url)}"),
) |> 
  mutate(group = case_when(
    str_detect(file, "_I1_") ~ "I1",
    str_detect(file, "_R1_") ~ "R1",
    str_detect(file, "_R2_") ~ "R2",
    TRUE ~ "Unknown"
  ))

filtered_df <- df |> 
  filter(!file.exists(file))

i1_urls <- filtered_df |> 
  filter(group =="I1") |> pull(url)
glue("axel -n10 {i1_urls}") |> 
  cat(file = glue("{outdir}/download_i.sh"), sep = "\n")

r1_urls <- filtered_df |> 
  filter(group =="R1") |> pull(url)
glue("wget {r1_urls}") |> 
  cat(file = glue("{outdir}/download_r1.sh"), sep = "\n")

list.files(outdir, pattern = "fastq", full.names = TRUE) |>
  fs::file_info() -> fi
fi |>
  filter(size ==0) |> pull(1) |> file.remove()

# DrakeWinslowWilliams2021 ------------------------------------------------
# Oral
# GSE164241
outdir <- sn_set_path("/home/sduan/data/DrakeWinslowWilliams2021/fastqs")


