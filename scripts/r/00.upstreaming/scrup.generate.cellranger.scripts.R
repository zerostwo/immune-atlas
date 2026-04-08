library(tidyverse)

r1 <- list.files("/home/sduan/data/scrup/fastqs",
                 pattern = "R1.fastq.gz$",
                 full.names = TRUE)

adata <- tibble(
  sample = str_split_fixed(basename(r1), "_", 2)[,1],
  r1 = r1,
  r2 = str_replace_all(r1, "_R1", "_R2")
)
table(file.exists(adata$r1))
table(file.exists(adata$r2))

outdir <- "/home/sduan/data/scrup/fastqs"

adata %>% 
  rowwise() %>%
  mutate(
    out_r1 = file.path(outdir, paste0(sample, "_S1_L001_R1_001.fastq.gz")),
    out_r2 = file.path(outdir, paste0(sample, "_S1_L001_R2_001.fastq.gz"))
  ) %>%
  ungroup() %>%
  rowwise() %>%
  do({
    file.symlink(.$r1, .$out_r1)
    file.symlink(.$r2, .$out_r2)
    tibble(sample = .$sample, r1_link = .$out_r1, r2_link = .$out_r2)
  })

library(glue)
cellranger <- "/home/sduan/software/cellranger/9.0.1/bin/cellranger"
reference <- "/home/sduan/reference_genomes/10xgenomics/GRCh38"
cmds <- glue("{cellranger} count --id={adata$sample} --transcriptome={reference} --fastqs={outdir} --sample={adata$sample} --localcores=32 --localmem=128 --output-dir=/home/sduan/data/scrup/alignments/cellranger/{adata$sample} --create-bam true --nosecondary")
cat(cmds, sep = "\n", file = "scripts/shell/scrup.cellranger_count.sh")
n_batches <- 20
cmd_tbl <- tibble(cmd = cmds) %>%
  mutate(batch = ceiling(row_number() / (n() / n_batches)))
outdir <- "scripts/shell/bash_batches/scrup"
dir.create(outdir, recursive = TRUE)

cmd_tbl %>%
  group_split(batch) %>%
  walk2(1:length(.), function(df, i) {
    file <- file.path(outdir, glue("batch_{i}.sh"))
    writeLines(c(
      "#!/bin/bash",
      "#SBATCH --time=14-00:00:00",
      glue("#SBATCH --job-name=scrup.batch_{i}"),
      "#SBATCH --cpus-per-task=32",
      "#SBATCH --mem=128G",
      glue("#SBATCH --output=/home/sduan/data/scrup/logs/cellranger/batch_{i}.%j.log"),
      "#SBATCH --mail-type=BEGIN,END,FAIL",
      "#SBATCH --mail-user=sduan@coh.org",
      "",
      df$cmd
    ), file)
    print(glue("bash /home/sduan/TabulaSapiens/{file}"))
    Sys.chmod(file, mode = "0755") 
    print(file)
  })



