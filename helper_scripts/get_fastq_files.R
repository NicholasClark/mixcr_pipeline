
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2) {
    stop("Usage: Rscript get_fastq_files.R <path>")
}

get_fastq_files = function(path) {
    #all_fastq = dir(path, pattern = "MVP*_L00[1-8]_R[12]_001.fastq.gz")
    df = tibble(r1_file = dir(path, pattern = ".*_R1_001.fastq.gz")) %>%
        mutate(
            r2_file = gsub("_R1_001.fastq.gz", "_R2_001.fastq.gz", r1_file, fixed = TRUE),
            base_path = path,
            prefix = gsub("_L00[1-8]_R1_001\\.fastq\\.gz", "", r1_file),
            r2_exists = file.exists(file.path(path, r2_file))
        ) %>%
        mutate(
            r1_file_long = file.path(base_path, r1_file),
            r2_file_long = file.path(base_path, r2_file)
        ) %>%
        filter(!grepl("Undetermined", r1_file))
    df_not_found = df %>% filter(!r2_exists)
    if(nrow(df_not_found) > 0) {
        for(i in 1:nrow(df_not_found)) print(paste("File not found:", df_not_found$r2_file[i]))
    }
    df_found = df %>% filter(r2_exists)
    return(df_found)
}

df = get_fastq_files(args[1])
write.csv(df, "fastq_files.csv", row.names = FALSE)