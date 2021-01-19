library("sleuth")
dir <- "/home/manager/course_data/WTAC_NGS_Chile2021_projectRNAseq_1b/data/"
s2c <- read.table(paste0(dir, "hiseq_info.txt"), header = TRUE, stringsAsFactors=FALSE)

sample_id <- c("WT1", "WT2","WT3", "KO1", "KO2", "KO3")
kal_dirs <- sapply(sample_id, function(id) file.path(".", id))
s2c <- dplyr::select(s2c, sample = sample, condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

t2g <- read.table(paste0(dir, "Pb.names.txt"), header=T, sep="\t")

so <- sleuth_prep(s2c, ~ condition, extra_bootstrap_summary = TRUE, target_mapping = t2g)
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
#so <- sleuth_lrt(so, 'reduced', 'full')
so <- sleuth_wt(so, 'conditionWT', 'full')

results_table <- sleuth_results(so, 'conditionWT', test_type = 'wt')

write.table(results_table, file=paste0(dir, "kallisto.results"), quote=FALSE, sep="\t", row.names=FALSE)

sleuth_live(so)
