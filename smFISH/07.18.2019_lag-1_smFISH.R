workspace <- ls()
wd <- getwd()
cell_diameter <- 25
smFISH <- "lag-1"
show_significance <- T
show_error_bar <- T
show_label <- F
plot_points <- F
pool_replicate <- T
threshold_pixels <- 100

threshold <- threshold_pixels/2

source("~/Dropbox/Acad/computational_work/R/all_functions/07.18.2019_smFISH_functions.R")

library(statsFunctions)

packagesToLoad = c("stringr", "reshape2", "ggplot2", "RSQLite", "dplyr", "xlsx", "ggsignif", "svDialogs", "XML", "LearnGeom", "ggthemes")
loadedPackages <- loadPackages(packagesToLoad)
loadedPackages <- c(loadedPackages, "statsFunctions")

mapping <- function(x){
  ifelse(x <= 0.0001, "***", ifelse(x <= 0.001, "**", ifelse(x <= 0.01, "*", "NS")))
}


dir_path <- "/run/media/areeba/8452B4A052B497FE/2019/06.19.2019_lag-1_smFISH/gld-2_gld-1/replicate_3/"
dir_path <- "/run/media/areeba/8452B4A052B497FE/2019/06.19.2019_lag-1_smFISH/"
dir_path <- "/run/media/fiji/8452B4A052B497FE/2019/smFISH/"
#dir_path <- "/run/media/areeba/8452B4A052B497FE/2019/smFISH/07.14.2019_gld-2_gld-1_lag-1_smFISH_L4/replicate_3/"


replicate_dirs <- list.files(path = dir_path, full.names = TRUE, recursive = TRUE, pattern = "replicate", include.dirs = TRUE)
replicate_dirs <- replicate_dirs[file.info(replicate_dirs)$isdir]
if (length(replicate_dirs)==0){
  replicate_dirs <- dir_path
}

print(replicate_dirs)

all_replicate_data <- data.frame()

# if (length(replicate_dirs) > 1){
#   plots <- NULL
# }

for (i in 1:length(replicate_dirs)){
  replicate_data <- getReplicateData(replicate_dir = replicate_dirs[i], cell_diameter)
  replicate_data$foci <- as.numeric(replicate_data$foci)
  if (pool_replicate == FALSE){
    all_replicate_data <- replicate_data
      p <- plotData(all_replicate_data, show_significance)
  }else{
    all_replicate_data <- rbind(all_replicate_data, replicate_data)
    if (i == length(replicate_dirs)){
      
      p <- plotData(all_replicate_data, show_significance)
    }
  }
}


writeData <- function() {
  all_replicate_data <- normalizeValues(all_replicate_data)
  all_replicate_data %>%
    dplyr::filter(staining == all_stainings[1]) %>% 
    group_by(genotype) %>%
    group_by(staining, add = TRUE) %>%
    group_by(cd, add = TRUE) %>%
    summarise(
      value = round(mean(value), 0),
      rel_value = round(mean(rel_value), 0),
      absolute_value = round(mean(absolute_value), 0)
    ) -> mean_gen_cd
  
  write.csv(x = abc, file = paste0(dir_path, "cd_values.csv"), sep = ",")
}

#writeData()



detachPackages(loadedPackages)
# setwd(wd)
#rm(list = setdiff(ls(), workspace))
 rm(list=ls())


#Sys.getenv("USER")
