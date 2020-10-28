rm(list=ls())
#rm(list=ls()[grep("dir_paths?[0-9]?", ls())])
color_list <- c("orange", "lawngreen", "blue", "indianred1", "black", "purple", "royalblue", "olive", "cyan", "darkturquoise", "aquamarine4", "darkgreen", "navajowhite4", "cornflowerblue", "magenta", "violet", "sienna", "red", "darkorange", "slategray", "mediumslateblue", "lightseagreen", "deeppink", "rosybrown", "darkgoldenrod", "dodgerblue", "olivedrab", "darkviolet", "forestgreen", "deepskyblue")
colors <- color_list
cell_diameter <- 35
staining <- "PROM-1_HA"
control_genotype <- staining
#control_genotype <- paste0(staining, "_LAG-1_HA")
#plot_genotype <- c("N2", staining, paste0(staining, "_bn18"))
plot_genotype <- NULL
#plot_genotype <- c("N2", "gld-2_gld-1_LAG-1_HA", "gld-2_gld-1_glp-1_LAG-1_HA", "gld-2_gld-1_lst-1_sygl-1_LAG-1_HA") #, "gld-1_LAG-1_HA", "gld-2_LAG-1_HA")  ## can specify the genotype to plot
include_row_data <- T
show_significance <- T
show_error_bar <- T
show_label <- F
plot_points <- T
normalize <- T
background_substraction <- F
pool_replicate <- T

dir_path1 <- "~/Dropbox/Acad/computational_work/R/protein_expression_analysis/PROM-1_HA_analysis/09.12.2020_TEMP/01.07.2020_bn18_PROM-1_HA_WAPL_nop1/"
#dir_path1 <- "/run/media/areeba/passport/image_analysis/PROM-1_HA/2020/done_vincent/01.21.2020_oz231_PROM-1_HA_WAPL_nop1/"

options(dplyr.summarise.inform=F)  # reduces some of the summary outputs by dplyr

if (normalize == FALSE) {
  pool_replicate <- FALSE
  background_substraction <- FALSE
}

script_dir <- getSrcDirectory(function(x) {x})  # to get path of the current script

functions_file  <- list.files(script_dir, pattern = ".*function.*", recursive = FALSE, full.names = TRUE)

source(functions_file[1])

packagesToLoad = c("stringr", "ggplot2", "dplyr", "ggsignif", "reshape2", "tools", "readr")
loadedPackages <- loadPackages(packagesToLoad)

dir_paths <- ls()[grep("dir_path[0-9]?", ls())]

dirs <- NULL
for (dir_path in dir_paths) dirs <- c(dirs, get(dir_path))

replicate_dirs <- getReplicateDirs(dirs)
print(replicate_dirs)

all_replicate_data <- data.frame()

all_replicate_rowcount <- data.frame()

for (i in 1:length(replicate_dirs)){
  
  y = get_y_to_plot(pool_replicate, normalize, background_substraction) # Determines what to plot raw value or narmalized etc...
  
  replicate_data_list <- getReplicateData(replicate_dir = replicate_dirs[i])
  replicate_data <- replicate_data_list[[1]]
  #stop("Stops Now")
  replicate_rowcount <- replicate_data_list[[2]]
  if (pool_replicate == FALSE){
    
    all_replicate_data <- replicate_data
    if (stringr::str_detect(string = replicate_dirs[i], "replicate_[0-9][0-9]?")){
      replicate_no <- stringr::str_extract(replicate_dirs[i], "replicate_[0-9]*")
      replicate_no <- stringr::str_replace(replicate_no, "replicate_", "")
    }

    all_replicate_rowcount <- replicate_rowcount
    p <- plotData(all_replicate_data, all_replicate_rowcount, include_row_data, y = y, show_significance, plot_genotype = plot_genotype, replicate_no = replicate_no)

#if pool_replicates  
  }else{
    all_replicate_data <- rbind(all_replicate_data, replicate_data)
    all_replicate_rowcount <- rbind(all_replicate_rowcount, replicate_rowcount)
  
    if (i == length(replicate_dirs)){
      all_replicate_data <- normalizeValues(all_replicate_data, control_genotype)
      p <- plotData(all_replicate_data, all_replicate_rowcount, include_row_data, y = y, show_significance, plot_genotype = plot_genotype)
    }
  }
}

# writeData <- function() {
#   all_replicate_data <- normalizeValues(all_replicate_data, control_genotype)
#   all_replicate_data %>%
# #    dplyr::filter(staining == all_stainings[1]) %>% 
#     group_by(genotype) %>%
# #    group_by(staining, add = TRUE) %>%
#     group_by(cd, add = TRUE) %>%
#     summarise(
#       value = round(mean(value), 0),
#       rel_value = round(mean(rel_value), 0),
#       absolute_value = round(mean(absolute_value), 0)
#     ) -> mean_gen_cd
#   
#   write.csv(x = mean_gen_cd, file = paste0(dir_path, "cd_values.csv"))
# }
# 
# writeData()

detachPackages(loadedPackages)
rm(list=ls())
