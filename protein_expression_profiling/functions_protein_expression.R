
# function loadPackages

loadPackages<- function(packagesToLoad){
  new.packages <- packagesToLoad[!(packagesToLoad %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0) install.packages(new.packages)
  loadedPackages <- NULL
  for (i in 1 : length(packagesToLoad)) {
    if (!packagesToLoad[i] %in% (.packages())) {
      library(packagesToLoad[i], character.only = TRUE)
      loadedPackages <- c(loadedPackages, packagesToLoad[i])
    }
  }
  return(loadedPackages)
}

packagesToLoad <- c("ggplot2", "Rcpp")


################################################
# function detachPackages

detachPackages <- function(loadedPackages) {
  if  (length(loadedPackages >0)) {
    for (i in 1:length(loadedPackages)){
      toDetach <- paste0("package:", loadedPackages[i])
      detach(toDetach, character.only = TRUE)
    }
  }
}

########################################################

# function mapping

mapping <- function(x){
  ifelse(x <= 0.0001, "***", ifelse(x <= 0.001, "**", ifelse(x <= 0.01, "*", "NS")))
}

########################################################

# function getReplicateDirs


getReplicateDirs <- function(dirs){
  replicate_dirs <- NULL
  for (dir_path in dirs){
    replicate_dir <- list.files(path = dir_path, full.names = TRUE, recursive = TRUE, pattern = "replicate_", include.dirs = TRUE)
    replicate_dir <- replicate_dir[file.info(replicate_dir)$isdir]
    if (length(replicate_dir)==0){
      if (stringr::str_detect(string = dir_path, "replicate_[0-9][0-9]?")){
        replicate_dirs <- c(replicate_dirs, stringr::str_extract(dir_path, paste0("^.*replicate.+?", .Platform$file.sep))) 
      }else{
        replicate_dirs <- c(replicate_dirs, dir_path)
        print("No replicates were found")
      }
    }else{
      replicate_dirs <- c(replicate_dirs, replicate_dir)
    }
  }
  return(replicate_dirs)
}

##########################################################

# Function get_y_to_plot()

##############################################################################################


get_y_to_plot <- function(pool_replicate, normalize, background_substraction){
  if(background_substraction){
    y <- "absolute_value"
  }else{
    if(pool_replicate){
      y <- "rel_value"
    }else{
      if (!normalize) y <- "value"
      if (normalize) y <- "rel_value"
    }
  }
  return(y)
}

#############################################################################################

## Function  getReplicateData

##########################################################

getReplicateData <- function(replicate_dir){
  if (stringr::str_detect(string = replicate_dir, "replicate_[0-9][0-9]?")){
    replicate_no <- stringr::str_extract(replicate_dir, "replicate_[0-9]")
    replicate_no <- as.numeric(stringr::str_replace(replicate_no, "replicate_", ""))
  }else{
    replicate_no <- 0
  }
  
  profile_dirs <- list.files(path = replicate_dir, full.names = TRUE, recursive = TRUE, include.dirs = TRUE, pattern = "plot_profile")
  
  profile_dirs <- profile_dirs[file.info(profile_dirs)$isdir]
  genotype <- ""
  raw_data <- data.frame()
  all_markers <- data.frame()
  
  for (j in 1:length(profile_dirs)) {
    x = unlist(str_split(profile_dirs[j], pattern = .Platform$file.sep))
    genotype[j] <- x[length(x)-1]
    profile_data_list <- getProfile(profile_dirs[j], staining, genotype = genotype[j], file_type = "csv", replicate_no = replicate_no) # get both intensity and marker data from every profile directory
    raw_data <- rbind(raw_data, profile_data_list[[1]])
    
    all_markers <- rbind(all_markers, profile_data_list[[2]])
  }
  
  #For every genotype, cd-wise mean values are calculated
  raw_data %>% 
    group_by(genotype, cd) %>%
    summarise(
      value = mean(value)
    ) %>% 
    ungroup() -> mean_gen_cd  
  
  # For Control genotype, peak mean value is calculated
  mean_gen_cd %>% 
  filter(genotype == control_genotype) %>% 
  summarise(
    value = max(value)
  ) -> max_staining_intensity
  
  # Normalized intensity value for each row
  raw_data$rel_value <- raw_data$value * 100 / max_staining_intensity$value

  background_ctrl_genotype <- genotype[grep(pattern = paste0(".*", staining, ".*"), x = genotype, invert = TRUE)]
  expt_genotype <- genotype[grep(pattern = paste0(".*", staining, ".*"), x = genotype)]

  raw_data %>% 
    dplyr::filter(raw_data$genotype %in% background_ctrl_genotype) %>% 
    group_by(cd) %>% 
    summarise(
      rel_value = mean(rel_value)
    ) %>% ungroup() -> background_ctrl_mean_cd
  

  raw_data$absolute_value <- raw_data$rel_value - background_ctrl_mean_cd$rel_value
  raw_data <- normalizeValues(raw_data, control_genotype)
  if (background_substraction == TRUE){
    raw_data <- dplyr::filter(raw_data, raw_data$genotype %in% expt_genotype)
  }
  all_markers$marker1 <- as.numeric(all_markers$marker1)
  all_markers$marker2 <- as.numeric(all_markers$marker2)
  raw_data_list <- list(raw_data, all_markers)
  return(raw_data_list)
}

#####################################################################################

## Function getProfile

######################################################################################

getProfile <- function(dir, staining, genotype, file_type = "csv", replicate_no = 0) {
  file_type <- paste0("*.", file_type)
  genotype_markers <- data.frame() #matrix(nrow = 0, ncol = 4))
  
  staining_dir <- list.files(path = dir, pattern = staining, full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
  
  df = data.frame(gonad = character(0), genotype = character(0), cd = numeric(0), cd_type = numeric(0), value = numeric(0))
  
  file_list = list.files(path = staining_dir, full.names = TRUE, pattern = file_type)

  for (k in 1:length(file_list)) {
    temp_data <- as.data.frame(suppressMessages(read_csv(file = file_list[k], col_names = TRUE, progress = F))) # from readr package
    gonad <- paste0(str_extract(basename(file_list[k]), "^[0-9-]*"), ".tiff")
    
    if (max(temp_data$cd) < cell_diameter) {
      next
    }else{
      markers_data <- getMarkersList(temp_data, genotype, gonad)
      genotype_markers <- rbind(genotype_markers, markers_data)
      temp_data <- temp_data[temp_data$cd <= cell_diameter, ]
    }
    #print(sapply(temp_data, class))
    #temp_data$value <- as.numeric(temp_data$value)
    temp_data %>% 
      ungroup() %>%
      group_by(cd) %>% 
      group_by(cd_type, .add = TRUE) %>%
      summarise(
        value = mean(value),
        n = n()
      ) %>%  ungroup() ->  temp_data
    temp_df <- data.frame(gonad = gonad, genotype = genotype, replicate_no = replicate_no, cd = temp_data$cd, cd_type = temp_data$cd_type, n = temp_data$n,  value = temp_data$value, stringsAsFactors = FALSE)
    df <- rbind(df, temp_df, stringsAsFactors = FALSE)
  }
  colnames(df) <- c("gonad", "genotype", "replicate", "cd", "cd_type", "n", "value")
  colnames(genotype_markers) <- c("genotype", "gonad", "marker1", "marker2")
  
  replicate_data_list <- list(df, genotype_markers)
  return(replicate_data_list)
}

##############################################################################################

normalizeValues <- function(df, control_genotype) {
  df %>% 
    group_by(genotype) %>% 
    group_by(cd, .add = TRUE) %>%
    summarise(
      rel_value = mean(rel_value),
      absolute_value = mean(absolute_value)
    ) %>% ungroup() -> mean_gen_cd
  
  max_abs_intensity <- max(mean_gen_cd$absolute_value[mean_gen_cd$genotype == control_genotype])
  max_rel_intensity <- max(mean_gen_cd$rel_value[mean_gen_cd$genotype == control_genotype])
  
  df$absolute_value <- df$absolute_value * 100 / max_abs_intensity
  df$rel_value <- df$rel_value * 100 / max_rel_intensity
  return(df)
}

#####################################################################################
# Function getMarkersList()

###############################################################################
#,X,Y,Ch,Slice,Counter,Count

getMarkersList <- function(markers_data, genotype, gonad){
  markers_data %>% 
    ungroup() %>% 
    group_by(cd) %>% 
     group_by(cd_type, .add = TRUE) %>%
    # group_by(genotype, .add = TRUE) %>%
    summarise(
      n = n()
    ) -> rows
  rows %>% 
    ungroup() %>% 
    group_by(cd_type) %>% 
    summarise(
      n = n()
    ) -> rows
    genotype_markers <- c(genotype, gonad, rows$n[rows$cd_type == 1], ifelse(max(rows$cd_type) > 2, rows$n[rows$cd_type == 2], NA))
  return(genotype_markers)
}

#############################################################################################

# Function plot_data()

##############################################################################################

plotData <- function(all_replicate_data, all_replicate_rowcount, include_row_data = FALSE, y = "value", show_significance = TRUE, plot_genotype = NULL, replicate_no = NULL) {

  if (include_row_data) {
    rowcount <- getRowcount(all_replicate_rowcount, all_replicate_data, y)
  }
  
  if (!is.null(plot_genotype)){
    all_replicate_data <- all_replicate_data[all_replicate_data$genotype %in% plot_genotype,]
  }
  
  all_genotype <- as.character(unique(all_replicate_data$genotype))

  #colors <- c("navajowhite4", "cornflowerblue", "darkred",  "orangered1")  

    
  mean_list_plot = ""
  for (j in 1:length(all_genotype)){
    mean_list_plot <- paste0(mean_list_plot, "geom_mean(data = all_replicate_data[all_replicate_data$genotype == all_genotype[", j,"], ], color = '", colors[j], "', y_nudge = 0.05* max(all_replicate_data[,y]))")

    mean_list_plot <- ifelse(j < length(all_genotype), paste0(mean_list_plot, " + "), mean_list_plot )
  }
  
  if (length(all_genotype) > 1) {
    comparisons <-  character(0)
    comparisons <- getComparisons(all_genotype, comparisons)
    signif_data <- getSignifData(comparisons, cell_diameter, all_genotype, all_replicate_data)

    maxY <- max(all_replicate_data[,y])
    scaling_factor <- maxY * 0.1
    maxY <- ceiling(maxY/scaling_factor)*scaling_factor + scaling_factor
    set.seed(2032)
    signif_colors <- sample(colors, ncol(signif_data)-1, replace = F)
    
    signif_list_plot = ""
    for (j in 1:length(comparisons)){
      signif_list_plot <- paste0(signif_list_plot, "geom_signif(annotations = mapping(signif_data[,", toString(j+1), "]), y_position = ", toString(maxY + 1.5*scaling_factor*(j-1)), ", xmin=c(1:", toString(cell_diameter), "), xmax=c(1:", toString(cell_diameter), "), color = '", signif_colors[j], "', tip_length = 0, textsize = 5)")
      signif_list_plot <- ifelse(j<length(comparisons), paste0(signif_list_plot, " + "), signif_list_plot )
    }
    
    annotation <- getAnnotation(comparisons, all_genotype, plot_points, cell_diameter, maxY, scaling_factor, signif_colors, colors)
  }
  
  for (genotype in all_genotype){
    if (genotype != staining  & str_detect(string = genotype, pattern = staining)){
      p_title <-  ifelse(!exists("p_title"), str_remove(genotype, staining), paste0(p_title, " & ", str_remove(genotype, staining)))
    }
  }
  if (!exists("p_title")) p_title <- staining
  
  replicate_no <- ifelse(!is.null(replicate_no), paste0(" replicate # ", replicate_no), "")
  p_title <- paste0(substr(staining, 1, str_locate(staining, "_")[1]-1), " staining in ", p_title, ifelse(replicate_no != "", " - ", ""), replicate_no)
  p_title <- str_replace_all(p_title, c("g1_g2" = "gld-2_gld-1", "_" = " "))
  
  
  n_numbers <- data.frame(genotype=all_genotype, n=0)
  for (j in 1:length(all_genotype)){
    n_numbers$genotype[j] <- all_genotype[j]
    n_numbers$n[j] <- nrow(all_replicate_data[all_replicate_data$genotype==all_genotype[j],])/cell_diameter
  }
  p <- "ggplot(all_replicate_data, aes(x = as.factor(cd), y = get(y), color = genotype))"
  
  if (plot_points) {
    p <- paste0(p, " + geom_jitter(width = 0.2, size = 2, alpha = 0.2, aes(shape = genotype))")
  }
  
  p <- paste0(p, " + ", mean_list_plot)
  if (exists("comparisons") && length(comparisons) > 0 && show_significance) {
    p <- paste0(p, " + ", signif_list_plot, " + ", annotation)
  }
  p <- eval(parse(text=p))
  labels <- paste0(str_replace_all(n_numbers$genotype, c("g1_g2" = "gld-2_gld-1", "_" = " ")), " (", n_numbers$n, ")")
  labels <- ifelse(include_row_data, c(labels, rowcount$variable), labels)  
  
  if (include_row_data){
    p = p + scale_color_manual(values = colors, name="Genotype", labels = c(paste0(str_replace_all(n_numbers$genotype, c("_" = " ")), " (", n_numbers$n, ")"), as.character(unique(rowcount$variable))))
  } else {
    p = p + scale_color_manual(values = colors, name="Genotype", labels = paste0(str_replace_all(n_numbers$genotype, c("_" = " ")), " (", n_numbers$n, ")"))
  }

  p = p + theme_bw(base_size = 9, base_family = "arial") +
    theme(panel.grid = element_blank()) +
    expand_limits(y=0) +

    ##  scale_y_continuous(limits = c(min(all_replicate_data$absolute_value)-5,maxY+scaling_factor*(length(comparisons)+0.5)), expand = expand_scale()) +
    labs(title = p_title, x = "Distance from distal end (gcd)", y = paste0(substr(staining[1], 1, str_locate(staining[1], "_")[1]-1), " levels")) +
    theme(plot.title = element_text(face = "italic", hjust = 0.5))
  
  if (include_row_data) {
    p = p + geom_point(data = rowcount, aes(value, y, color = variable), size = 7, color = colors[1]) # 
  }
  p <- p + scale_x_discrete(breaks = c(1, seq(5, cell_diameter, by = 5)))

  print(p)
  return(p)
}


###############################################################################

# Function getRowcount()
# for prom-1_HA, start of protein and progenitor zone size
##########################################################################################

getRowcount <- function(all_replicate_rowcount, all_replicate_data, y){

  
  all_replicate_rowcount %>% 
    group_by(genotype) %>% 
    summarise(
      marker1 = mean(marker1, na.rm = TRUE),
      marker2 = mean(marker2, na.rm = TRUE)
    ) -> rowcount

  rowcount[is.na(rowcount)] <- 0
  rowcount$start_HA <- ifelse(rowcount$marker2 == 0, 0, rowcount$marker1)
  rowcount$end_prog <- rowcount$marker1 + rowcount$marker2
  rowcount <- melt(data = rowcount, c(1:3))
  all_replicate_data %>%
    group_by(genotype) %>%
    group_by(cd, .add = TRUE) %>%
    summarise(
      value = mean(value),
      rel_value = mean(rel_value),
      absolute_value = mean(absolute_value)
    ) -> mean_gen_cd
  
  for (i in 1:nrow(rowcount)){
    if (rowcount$value[i] == 0){
      rowcount$y[i] <- 0
    }else{
      x <- c(floor(rowcount$value[i]), ceiling(rowcount$value[i]), rowcount$value[i] - floor(rowcount$value[i]))
      col_num <- which( colnames(mean_gen_cd) == y )
      z <- unlist(c(mean_gen_cd[mean_gen_cd$genotype == rowcount$genotype[i] & mean_gen_cd$cd == x[1], col_num],
                    mean_gen_cd[mean_gen_cd$genotype == rowcount$genotype[i] & mean_gen_cd$cd == x[2], col_num]
      ))
      z = z[1] + (z[2] - z[1]) * x[3]
      rowcount$y[i] <- z
      
    }
  }
  rowcount <- rowcount[!(rowcount$y==0),]
  return(rowcount)
  
}


#############################################################################################

# Function getComparisions

#############################################################################################

getComparisons <- function(all_genotype, compare_list){
  while(TRUE){
    for (j in 1:(length(all_genotype)-1)){
      if ((j+1) < length(all_genotype)) {
        for (k in (j+1):length(all_genotype)) {
          compare_list <- c(compare_list, paste0(as.character(letters[j]), "_", as.character(letters[k])))
        }
      }else{
        compare_list <- c(compare_list, paste0(as.character(letters[j]), "_", as.character(letters[j+1])))
      }
    }
    return(compare_list)
    break
  }
}

#############################################################################################

##Function getSignifData()

#############################################################################################

getSignifData <- function(comparisons, cell_diameter, all_genotype, all_replicate_data){
  signif_data <- data.frame(matrix(ncol = length(comparisons) + 1, nrow = 0))
  colnames(signif_data) <- c("cd", comparisons)
  for (j in 1:cell_diameter) {
    for (k in 1:length(all_genotype)){
      assign(letters[k], all_replicate_data$absolute_value[all_replicate_data$genotype == all_genotype[k] & all_replicate_data$cd == j])
    }
    cd_data <-  j
    for (k in 1:length(comparisons)){
      x <- unlist(strsplit(comparisons[k], "_"))
      cd_data <- c(cd_data, t.test(x = get(x[1]), y = get(x[2]))$p.value)
    }
    signif_data <- rbind(signif_data, cd_data)
  }
  colnames(signif_data) <- c("cd", comparisons)
  signif_data$cd <- as.integer(signif_data$cd)
  return(signif_data)
}



#######################################################################################

# Function getAnnotation()

###########################################################################################

getAnnotation <- function(comparisons, all_genotype, plot_points, cell_diameter, maxY, scaling_factor, signif_colors, colors){
  annotation <- ""
  for (j in 1:length(comparisons)){
    x <- unlist(strsplit(comparisons[j], "_"))
    toCompare <- paste0(all_genotype[which(letters == x[1])], " vs ", all_genotype[which(letters == x[2])])
    annotation <- paste0(annotation, "annotate('text', x = 1, y = ", toString(maxY+scaling_factor+1.5*scaling_factor*(j-1)), ", label = '", toCompare, "', hjust = 'left', color = '", signif_colors[j],"', fontface = 'italic' )")
    annotation <- ifelse(j<length(comparisons), paste0(annotation, " + "), annotation)
  }
  
  if (!plot_points) {
    for (j in 1:length(all_genotype)){
      annotation <- paste0(annotation, " + annotate('text', x = ", toString(cell_diameter-6), ", y = ", toString(maxY - j * scaling_factor), ", label = '", all_genotype[j], "', hjust = 'left', color = '", colors[j],"', fontface = 'italic' )")
    }
  }
  annotation <- str_replace_all(annotation, c("g1_g2" = "gld-2_gld-1", "_" = " "))
  return(annotation)
}

############################################################################################

## Function geom_mean()

############################################################################################

geom_mean <- function(data = NULL, size = 8, width = 0.2, color = "blue", x_nudge = 0, y_nudge = 0) {
  mean_list <- list(
    #stat_summary(data=data, fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = 0.3, linetype = "solid", size = 2, color = color),
    stat_summary(data=data, fun=mean, color = color, geom="line", aes(group = 1), size = 1.25),
    stat_summary(data=data, fun=mean, color = color, geom="point", size=2, fill = color, aes(group = 1))
  )
  
  if (show_label){
    mean_list <- c(mean_list, stat_summary(data=data, fun=mean, color = "blue", geom="text", aes(group = 1, label=sprintf("%1.0f", ..y..)), size=5, show.legend=FALSE, position_nudge(x = 0, y= y_nudge)))
  }
  
  if (show_error_bar){
    mean_list <- c(mean_list, stat_summary(data=data, fun.data = "data_summary", geom = "errorbar", width = 0.2, color = color, size = 0.2))
  }
  
  return(mean_list)
}
#############################################################################################

############################################################################################

getNudge <- function(x){
  nudge <- 1.1 * data_summary(x)[3]
  return(nudge)
}

##############################################################################################

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

