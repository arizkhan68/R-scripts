## Function  getReplicateData

##########################################################

getReplicateData <- function(replicate_dir, cell_diameter){
  
  markers_dirs <- list.files(path = replicate_dir, full.names = TRUE, recursive = TRUE, include.dirs = TRUE, pattern = "markers")
  genotype <- ""
  raw_data <- data.frame()
  
  for (j in 1:length(markers_dirs)) {
    cd_dir <- paste0(dirname(markers_dirs[j]), "/cd/")
    x = unlist(str_split(markers_dirs[j], pattern = "/"))
    genotype[j] <- x[length(x)-2]
    tempdata <- getMarkers(markers_dirs[j], cd_dir, cell_diameter, file_type = "csv")
    
    tempdata$genotype <- as.factor(genotype[j])
    raw_data <- rbind(raw_data, tempdata)
    
  }

  return(raw_data)
}
#####################################################################################

## Function getMarkers

######################################################################################

getMarkers <- function(dir, cd_dir, cell_diameter, file_type = "csv") {
  df = data.frame(gonad = numeric(0), cd = numeric(0), foci = numeric(0) )
  file_list = list.files(path = dir, full.names = TRUE, pattern = paste0("*.", file_type))
  for (k in 1:length(file_list)) {
    gonad <- basename(file_list[k])
    gonad <- gsub(pattern = paste0(".", file_type), replacement = "", gonad)
    print(gonad)
    temp_data <- read.csv2(file_list[k], header = TRUE, sep = ",", stringsAsFactors = FALSE)
    temp_data <- temp_data[,c("BX", "BY")]
    
    colnames(temp_data) <- c("x", "y")
    #temp_data$x <- as.integer(temp_data$x)
    #temp_data$y <- as.integer(temp_data$y)
    cd_file <- basename(file_list[k])  # str_remove(basename(file_list[k]), ".csv")
    cd_file <- paste0(cd_dir, cd_file)
    cd_data <- read.csv(cd_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
    if (ncol(cd_data) == 9){
      colnames(cd_data) <- c("cd", "area", "mean", "min", "max", "x", "y", "ch", "slice")
    }else{
      colnames(cd_data) <- c("cd", "area", "mean", "x", "y", "ch", "slice")
    }
    
    
    cd_data <- cd_data[c(1:(1 + cell_diameter)), c("cd", "x", "y")]
    cd_data$cd <- as.integer(cd_data$cd)
    cd_data$x <- as.integer(cd_data$x)
    cd_data$y <- as.integer(cd_data$y)
    
    temp_data$cd <- NA
    for (l in 1:nrow(temp_data)){
      p1 <- c(temp_data[l, "x"], temp_data[l, "y"])
      
      cd_data1 <- cd_data
      cd_data1$distance <- NA
      for (n in 1:nrow(cd_data1)){
        cd_data1$distance[n] = DistancePoints(p1, c(cd_data1$x[n], cd_data1$y[n]))
      }

      cd_data1 <- cd_data1[order(cd_data1$distance),]

      if (cd_data1$distance[1] < threshold){
        if (gonad == "005-006" & l == 337) print(paste(l, gonad))
        temp_data$cd[l] <- cd_data1$cd[1]
      }else if (abs(cd_data1$cd[1] - cd_data1$cd[2]) == 1){
        
        b <- c(cd_data1$x[1], cd_data1$y[1])
        c <- c(cd_data1$x[2], cd_data1$y[2])
        distance <- getDistance(p1,b,c)
        if (distance <= threshold){
          #print(c(cd_data1$distance[1], distance, threshold))
          temp_data$cd[l] <- cd_data1$cd[1]

        }
      }
    }
    
    temp_data <- temp_data[!is.na(temp_data$cd),]
    
    gonad <- basename(file_list[k])
    gonad <- gsub(pattern = paste0(".", file_type), replacement = "", gonad)
    x <- table(temp_data$cd)
    x <- as.data.frame(x)
    print(gonad)

    colnames(x) <- c("cd", "freq")

    for (m in (1 : cell_diameter)) {
      freq = x$freq[x$cd==m]
      if (length(freq) == 0) freq <- 0

      temp_df <- data.frame(gonad = gonad, cd = m, foci = freq)
      df <- rbind(df, temp_df, stringsAsFactors = FALSE)
    }
  }

  colnames(df) <- c("gonad", "cd", "foci")
  df$gonad <- as.factor(df$gonad)
  df$cd <- as.factor(df$cd)
  return(df)
}
###########################################################################################

getDistance <- function(a,b,c){
  line <- CreateLinePoints(b, c)
  p <- ProjectPoint(a, line)
  if (DistancePoints(b,c) > DistancePoints(b,p)){
    d <- DistancePointLine(a, line)
  }else{
    d <- DistancePoints(a, b)
  }
  
  return(d)
}

#############################################################################################

# Function plot_data()

##############################################################################################

plotData <- function(all_replicate_data, show_significance = TRUE, y = "foci") {

 
  all_genotype <- as.character(unique(all_replicate_data$genotype))
  
  n_colors <- length(all_genotype)
  colors <- darkColorPicker(n_colors)
  colors <- c("darkorchid4", "pink")
  colors <- c("lawngreen", "orange", "gray30")
  #colors <- c("gray30", "orange")
  print(colors)
  
  mean_list_plot = ""
  n = 0
  for (j in 1:length(all_genotype)){
      n = n + 1
      mean_list_plot <- paste0(mean_list_plot, "geom_mean(data = all_replicate_data[all_replicate_data$genotype == all_genotype[", j,"], ], color = '", colors[n], "', y_nudge = 0.05* max(all_replicate_data[,y]))")
      mean_list_plot <- ifelse(j < length(all_genotype), paste0(mean_list_plot, " + "), mean_list_plot )
  }
  
  if (length(all_genotype) > 1) {
    comparisons <-  character(0)
    comparisons <- getComparisons(all_genotype, comparisons)
    signif_data <- getSignifData(comparisons, cell_diameter, all_genotype, all_replicate_data)
    maxY <- max(all_replicate_data[,y], na.rm = TRUE)
    scaling_factor <- maxY * 0.1
    maxY <- ceiling(maxY/scaling_factor)*scaling_factor + scaling_factor
    signif_list_plot = ""
    signif_colors <- darkColorPicker(ncol(signif_data)-1)
    
    # signif_colors <- c("mediumorchid",  "midnightblue", "mediumpurple2", "blue4", "paleturquoise3", "chartreuse3")   
    #signif_colors <- "dodgerblue"

    for (j in 1:length(comparisons)){
      signif_list_plot <- paste0(signif_list_plot, "geom_signif(annotations = mapping(signif_data[,", toString(j+1), "]), y_position = ", toString(maxY + 1.5*scaling_factor*(j-1)), ", xmin=c(1:", toString(cell_diameter), "), xmax=c(1:", toString(cell_diameter), "), color = '", signif_colors[j], "', tip_length = 0, textsize = 5)")
      signif_list_plot <- ifelse(j<length(comparisons), paste0(signif_list_plot, " + "), signif_list_plot )
    }
    
    annotation <- getAnnotation(comparisons, all_genotype, plot_points, cell_diameter, maxY, scaling_factor, signif_colors, colors)
  }
  
  for (genotype in all_genotype){
    p_title <-  ifelse(!exists("p_title"), genotype, paste0(p_title, " & ", genotype))
  }
  
  n_numbers <- data.frame(genotype=all_genotype, n=numeric(length(all_genotype)))
  for (j in 1:length(all_genotype)){
    n_numbers$genotype[j] <- all_genotype[j]
    n_numbers$n[j] <- nrow(all_replicate_data[all_replicate_data$genotype==all_genotype[j],])/cell_diameter
  }

  p <- "ggplot(all_replicate_data, aes(x = cd, y = get(y), color = genotype))"
  
  if (plot_points) {
    p <- paste0(p, " + geom_jitter(width = 0.2, size = 2, alpha = 0.4)")
  }
  
  p <- paste0(p, " + ", mean_list_plot)
  if (exists("comparisons") && length(comparisons) > 0 && show_significance) {
    p <- paste0(p, " + ", signif_list_plot, " + ", annotation)
  }
  print(p)
  
  p <- eval(parse(text=p))
  
  labels <- paste0(str_replace_all(n_numbers$genotype, c("g1_g2" = "gld-2_gld-1", "_" = " ")), " (", n_numbers$n, ")")
  
  

  p = p + scale_color_manual(values = colors, name="Genotype", labels = paste0(str_replace_all(n_numbers$genotype, c("g1_g2" = "gld-2_gld-1", "_" = " ")), " (", n_numbers$n, ")"))
  
  
  p = p + theme_few() +
    expand_limits(y=0)  +

    ##  scale_y_continuous(limits = c(min(all_replicate_data$absolute_value)-5,maxY+scaling_factor*(length(comparisons)+0.5)), expand = expand_scale()) +
    labs(title = p_title, x = "Distance from distal end (gcd)", y = "No of foci" )  +
    theme(plot.title = element_text(face = "italic", hjust = 0.5))
  
  p <- p + scale_x_discrete(breaks = c(1, seq(5, cell_diameter, by = 5)))
  
  print(p)
  
  return(p)
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
      #  annotation <- ifelse(j<length(all_genotype), paste0(annotation, " + "), annotation)
    }
  }
  annotation <- str_replace_all(annotation, c("g1_g2" = "gld-2_gld-1", "_" = " "))
  return(annotation)
}

############################################################################################


## Function listDirs

#############################################################################################

listDirs <- function(path=".", pattern=NULL, all.dirs=FALSE, full.names=TRUE, ignore.case=FALSE, recursive=FALSE) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path=path, pattern=pattern, all.files=all.dirs,
                    full.names=full.names, recursive=recursive, ignore.case=ignore.case)
  dirs <- all[file.info(all)$isdir]
  # determine whether to return full names or just dir names
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
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
      assign(letters[k], all_replicate_data$foci[all_replicate_data$genotype == all_genotype[k] & all_replicate_data$cd == j])
    }
    cd_data <-  j
    for (k in 1:length(comparisons)){
      x <- unlist(strsplit(comparisons[k], "_"))
      cd_data <- c(cd_data, t.test(x = get(x[1]), y = get(x[2]))$p.value)
    }
    signif_data <- rbind(signif_data, cd_data)
  }
  
  colnames(signif_data) <- c("cd", comparisons)
  
  signif_data$cd <- as.factor(signif_data$cd)
  return(signif_data)
}



###########################################################################
#theme theme_bluewhite

########################################################################


theme_bluewhite <- function (base_size = 11, base_family = "") {
  theme_bw() %+replace% 
    theme(
      panel.grid.major  = element_line(color = "white"),
      panel.background = element_rect(fill = "lightblue"),
      panel.border = element_rect(color = "lightblue", fill = NA),
      axis.line = element_line(color = "lightblue"),
      axis.ticks = element_line(color = "lightblue"),
      axis.text = element_text(color = "steelblue")
    )
}



#########################################################################################

## Function geom_mean()

############################################################################################

geom_mean <- function(data = NULL, size = 8, width = 0.2, color = "blue", x_nudge = 0, y_nudge = 0) {
  mean_list <- list(
    #stat_summary(data=data, fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = 0.3, linetype = "solid", size = 2, color = color),
    stat_summary(data=data, fun.y=mean, color = color, geom="line", aes(group = 1), size = 1.25),
    stat_summary(data=data, fun.y=mean, color = color, geom="point", size=4, fill = color, aes(group = 1))
  )
  
  if (show_label){
    mean_list <- c(mean_list, stat_summary(data=data, fun.y=mean, color = "blue", geom="text", aes(group = 1, label=sprintf("%1.0f", ..y..)), size=5, show.legend=FALSE, position_nudge(x = 0, y= y_nudge)))
  }
  
  if (show_error_bar){
    mean_list <- c(mean_list, stat_summary(data=data, fun.data = "data_summary", geom = "errorbar", width = 0.2, color = color, size = 0.5))
  }
  
  return(mean_list)
}

#############################################################################################

############################################################################################

getNudge <- function(x){
  nudge <- 1.1 * data_summary(x)[3]
  return(nudge)
}



###################################################################################################

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
