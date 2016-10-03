forest.robu <- function(x, es.lab, study.lab, ...){
  
  if (paste(x$ml[3]) != 1){
    stop("Requires an intercept-only model.")
  }  
  ellipsis        <- lapply(substitute(list(...))[-1L], deparse)
  n_user_cols     <- length(ellipsis) # num. of additional columns
  reg_table       <- x$reg_table    
  N               <- x$N # num. of studies
  M               <- x$M # num. of clusters 
  n_rows          <- M + (2 * N) + 4
  data            <- as.data.frame(x$data) 
  data.full       <- as.data.frame(x$data.full) 
  data$orig.study <- as.factor(x$study_orig_id)
  data            <- data[order(data$orig.study),]            
  data$r.weights  <- data.full$r.weights
  data$effect.size  <- data.full$effect.size
  data$var.eff.size <- data.full$var.eff.size
  data$study.num  <- data.full$study 
  add_col_titles  <- as.list(names(ellipsis)) # user supplied titles
  add_col_values  <- as.list(data[, unlist(ellipsis, use.names = FALSE)]) 
  id_col_title    <- "Studies" 
  id_col_study_values <- unique(data[,study.lab]) 
  id_col_es_values    <- as.character(data[,es.lab]) 
  data$obs_num    <- seq(1, M)
  data$study_num  <- data$study.num 
  data$es_rows    <- as.numeric(data$obs_num + (2 * data$study_num) + 1) 
  data$study_rows <- as.numeric(ave(data$es_rows, data$study_num, FUN = min)- 1)
  es_rows         <- data$es_rows 
  study_rows      <- unique(data$study_rows) 
  total_row       <- max(n_rows)
  title_row       <- min(n_rows)
  data_col_values <- data[, c("r.weights", "effect.size", "var.eff.size")]
  data_col_values <- cbind(data_col_values, es_rows)
  grand.ES        <- reg_table$b.r
  grand.CI.L      <- reg_table$CI.L 
  grand.CI.U      <- reg_table$CI.U
  
  is.numeric_df   <- function(x) all(sapply(x, is.numeric))
  specify_decimal <- function(x, k) format(round(x, k), nsmall = k)
  
  makeTextGrob <- function(values, rows, just = "left", bold = FALSE){  
    if (is.numeric_df(values)) 
      values <- lapply(values, function (x) specify_decimal(x, 3))
    if (!bold){
      t <- lapply(values, function(x) textGrob(paste(x), x = 0, just = just))
    } else {
      t <- lapply(values, function(x) textGrob(paste(x), x = 0, just = just,
                                               gp = gpar(fontface = "bold")))
    }
    return(list(values = t, rows = rows)) 
  }
  
  addTitleToGrob <- function(col, title){
    titleGrob  <- makeTextGrob(values = title, rows = 1, bold = TRUE)
    values     <- c(col$values, titleGrob$values)
    rows       <- c(col$rows, titleGrob$rows)
    return(list(values = values, rows = rows)) 
  }
  
  addGrobToGrob <- function(col1, col2){
    values <- c(col1$values, col2$values)
    rows   <- c(col1$rows, col2$rows)
    return(list(values = values, rows = rows)) 
  } 
  
  makeDataGrob <- function(x){  
    ES      <- x$effect.size
    size    <- x$r.weights / max(x$r.weights)
    CI.U    <- x$effect.size + (1.96 * sqrt(x$var.eff.size))
    CI.L    <- x$effect.size - (1.96 * sqrt(x$var.eff.size))
    type    <- rep("n", M)   
    rows    <- x$es_rows
    return(list(type = type, rows = rows, size = size, CI.L = CI.L, CI.U = CI.U, 
                ES = ES)) 
  }
  
  addSummaryToDataGrob <- function(x){ 
    type  <- c(x$type, "s")
    rows  <- c(x$rows, max(x$rows) + 2)
    size  <- as.numeric(x$size)
    size  <- x$size / max(x$size)
    ES    <- c(x$ES, grand.ES)
    CI.L  <- c(x$CI.L, grand.CI.L)
    CI.U  <- c(x$CI.U, grand.CI.U)
    min   <- floor(as.numeric(min(CI.L)))
    max   <- ceiling(as.numeric(max(CI.U)))
    range <- c(min, max)
    return(list(type = type, rows = rows, size = size, CI.L = CI.L, CI.U = CI.U, 
                ES = ES, min = min, max = max, range = range)) 
  }
  
  if (n_user_cols > 1){
    add_col <- lapply(add_col_values, function(x) makeTextGrob(x, es_rows))
    add_col <- Map(function(x, y) addTitleToGrob(x, y), add_col, add_col_titles)
  }
  
  if (n_user_cols == 1){
    add_col <- makeTextGrob(add_col_values, es_rows)
    add_col <- addTitleToGrob(add_col, add_col_titles)
  }
  
  id_col_study_grob <- makeTextGrob(id_col_study_values, study_rows, bold =TRUE)
  id_col_es_grob    <- makeTextGrob(id_col_es_values, es_rows)
  id_col            <- addGrobToGrob(id_col_study_grob, id_col_es_grob)
  id_col            <- addTitleToGrob(id_col, id_col_title) 
  data_col          <- makeDataGrob(data_col_values)
  data_col          <- addSummaryToDataGrob(data_col)
  
  drawLabelCol <- function(col, j) {
    for (i in 1:length(col$rows)) {
      pushViewport(viewport(layout.pos.row = col$rows[i], layout.pos.col = j))
      grid.draw(col$values[[i]])
      popViewport()
    }
  }
  
  drawNormalCI <- function(CI.L, ES, CI.U, size) {
    grid.rect(x = unit(ES, "native"), 
              width = unit(size, "snpc"), 
              height = unit(size, "snpc"),
              gp = gpar(fill = "black"))
    
    if (convertX(unit(CI.U, "native"), "npc", valueOnly = TRUE) > 1)
      grid.lines(x = unit(c(CI.L, 1), c("native", "npc")), 
                 y = .5,
                 arrow = arrow(length = unit(0.05, "inches")))
    else { 
      lineCol <- "black"
      grid.lines(x = unit(c(CI.L, CI.U), "native"), 
                 y = 0.5,
                 gp = gpar(col = lineCol))
    }
  }
  
  drawSummaryCI <- function(CI.L, ES, CI.U) {
    grid.polygon(x=unit(c(CI.L, ES, CI.U, ES), "native"),
                 y=unit(0.5 + c(0, 0.25, 0, -0.25), "npc"))
  }
  
  drawDataCol <- function(col, j) { # j = col_place
    pushViewport(viewport(layout.pos.col = j, xscale = col$range))
    grid.lines(x = unit(col$ES[length(col$ES)], "native"),
               y = unit(0:(n_rows - 2), "lines"), 
               gp = gpar(lty = "dashed"))
    grid.xaxis(gp=gpar(cex = 1))
    grid.text("Effect Size",                  
              y = unit(-3, "lines"), 
              x = unit(0.5, "npc"), 
              just = "centre", 
              gp = gpar(fontface = "bold"))
    popViewport()
    x = unit(0.5, "npc")
    
    for (i in 1:length(col$rows)) {
      pushViewport(viewport(layout.pos.row = col$rows[i], 
                            layout.pos.col = j,
                            xscale = col$range))
      if (col$type[i] == "n")
        drawNormalCI(col$CI.L[i], col$ES[i], col$CI.U[i], col$size[i])
      else
        drawSummaryCI(col$CI.L[i], col$ES[i], col$CI.U[i])
      popViewport()
    }
  }
  
  id_col_width   <- max(unit(rep(1, length(id_col$values)), "grobwidth", 
                             id_col$values))
  data_col_width <- unit(3, "inches")
  gap_col        <- unit(10, "mm")
  cols           <- unit.c(id_col_width, gap_col, data_col_width, gap_col)
  
  add_col_widths <- c()
  if (n_user_cols > 1){
    for (i in 1:n_user_cols) {
      add_col_widths[[i]] <-  max(unit(rep(1, length(add_col[[i]]$values)), 
                                       "grobwidth", add_col[[i]]$values))
      cols <- unit.c(cols, add_col_widths[[i]])
      cols <- unit.c(cols, gap_col)
    }
  }
  if (n_user_cols == 1){
    add_col_widths <-  max(unit(rep(1, length(add_col[1]$values)), 
                                "grobwidth", add_col[1]$values))
    cols <- unit.c(cols, add_col_widths[1])
    cols <- unit.c(cols, gap_col)
  }
  
  pushViewport(viewport(layout = grid.layout(n_rows, (4 + (2 * n_user_cols)),
                                             widths = cols,
                                             heights = unit(c(1, rep(1, n_rows)), "lines"))))
  
  pushViewport(viewport(layout.pos.row = 1))
  
  grid.text("Forest Plot", 
            y = unit(+3, "lines"),
            just = "center", 
            gp = gpar(fontface = "bold"))
  
  grid.text(paste(x$model.lab1), 
            y = unit(+2, "lines"),
            just = "center", 
            gp = gpar(fontface = "italic"))
  
  popViewport()
  
  drawLabelCol(id_col, 1)
  
  if (n_user_cols > 1){
    for (i in 1:n_user_cols) {
      drawLabelCol(add_col[[i]], ((i * 2) + 3))
    }
  }
  if (n_user_cols == 1){
    for (i in 1:n_user_cols) {
      drawLabelCol(add_col, 5)
    }
  }
  
  drawDataCol(data_col, 3)
  popViewport()
}