#' Forest Plots for Robust Variance Estimation Meta-Analysis
#' 
#' \code{forest.robu} In meta-analysis, forest plots provide a graphical
#' depiction of effect size estimates and their corresponding confidence
#' intervals. The \code{forest.robu()} function in \pkg{robumeta} can be used
#' to produce forest plots for RVE meta-analyses. The function requires the
#' \pkg{grid} package and is based on examples provided in (Murrell, 2011). As
#' is the case with traditional forest plots, point estimates of individual
#' effect sizes are plotted as boxes with areas proportional to the weight
#' assigned to that effect size.  Importantly, here the weight is not
#' necessarily proportional to the effect size variance or confidence
#' intervals, since the combined study weight is divided evenly across the
#' study effect sizes. Two-sided 95\% confidence intervals are calculated for
#' each effect size using a standard normal distribution and plotted along with
#' each block. The overall effect is included at the bottom of the plot as a
#' diamond with width equivalent to the confidence interval for the estimated
#' effect. The RVE forest function is designed to provide users with forest
#' plots which display each individual effect size used in the meta-analysis,
#' while taking into account the study- or cluster-level properties inherent to
#' the RVE analysis. As such, the user must specify columns from their original
#' dataset that contain labels for the study or cluster and for the individual
#' effect sizes.
#' 
#' 
#' @param x An intercept-only RVE model previously fit using the \code{robu()}
#' function..
#' @param study.lab A vector of labels to be used to identify study (or
#' cluster) level groupings in the forest plot. For instance, labels for the
#' study column might be author names with corresponding publication years.
#' @param es.lab A vector of labels to be used to individual effect sizes in
#' the forest plot. Labels for individual effect sizes might be ``Math Score''
#' or ``Reading Score'' for a meta-analysis that included such measures or as
#' simple as ``Effect Size 1'' and ``Effect Size 2.''
#' @param ...  Additional arguments to be passed to the forest function. Any
#' number of additional columns can be specified to be plotted along side the
#' confidence interval column and can be specified with the following syntax
#' \code{``arg1'' = ``arg2''} where \code{``arg1''} is the title of the column
#' on the forest plot, and \code{``arg2''} is the name of the column from the
#' original \code{data.frame} that contains the information to be displayed
#' alongside the estimates and confidence intervals.
#' @references
#' 
#' Hedges, L.V., Tipton, E., Johnson, M.C. (2010) Robust variance estimation in
#' meta-regression with dependent effect size estimates. \emph{Research
#' Synthesis Methods}. \bold{1}(1): 39--65. Erratum in \bold{1}(2): 164--165.
#' DOI: 10.1002/jrsm.5
#' 
#' Murrell P (2011). R Graphics. CRC/Taylor & Francis. ISBN 9781439831762.
#' 
#' Tipton, E. (in press) Small sample adjustments for robust variance
#' estimation with meta-regression. \emph{Psychological Methods}.
#' @keywords forest.robu
#' @examples
#' 
#' 
#' # Load data
#' data(oswald2013.ex1)
#' 
#' # Run intercept only model.
#' oswald_intercept <- robu(formula = effect.size ~ 1, data = oswald2013.ex1, 
#'                          studynum = Study, var.eff.size = var.eff.size, 
#'                          rho = 0.8, small = TRUE)
#' 
#' # Create forest plot. 
#' forest.robu(oswald_intercept, es.lab = "Crit.Cat", study.lab = "Study", 
#'             "Effect Size" = effect.size, # optional column
#'             "Weight" = r.weights)        # optional column
#' 
#' @export 
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
  data$study_rows <- as.numeric(stats::ave(data$es_rows, data$study_num, FUN = min)- 1)
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
      t <- lapply(values, function(x) grid::textGrob(paste(x), x = 0, just = just))
    } else {
      t <- lapply(values, function(x) grid::textGrob(paste(x), x = 0, just = just,
                                               gp = grid::gpar(fontface = "bold")))
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
      grid::pushViewport(grid::viewport(layout.pos.row = col$rows[i], layout.pos.col = j))
      grid::grid.draw(col$values[[i]])
      grid::popViewport()
    }
  }
  
  drawNormalCI <- function(CI.L, ES, CI.U, size) {
    grid::grid.rect(x = grid::unit(ES, "native"), 
              width = grid::unit(size, "snpc"), 
              height = grid::unit(size, "snpc"),
              gp = grid::gpar(fill = "black"))
    
    if (grid::convertX(grid::unit(CI.U, "native"), "npc", valueOnly = TRUE) > 1)
      grid::grid.lines(x = grid::unit(c(CI.L, 1), c("native", "npc")), 
                 y = .5,
                 arrow = grid::arrow(length = grid::unit(0.05, "inches")))
    else { 
      lineCol <- "black"
      grid::grid.lines(x = grid::unit(c(CI.L, CI.U), "native"), 
                 y = 0.5,
                 gp = grid::gpar(col = lineCol))
    }
  }
  
  drawSummaryCI <- function(CI.L, ES, CI.U) {
    grid::grid.polygon(x=grid::unit(c(CI.L, ES, CI.U, ES), "native"),
                 y=grid::unit(0.5 + c(0, 0.25, 0, -0.25), "npc"))
  }
  
  drawDataCol <- function(col, j) { # j = col_place
    grid::pushViewport(grid::viewport(layout.pos.col = j, xscale = col$range))
    grid::grid.lines(x = grid::unit(col$ES[length(col$ES)], "native"),
               y = grid::unit(0:(n_rows - 2), "lines"), 
               gp = grid::gpar(lty = "dashed"))
    grid::grid.xaxis(gp=grid::gpar(cex = 1))
    grid::grid.text("Effect Size",                  
              y = grid::unit(-3, "lines"), 
              x = grid::unit(0.5, "npc"), 
              just = "centre", 
              gp = grid::gpar(fontface = "bold"))
    grid::popViewport()
    x = grid::unit(0.5, "npc")
    
    for (i in 1:length(col$rows)) {
      grid::pushViewport(grid::viewport(layout.pos.row = col$rows[i], 
                            layout.pos.col = j,
                            xscale = col$range))
      if (col$type[i] == "n")
        drawNormalCI(col$CI.L[i], col$ES[i], col$CI.U[i], col$size[i])
      else
        drawSummaryCI(col$CI.L[i], col$ES[i], col$CI.U[i])
      grid::popViewport()
    }
  }
  
  id_col_width   <- max(grid::unit(rep(1, length(id_col$values)), "grobwidth", 
                             id_col$values))
  data_col_width <- grid::unit(3, "inches")
  gap_col        <- grid::unit(10, "mm")
  cols           <- grid::unit.c(id_col_width, gap_col, data_col_width, gap_col)
  
  add_col_widths <- c()
  if (n_user_cols > 1){
    for (i in 1:n_user_cols) {
      add_col_widths[[i]] <-  max(grid::unit(rep(1, length(add_col[[i]]$values)), 
                                       "grobwidth", add_col[[i]]$values))
      cols <- grid::unit.c(cols, add_col_widths[[i]])
      cols <- grid::unit.c(cols, gap_col)
    }
  }
  if (n_user_cols == 1){
    add_col_widths <-  max(grid::unit(rep(1, length(add_col[1]$values)), 
                                "grobwidth", add_col[1]$values))
    cols <- grid::unit.c(cols, add_col_widths[1])
    cols <- grid::unit.c(cols, gap_col)
  }
  
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n_rows, (4 + (2 * n_user_cols)),
                                             widths = cols,
                                             heights = grid::unit(c(1, rep(1, n_rows)), "lines"))))
  
  grid::pushViewport(grid::viewport(layout.pos.row = 1))
  
  grid::grid.text("Forest Plot", 
            y = grid::unit(+3, "lines"),
            just = "center", 
            gp = grid::gpar(fontface = "bold"))
  
  grid::grid.text(paste(x$model.lab1), 
            y = grid::unit(+2, "lines"),
            just = "center", 
            gp = grid::gpar(fontface = "italic"))
  
  grid::popViewport()
  
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
  grid::popViewport()
}
