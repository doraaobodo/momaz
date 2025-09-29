# TODO: Fix examples using simulation results 

#' Bivariate scatter with SSz, Fisher, and cellMCD contours
#'
#' Draws a bivariate scatterplot for a chosen pair of variables from a cwMCD result
#' and overlays: the SSz circular boundary, Fisher’s method boundary, and the
#' cellMCD ellipse at the given inclusion level (default 0.99).
#'
#' @param obj A cwMCD result. The full list returned by \code{cwMCD()} with
#'   elements \code{$cellmcd.obj} and \code{$res.tbl}.
#' @param vars Character vector of length 2 giving the two \strong{input} column names. 
#'   You may also pass integer column
#'   indices (referring to the \code{input.*} columns inside \code{res.tbl}).
#' @param qnt Numeric in (0,1). Inclusion level for contours (default \code{0.99}).
#' @param xlim,ylim Numeric length-2 vectors for plot limits (default \code{c(-5,5)}).
#' @param pch,col,cex Base plot point style (passed to \code{points()}).
#' @param ... Additional arguments passed to \code{plot()}.
#'
#' @details
#' SSz boundary is the circle of radius \eqn{\sqrt{qchisq(qnt, 2)}}.
#' Fisher boundary is the level set of \eqn{-2\sum\log(p_i)} at the same quantile
#' (converted back to z via \eqn{z = \Phi^{-1}(p/2)}), drawn symmetrically in all quadrants.
#' The cellMCD boundary is an ellipse based on the robust covariance \code{S}
#' returned by \code{cellwise::cellMCD}, drawn with \code{ellipse::ellipse}.
#'
#' @return (Invisibly) a list with the x/y vectors used to draw the three contours.
#'
#' @examples
#' \dontrun{
#' # mcd <- cwMCD(X)           # your matrix/data as in the manuscript
#' plot.bivariate(mcd, vars = c(1,2), qnt = 0.99)
#' }
#'
#' @importFrom ellipse ellipse
#' @importFrom stats qchisq qnorm
#' @importFrom graphics plot points lines polygon
#' @export
plot.bivariate <- function(
    obj, vars, qnt = 0.99,
    xlim = c(-5, 5), ylim = c(-5, 5),
    pch = 19, col = "darkgray", cex = 1.05,
    ...
) {
  # Accept either cwMCD list or just res.tbl
  res.tbl <- if (is.list(obj) && !is.null(obj$res.tbl)) obj$res.tbl else obj
  cw.obj  <- if (is.list(obj) && !is.null(obj$cellmcd.obj)) obj$cellmcd.obj else NULL
  
  # Locate input.* columns
  input_idx <- grep("^input\\.", colnames(res.tbl))
  stopifnot(length(input_idx) >= 2)
  
  # Resolve vars (names or integer positions within input columns)
  if (is.numeric(vars)) {
    vars <- colnames(res.tbl)[input_idx][vars]
  }
  stopifnot(all(vars %in% colnames(res.tbl)))
  
  x <- res.tbl[[vars[1]]]
  y <- res.tbl[[vars[2]]]
  
  # ---- Base scatter ----
  plot(x, y, type = "n", xlim = xlim, ylim = ylim,
       xlab = sub("^input\\.", "z ", vars[1]),
       ylab = sub("^input\\.", "z ", vars[2]),
       las = 1, ...)
  points(x, y, pch = pch, col = col, cex = cex)
  
  # ---- SSz circle ----
  th <- seq(0, 2 * pi, by = 0.01)
  R  <- sqrt(qchisq(qnt, df = 2))
  x.circle <- R * sin(th)
  y.circle <- R * cos(th)
  polygon(x.circle, y.circle, border = "goldenrod", col = NA, lwd = 3)
  
  # ---- Fisher bounds ----
  lg1 <- seq(from = 0, to = -0.5 * qchisq(qnt, df = 4), length.out = 10000)
  lg2 <- -0.5 * qchisq(qnt, df = 4) - lg1
  p1  <- exp(lg1); p2 <- exp(lg2)
  z1  <- qnorm(p1 / 2); z2 <- qnorm(p2 / 2)
  lines(z1,  z2,  col = "darkblue", lwd = 3)
  lines(z1, -z2,  col = "darkblue", lwd = 3)
  lines(-rev(z1),  rev(z2),  col = "darkblue", lwd = 3)
  lines(-rev(z1), -rev(z2),  col = "darkblue", lwd = 3)
  
  # ---- cellMCD ellipse ----
  # If we have the cw object, use its S; otherwise estimate from the plotted pair
  if (!is.null(cw.obj)) {
    # figure out indices in the original S that correspond to the two inputs
    all_in <- colnames(res.tbl)[input_idx]
    idx1 <- match(vars[1], all_in)
    idx2 <- match(vars[2], all_in)
    # robust variances & correlation from S
    S <- cw.obj$S
    v <- diag(S)[c(idx1, idx2)]
    rho <- S[idx1, idx2] / sqrt(prod(v))
    ell <- ellipse::ellipse(x = rho, scale = sqrt(v), level = qnt)
  } else {
    # fallback: assume unit marginals and zero rho (safe default)
    ell <- ellipse::ellipse(x = 0, level = qnt)
  }
  polygon(ell, border = "darkred", lwd = 3)
  
  invisible(list(
    ssz = list(x = x.circle, y = y.circle),
    fisher = list(x = c(z1, z1, -rev(z1), -rev(z1)),
                  y = c(z2, -z2, rev(z2), -rev(z2))),
    mcd = list(x = ell[, 1], y = ell[, 2])
  ))
}

#' Uniform QQ plot for cwMCD results
#'
#' Produces a uniform QQ plot (EDF vs p-value) highlighting the range near zero.
#' Accepts either the full \code{cwMCD()} list (\code{$res.tbl} inside) or the table itself.
#'
#' @param obj A cwMCD result list or its \code{res.tbl}.
#' @param thresh Numeric; right edge of x-range for plotting (default \code{0.20}).
#' @param p Numeric; p-value threshold for “non-sig” filtering in the scatter
#'   (default \code{0.01}).
#' @param ... Passed to \code{plot()}.
#'
#'
#' @examples
#' \dontrun{
#' mcd = cwMCD(X)
#' plot.qq.uniform(mcd, p = 0.01, thresh = 0.20)
#' }
#'
#' @importFrom graphics plot points lines mtext
#' @export
plot.qq.uniform <- function(obj,  p = 0.01, thresh = 0.20, ...) {
  
  res.tbl <- if (is.list(obj) && !is.null(obj$res.tbl)) obj$res.tbl else obj
  
  fsig=(res.tbl$fisher.pval<p)
  cwsg=c(res.tbl$MCD.pval<p)
  sssg=(res.tbl$SSz.pval<p)
  
  a=thresh
  plot(c(0,a),c(0,a),type="n",xlab="p",ylab="rank/N",las=1, ...)
  points(res.tbl$MCD.pval[(res.tbl$MCD.pval<a)&!cwsg],rank(res.tbl$MCD.pval[(res.tbl$MCD.pval<a)&!cwsg])/sum(!cwsg),col="darkred",pch=19)
  points(res.tbl$fisher.pval[(res.tbl$fisher.pval<a)&!fsig],rank(res.tbl$fisher.pval[(res.tbl$fisher.pval<a)&!fsig])/sum(!fsig),col="darkblue",pch=19)
  points(res.tbl$SSz.pval[(res.tbl$SSz.pval<a)&!sssg],rank(res.tbl$SSz.pval[(res.tbl$SSz.pval<a)&!sssg])/sum(!sssg),col="goldenrod",pch=19,cex=0.75)
  lines(c(p,a),c(p,a),col="darkgray",lwd=2)
  
  rmse.vec = cbind(
    'Fisher' = uqq.rmse(res.tbl$fisher.pval),
    'SSz'    = uqq.rmse(res.tbl$SSz.pval),
    'cellMCD'    = uqq.rmse(res.tbl$MCD.pva)
  )
  
  
  return(rmse.vec)
  
}

#' Compute RMSE of uniform QQ deviation
#'
#' Calculates the root mean squared error (RMSE) between observed p-values
#' and the expected uniform(0,1) distribution in a QQ plot.
#'
#' @param p Numeric vector of p-values.
#' @param lwr Lower bound to truncate the p-value range (default = 0).
#'
#' @return A numeric scalar giving the RMSE.
#' @examples
#' set.seed(1)
#' pvals <- runif(100)
#' uqq.rmse(pvals)
#' uqq.rmse(pvals, lwr = 0.1)
#'
#' @export
uqq.rmse=function(p,lwr=0)
{
  na=is.na(p)
  p=p[!na]
  ok=(p>=lwr)
  p=p[ok]
  p=(p-lwr)/(1-lwr)
  m=length(p)
  res=sqrt(mean((p-rank(p)/(m+1))^2))
  return(res)
}

#' Heatmap of top rows under selected methods (panel F style)
#'
#' Draws adjacent heatmaps of the top \code{top_n} rows according to the chosen
#' method(s), using the \code{input.*} columns from a cwMCD result. Rows within each
#' block are ordered by the projection score used in the manuscript.
#'
#' @param obj A cwMCD result list or its \code{res.tbl}.
#' @param methods Character vector among \code{c("MCD","Fisher","SSz")}.
#'   Order controls the left-to-right arrangement (default: all three).
#' @param top_n Integer; number of top rows per method (default \code{25}).
#' @param palette Character; name of \code{hcl.colors()} palette
#'   (default \code{"Green-Orange"}).
#' @param show_rownames Logical; whether to draw row names next to each block (default \code{FALSE}).
#' @param col_labels Optional character vector to draw under each block (default \code{NULL}).
#'
#' @return Invisibly returns a list with the matrices used for each method block.
#'
#' @examples
#' \dontrun{
#' plot.heatmap(mcd, methods = c("MCD","Fisher","SSz"), top_n = 25)
#' }
#'
#' @importFrom graphics image text mtext par
#' @export
plot.heatmap <- function(
    obj,
    methods = c("MCD", "Fisher", "SSz"),
    top_n = 25,
    palette = "Green-Orange",
    show_rownames = FALSE,
    col_labels = NULL
) {
  
  res.tbl <- if (is.list(obj) && !is.null(obj$res.tbl)) obj$res.tbl else obj
  
  # columns
  input_clms <- grep("^input\\.", colnames(res.tbl))
  stopifnot(length(input_clms) >= 2)
  
  # helper to get top rows & order
  prep_block <- function(score_col, decreasing = TRUE) {
    idx <- if (decreasing) order(res.tbl[[score_col]])[1:top_n]
    else rev(order(res.tbl[[score_col]])[1:top_n])
    M <- as.matrix(res.tbl[idx, input_clms, drop = FALSE])
    dp <- rowSums(M) / sqrt(ncol(M) * rowSums(M^2)) # order val
    M[order(dp), , drop = FALSE]
  }
  
  mats <- list()
  if ("MCD" %in% methods)    mats$MCD    <- prep_block("MCD.pval", decreasing = TRUE)
  if ("Fisher" %in% methods) mats$Fisher <- prep_block("fisher.pval", decreasing = FALSE) 
  if ("SSz" %in% methods)    mats$SSz    <- prep_block("SSz.pval", decreasing = FALSE)
  mats <- mats[names(mats) %in% methods]  # preserve requested order
  if (!length(mats)) return(invisible(mats))
  
  breaks = c(-10, -2.57, 0, 2.57, 10)
  pal <- hcl.colors(length(breaks) - 1, palette = palette)
  
  
  # layout: one column per method + 1 narrow colorbar
  k <- length(mats)
  lay <- matrix(seq_len(k + 1), nrow = 1)
  widths <- c(rep(4, k), 1)  # last column is legend
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  layout(lay, widths = widths)

  # draw each method panel
  for (m in names(mats)) {
    M <- mats[[m]]
    nr <- nrow(M); nc <- ncol(M)
    
    par(mar = c(3, 3, 3, 1))
    plot.new()
    # create panel coordinates [0,1] x [0,nr]
    plot.window(xlim = c(0, 1), ylim = c(0, nr))
    # x grid like manuscript (equal-width tiles)
    x_grid <- seq(0, 1, length.out = nc + 1)
    image(z = t(M), x = x_grid, y = 0:nr, col = pal, breaks = breaks, add = TRUE)
    
    # title and bottom labels
    text(0.5, nr, labels = m, pos = 3, cex = 1.2)
    if (!is.null(col_labels) && length(col_labels) == 3) {
      # place A/D/E near the bottom thirds
      thirds <- (x_grid[-1] + x_grid[-(nc + 1)]) / 2
      pos <- thirds[c(1, floor(nc/2), nc)]
      text(pos, 0, labels = col_labels, pos = 1, cex = 1.2)
    }
    
    # optional rownames
    if (show_rownames) {
      rn <- rownames(M); rn[is.na(rn)] <- ""
      text(x = -0.01, y = seq_len(nr) - 0.5, labels = rn, adj = 1, xpd = NA, cex = 1)
    }
    
  }
  
  breaks.cut = seq(-2, 2, by = 1)*2.57
  pal.cut <- hcl.colors(length(breaks) - 1, palette = palette)
  
  par(mar = c(3, 1, 3, 3))
  plot.new(); plot.window(xlim = c(0,1), ylim = range(breaks.cut))
  rect(0, head(breaks.cut, -1), 1, tail(breaks.cut, -1), 
       col = pal.cut, border = NA)
  axis(4, at = breaks.cut[2:4], tick = F)

  invisible(mats)
}

