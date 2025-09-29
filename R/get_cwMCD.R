#' Cellwise Minimum Covariance Determinant (cwMCD) analysis
#'
#' Runs the cellwise MCD procedure on a numeric data matrix \code{X}, returning
#' per-row outlier diagnostics including Fisher’s method, sum-of-squares
#' statistics, and MCD-based Mahalanobis distances.
#'
#' @param X A numeric matrix or data frame. Rows are observations and columns
#'   are variables. If row or column names are missing, they are automatically
#'   assigned.
#' @param alpha Numeric scalar (default \code{0.75}). Robustness parameter passed
#'   to \code{\link{cellMCD}}.
#' @param quant Numeric scalar (default \code{0.99}). Quantile cutoff passed to
#'   \code{\link{cellMCD}}.
#' @param crit Numeric scalar (default \code{1e-4}). Convergence criterion.
#' @param noCits Integer (default \code{100}). Maximum number of iterations.
#' @param lmin Numeric scalar (default \code{1e-4}). Minimum eigenvalue bound.
#' @param checkPars List of optional arguments passed to \code{\link{cellMCD}}
#'   (default \code{list()}).
#'
#' The output table includes original input columns (prefixed with
#' \code{"input."}), row-level summary statistics, and per-variable MCD Z-scores
#' (prefixed with \code{"MCD.Z."}).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{cellmcd.obj}{The raw object returned by \code{\link{cellMCD}},
#'     containing robust location, scatter, and residuals.}
#'   \item{res.tbl}{A \code{data.frame} with one row per observation, containing: 
#'     \describe{
#'       \item{input.*}{Original variables.}
#'       \item{max.abs.X.column}{Column with largest absolute input per row.}
#'       \item{max.abs.X.value}{Value of largest absolute input per row.}
#'       \item{fisher.stat, fisher.crit, fisher.pval}{Fisher’s combined test results.}
#'       \item{SSz.stat, SSz.crit, SSz.pval}{Sum-of-squares Z-statistics.}
#'       \item{MCD.stat, MCD.crit, MCD.pval}{Mahalanobis distances and p-values.}
#'       \item{max.MCD.Z.column, max.MCD.Z.value}{Maximum cwMCD Z-score per row.}
#'       \item{MCD.Z.*}{Per-variable cwMCD Z-scores.}
#'     }
#'   }
#' }
#'
#' @seealso \code{\link[cellwise]{cellMCD}}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' X <- matrix(rnorm(100), nrow = 20, ncol = 5)
#' res <- cwMCD(X)$res.tbl
#' head(res)
#' }
#'
#' @importFrom cellwise cellMCD
#' @export
#' 
cwMCD=function(X,alpha=0.75,quant=0.99,
               crit=1e-4,noCits=100,lmin=1e-4,
               checkPars=list())
  
{
  if (is.null(rownames(X))) rownames(X)=paste0("row_",1:nrow(X))
  if (is.null(colnames(X))) colnames(X)=paste0("clm_",1:ncol(X))
  
  # cellMCD
  cw.res=cellMCD(X,alpha,quant,crit,noCits,lmin,checkPars)
  mhd=mahalanobis(X,cw.res$mu,cw.res$S)
  p.mhd=pchisq(mhd,ncol(X),lower.tail=F)
  mhc=sqrt(mhd/ncol(X))
  
  # maximum MCD Z stat
  which.max.Z=apply(abs(cw.res$Zres),1,which.max)
  max.mcd.Z.clm=colnames(X)[which.max.Z]
  max.mcd.Z.value=rep(NA,nrow(X))
  for (i in 1:nrow(X))
    max.mcd.Z.value[i]=cw.res$Zres[i,which.max.Z[i]]
  
  
  # maximum input X
  which.max.X=apply(abs(X),1,which.max)
  max.input.X.clm=colnames(X)[which.max.X]
  max.input.X.value=rep(NA,nrow(X))
  for (i in 1:nrow(X))
    max.input.X.value[i]=X[i,max.input.X.clm[i]]
  
  
  # Fisher's test
  P=2*pnorm(-abs(X))
  fisher.stat=rowSums(-2*log(P))
  fisher.crit=sqrt(fisher.stat/(2*ncol(P)))
  p.fisher=pchisq(fisher.stat,2*ncol(P),lower.tail=F)
  
  # sum of squared z-statistics
  ssZ.stat=rowSums(X^2)
  ssZ.crit=sqrt(ssZ.stat/ncol(P))
  ssZ.pval=pchisq(ssZ.stat,ncol(P),lower.tail=F)
  
  
  res.tbl=cbind.data.frame(X,
                           max.abs.X.column=max.input.X.clm,
                           max.abs.X.value=max.input.X.value,
                           fisher.stat=fisher.stat,
                           fisher.crit=fisher.crit,
                           fisher.pval=p.fisher,
                           SSz.stat=ssZ.stat,
                           SSz.crit=ssZ.crit,
                           ssZ.pval=ssZ.pval,
                           mcd.dist=mhd,
                           mcd.crit=mhc,
                           mcd.pval=p.mhd,
                           max.MCD.Z.column=max.mcd.Z.clm,
                           max.MCD.Z.value=max.mcd.Z.value,
                           Z.mcd=cw.res$Zres)
  
  rownames(res.tbl)=rownames(X)
  colnames(res.tbl)=c(paste0("input.",colnames(X)),
                      "max.abs.X.column","max.abs.X.value",
                      "fisher.stat","fisher.crit","fisher.pval",
                      "SSz.stat","SSz.crit","SSz.pval",
                      "MCD.stat","MCD.crit","MCD.pval",
                      "max.MCD.Z.column","max.MCD.Z.value",
                      paste0("MCD.Z.",colnames(X)))

  res.list = list(
    cellmcd.obj = cw.res,
    res.tbl = res.tbl
  )
  
  return(res.list)
  
}


#' Draw comparison of elliptical contours for outlier methods
#'
#' Plots the theoretical rejection contours from three outlier-detection methods
#' (Fisher’s method, sum-of-squared Z statistics, and cellMCD) in the
#' \eqn{z_1}-\eqn{z_2} plane, optionally overlaying simulated bivariate normal
#' data points with correlation \code{rho}.
#'
#' @param alpha Significance level for contours (default: 0.05).
#' @param rho Correlation parameter for the bivariate normal used in both the
#'   simulation and the theoretical cellMCD ellipse (default 0.5).
#' @param n Number of simulated points to generate from a bivariate normal with
#'   correlation \code{rho} (default \code{100}). Set to \code{0} to suppress
#'   simulation.
#'
#' @return A numeric matrix with two columns (\code{x1}, \code{x2}) containing
#'   the simulated data points if \code{n > 0}.
#'
#' @importFrom ellipse ellipse
#' @importFrom mvtnorm rmvnorm
#' @seealso \code{\link[cellwise]{cellMCD}}, \code{\link{qq.uniform}}
#'
#' @examples
#' \dontrun{
#' # Default with n=100 simulated points
#' sim.ellipse.figure()
#'
#' # Only draw contours, no points
#' sim.ellipse.figure(n = 0)
#'
#' # Increase correlation
#' sim.ellipse.figure(rho = 0.8, n = 200)
#' }
#'
#' @export
#' 
sim.ellipse.figure=function(alpha=0.05,rho=0.50,n=100, ...)
{
  
  # elliptical contour
  
  # 1-alpha coordinates for sum of squared z-stats method
  z1=seq(from=-sqrt(qchisq(alpha,2,lower.tail=F)),
         to=sqrt(qchisq(alpha,2,lower.tail=F)),
         length=10000)
  z2=sqrt(qchisq(alpha,2,lower.tail=F)-z1^2)
  z2[1]=0
  z2[10000]=0
  z1.st=z1
  z2.st=z2
  
  
  # 1-alpha coordinates for Fisher's method
  lg1=seq(from=0,to=-0.5*qchisq(alpha,4,lower.tail=F),length=10000)
  lg2=-0.5*qchisq(alpha,4,lower.tail=F)-lg1
  p1=exp(lg1)
  p2=exp(lg2)
  z1=qnorm(p1/2)
  z2=qnorm(p2/2)
  z1.fs=c(z1)
  z2.fs=c(z2)
  
  x1=x2=X=NULL
  if (n>0)
  {
    S=matrix(0,2,2)
    diag(S)=1
    S[1,2]=S[2,1]=rho
    X=rmvnorm(n,rep(0,2),S)
    x1=X[,1]; x2=X[,2]
    
  }
  
  
  cw.res=cwMCD(cbind(x1,x2))$mcd.obj
  cw.rho=cw.res$S[1,2]
  
  ell=ellipse::ellipse(x=rho,
                       level=1-alpha)
  
  plot(1.1*range(z1.st,z1.fs,ell[,1],x1),
       1.1*range(z2.st,z2.fs,ell[,2],x2),
       xlim=c(-4,4),ylim=c(-4,4),
       type="n",xlab="z-stat 1",ylab="z-stat 2",las=1)
  if (n>0) points(x1,x2,pch=19,col="darkgray",cex=0.50)
  lines(ell,col="darkred",lwd=3)
  lines(z1.fs,z2.fs,col="darkblue",lwd=3)  
  lines(z1.fs,-z2.fs,col="darkblue",lwd=3)
  lines(-rev(z1.fs),rev(z2.fs),col="darkblue",lwd=3)
  lines(-rev(z1.fs),-rev(z2.fs),col="darkblue",lwd=3)
  lines(z1.st,z2.st,col="goldenrod",lwd=3)
  lines(z1.st,-z2.st,col="goldenrod",lwd=3)
  
  return(cbind(x1,x2))
  
}

#' 
#' 
#' #' Uniform QQ plot for cellMCD results
#' #'
#' #' Produces a uniform QQ plot (EDF vs p-value) for a `cwMCD()` result, 
#' #' and returns plot data and RMSE.
#' #' 
#' #' @param mcd   Result from cwMCD(), expected to contain p-value columns:
#' #'              "fisher.pval", "SSz.pval", and "MCD.pval".
#' #' @param p     Threshold for inset region (default 0.01).
#' #' @param ilim  Range for inset region (default c(0, 0.20)).
#' #' @param main  Optional plot title.
#' #' @return A list with:
#' #'   - data: a data.frame of QQ points (method, pval, edf)
#' #'   - rmse: RMSE vs uniform for each method (full range, lwr = 0)
#' #'   - params: list of p and ilim used
#' #'   
#' #'   #' @examples
#' #' \dontrun{
#' #' # mcd1 <- cwMCD(ellA)
#' #' out <- qq.uniform(mcd1, p = 0.01, ilim = c(0, 0.20))
#' #' head(out$data)
#' #' out$rmse
#' #' }
#' #'
#' #' @export
#' #' 
#' qq.uniform <- function(mcd, p = 0.01, ilim = c(0, 0.20), main = NULL, ...) {
#'   
#'   # extract available series
#'   has <- colnames(mcd)
#'   stopifnot("fisher.pval" %in% has)
#'   stopifnot("SSz.pval"     %in% has)
#'   stopifnot("MCD.pval"     %in% has)
#'   stopifnot(p > ilim[1] & p <= ilim[2])
#'   # colors
#'   col_fish <- "darkblue"; col_ssz <- "goldenrod"; col_mcd <- "darkred"
#'   
#'   # Uniform quantile-quantile plot
#'   plot(c(0,1),c(0,1),
#'        xlab="row-outlier p-value",
#'        ylab="EDF",las=1,type="n", main = main)
#'   
#'   # shaded lower-left box up to min(p, ilim max)
#'   
#'   pdata = data.frame(fisher.pval = mcd$fisher.pval,
#'                      fisher.edf  = rank(mcd$fisher.pval)/nrow(mcd),
#'                      MCD.pval    = mcd$MCD.pval,
#'                      MCD.edf  = rank(mcd$MCD.pval)/nrow(mcd),
#'                      SSz.pval = mcd$SSz.pval,
#'                      SSz.edf  = rank(mcd$SSz.pval)/nrow(mcd)
#'   )
#'   rownames(pdata) = rownames(mcd)
#'   
#'   lower_box <- min(p, max(ilim))
#'   rect(xleft = ilim[1], ybottom = ilim[1], xright = lower_box, ytop = lower_box, col = "gainsboro")
#'   points(pdata$fisher.pval,pdata$fisher.edf,col=col_fish,cex=0.5,pch=19)
#'   points(pdata$SSz.pval,pdata$SSz.edf,col=col_ssz,cex=0.5,pch=19)
#'   points(pdata$MCD.pval,pdata$MCD.edf,col=col_mcd,pch=19,cex=0.5)
#'   lines(c(0,1),c(0,1),col="black",lwd=2,lty=2)
#'   
#'   rect(-0.02,0.58,0.42,1.02,col="gainsboro")
#'   a=lower_box
#'   fish.ok=(mcd$fisher.pval<a)&(rank(mcd$fisher.pval)<(a*nrow(mcd)))
#'   mcd.ok=(mcd$MCD.pval<a)&(rank(mcd$MCD.pval)<(a*nrow(mcd)))
#'   ssz.ok=(mcd$SSz.pval<a)&(rank(mcd$SSz.pval)<(a*nrow(mcd)))
#'   points(mcd$fisher.pval[fish.ok]*0.4/a,
#'          0.6+rank(mcd$fisher.pval)[fish.ok]/nrow(mcd)*0.4/a,
#'          col=col_fish,pch=19,cex=0.5)
#'   points(mcd$SSz.pval[ssz.ok]*0.4/a,
#'          0.6+0.4/a*rank(mcd$SSz.pval)[ssz.ok]/nrow(mcd),
#'          col=col_ssz,cex=0.5)
#'   points(mcd$MCD.pval[mcd.ok]*0.4/a,
#'          0.6+0.4/a*rank(mcd$MCD.pval)[mcd.ok]/nrow(mcd),
#'          col=col_mcd,pch=19,cex=0.5)
#'   lines(c(0,0.4),c(0.6,1),col="black",lwd=2,lty=2)
#'   
#'   legend(0.75,0.20,
#'          lwd=2,lty=c(1,1,1,2),
#'          col=c("darkblue","goldenrod","darkred","black"),
#'          legend=c("Fisher","SSz","cellMCD","Uniform"),
#'          xjust=0,yjust=1,cex=0.75)
#'   
#'   rmse.vec = cbind(
#'     'Fisher' = uqq.rmse(mcd$fisher.pval),
#'     'SSz'    = uqq.rmse(mcd$SSz.pval),
#'     'cellMCD'    = uqq.rmse(mcd$MCD.pval)
#'   )
#'   
#'   rlist = list(
#'     data = pdata,
#'     rmse = rmse.vec,
#'     params = list(p = p, inset.limits = ilim)
#'   )
#'   
#'   return(rlist)
#'   
#'   
#' }
#' 
#' #' Compute RMSE of uniform QQ deviation
#' #'
#' #' Calculates the root mean squared error (RMSE) between observed p-values
#' #' and the expected uniform(0,1) distribution in a QQ plot.
#' #'
#' #' @param p Numeric vector of p-values.
#' #' @param lwr Lower bound to truncate the p-value range (default = 0).
#' #'
#' #' @return A numeric scalar giving the RMSE.
#' #' @examples
#' #' set.seed(1)
#' #' pvals <- runif(100)
#' #' uqq_rmse(pvals)
#' #' uqq_rmse(pvals, lwr = 0.1)
#' #'
#' #' @export
#' uqq.rmse=function(p,lwr=0)
#' {
#'   na=is.na(p)
#'   p=p[!na]
#'   ok=(p>=lwr)
#'   p=p[ok]
#'   p=(p-lwr)/(1-lwr)
#'   m=length(p)
#'   res=sqrt(mean((p-rank(p)/(m+1))^2))
#'   return(res)
#' }
#' 
#' 

