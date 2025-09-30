#' Simulation study driver across parameter grids
#'
#' Runs a series of Monte Carlo simulations across all combinations of user-specified
#' parameters, calling [sim.metaZ()] for each setting, and aggregates AUCs and ROC curves.
#'
#' @param r Integer. Number of simulation replicates per setting (default \code{100}).
#' @param m Integer. Number of markers (rows / observations) per replicate (default \code{1000}).
#' @param d Integer vector. Dimensions to evaluate (default \code{c(5, 10)}).
#' @param perout Numeric vector in \[0,1\]. Proportion of outliers (default \code{c(0.01, 0.05)}).
#' @param farout Numeric vector. Magnitude/severity of outliers (default \code{c(0.25, 0.50)}).
#' @param outlierType Character vector of outlier mechanisms:
#'   \code{"casewise"}, \code{"cellwiseStructured"}, or \code{"both"} (default all three).
#' @param corrType Character vector of correlation structures to use when generating
#'   the covariance matrix via \code{\link[cellWise]{generateCorMat}}
#'   (default \code{c("independent","ALYZ","A09")}).
#' @param u Numeric vector in \[0,1\] for the EDF grid used in ROC construction (default \code{0:10000/10000}).
#' @param seed Optional integer. Random seed passed to \code{set.seed} (default \code{NULL}).
#'
#' @details
#' For each row of the parameter grid (formed by \code{d}, \code{perout}, \code{farout},
#' \code{outlierType}, and \code{corrType}), the function calls [sim.metaZ()], collects
#' the AUCs and ROC matrices, and returns both the per-setting summary table and the
#' list of full per-setting results.
#'
#' @return A list with:
#' \describe{
#'   \item{res.mtx}{\code{data.frame} with one row per parameter setting and columns for
#'     the setting plus AUCs returned by [sim.metaZ()].}
#'   \item{res.list}{List of length equal to the number of settings; each element contains
#'     the covariance matrix \code{S} used and the averaged ROC matrix \code{roc.mtx}.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' out = sim.series(r = 2, m = 200, d = c(5, 10), perout = 0.05, farout = 0.5,
#'                   outlierType = "casewise", corrType = c("independent", "ALYZ"))
#' str(out$res.mtx)
#' }
#'
#' @importFrom cellWise generateCorMat
#' @export
sim.series=function(r=100,
                    m=1000,
                    d=c(5,10),
                    perout=c(0.01,0.05),
                    farout=c(0.25,0.50),
                    outlierType=c("casewise",
                                  "cellwiseStructured",
                                  "both"),
                    corrType=c("independent","ALYZ","A09"),
                    u=0:10000/10000,
                    seed=NULL)

{
  sim.settings=expand.grid(d=d,
                           perout=perout,
                           farout=farout,
                           outlierType=outlierType,
                           corrType=corrType)
  set.seed(seed)
  nsett=nrow(sim.settings)
  res.list=vector("list",nsett)
  auc.mtx=matrix(NA,nsett,3)
  for (i in 1:nsett)
  {
    sim.res=try(sim.metaZ(r=r,m=m,
                          d=sim.settings$d[i],
                          perout=sim.settings$perout[i],
                          farout=sim.settings$farout[i],
                          outlierType=sim.settings$outlierType[i],
                          corrType=sim.settings$corrType[i],
                          u=u))
    if (!inherits(sim.res, "try-error"))
    {
      res.list[[i]]=list(S=sim.res$S,
                         roc.mtx=sim.res$roc.mtx)
      auc.mtx[i,]=sim.res$auc
    }

  }
  colnames(auc.mtx)=names(sim.res$auc)

  res.mtx=cbind.data.frame(sim.settings,auc.mtx)

  full.res=list(res.mtx=res.mtx,
                res.list=res.list)

  return(full.res)
}


#' Single-setting simulation with cwMCD evaluation
#'
#' Simulates data under a chosen correlation structure and outlier mechanism,
#' applies [cwMCD()] to each replicate, builds ROC curves (via EDF grid \code{u}),
#' and averages AUCs/ROCs over \code{r} replicates.
#'
#' @param r Integer. Number of simulation replicates (default \code{100}).
#' @param m Integer. Number of markers (rows / observations) per replicate (default \code{1000}).
#' @param d Integer. Dimension (number of variables) (default \code{5}).
#' @param perout Numeric in \[0,1\]. Proportion of outliers (default \code{0.05}).
#' @param farout Numeric. Magnitude/severity of outliers (default \code{0.5}).
#' @param outlierType Character. One of \code{"casewise"}, \code{"cellwiseStructured"}, or \code{"both"} (default \code{"casewise"}).
#' @param corrType Character. Correlation structure for \code{\link[cellWise]{generateCorMat}}
#'   (default \code{"ALYZ"}). Use \code{"independent"} for identity covariance.
#' @param u Numeric vector in \[0,1\] giving the EDF grid for ROC curve interpolation
#'   (default \code{0:10000/10000}).
#' @param seed Optional integer. Initial random seed (default \code{NULL}).
#'
#' @details
#' The covariance matrix \code{S} is formed by \code{\link[cellWise]{generateCorMat}} unless
#' \code{corrType == "independent"}, in which case \code{S = I_d}. Data are simulated by
#' \code{\link[cellWise]{generateData}}.
#'
#' @return A list with:
#' \describe{
#'   \item{auc}{Named numeric vector of mean AUCs across replicates (one per criterion).}
#'   \item{roc.mtx}{Matrix with columns \code{lvl} (EDF grid) and per-criterion mean ROC ordinates.}
#'   \item{pm}{The last replicate's performance matrix (includes attributes \code{auc}).}
#'   \item{S}{The covariance matrix used.}
#'   \item{m, d, perout, farout, outlierType, seed}{Input parameters echoed back.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' res = sim.metaZ(r = 3, m = 300, d = 4, perout = 0.05, farout = 0.5,
#'                  outlierType = "both", corrType = "A09",
#'                  u = seq(0, 1, length.out = 51), seed = 42)
#' names(res$auc)
#' }
#'
#' @importFrom cellWise generateCorMat generateData
#' @importFrom stats approx
#' @export

sim.metaZ=function(r=100,                  # number of simulation reps
                   m=1000,                 # number of markers
                   d=5,                    # dimension
                   perout=0.05,            # percentage outliers
                   farout=0.5,             # how far out are the outliers
                   outlierType="casewise", # type of outlier
                   corrType="ALYZ",        # correlation type, see generateCorMat of cellWise package
                   u=0:10000/10000,        # u-grid for EDFs
                   seed=NULL)              # initial random seed

{
  seed0=seed
  S=matrix(0,d,d)
  diag(S)=1
  if (corrType!="independent")
    S=cellWise::generateCorMat(d,corrType,seed=seed)

  Zsim=gen.Zmtx(m,d,S,perout,farout,outlierType)
  h=attr(Zsim,"h")
  # mZ.result=metaZ(Zsim$X)
  mZ.result=cwMCD(Zsim$X)
  mZ.result = mZ.result$res.tbl
  pm=evaluate.performance(mZ.result,h)
  this.auc=attr(pm,"auc")

  auc.lbls=names(this.auc)
  lvl.lbls=paste0("lvl_",auc.lbls)
  pwr.lbls=paste0("pwr_",auc.lbls)

  this.roc=matrix(NA,length(u),length(this.auc))
  colnames(this.roc)=auc.lbls
  for (j in 1:length(this.auc))
  {

    this.roc[,auc.lbls[j]]=stats::approx(pm[,lvl.lbls[j]],
                                  pm[,pwr.lbls[j]],
                                  xout=u)$y
  }

  auc=this.auc
  roc.mtx=this.roc

  for (i in 1:(r-1))
  {
    Zsim=gen.Zmtx(m,d,S,perout,farout,outlierType)
    h=attr(Zsim,"h")
    mZ.result=cwMCD(Zsim$X)
    mZ.result=mZ.result$res.tbl
    pm=evaluate.performance(mZ.result,h)
    this.auc=attr(pm,"auc")
    for (j in 1:length(auc))
    {

      this.roc[,j]=approx(pm[,lvl.lbls[j]],
                          pm[,pwr.lbls[j]],
                          xout=u)$y
    }

    auc=auc+this.auc
    roc.mtx=roc.mtx+this.roc


  }

  auc=auc/r
  roc.mtx=roc.mtx/r
  roc.mtx=cbind(u,roc.mtx)
  colnames(roc.mtx)=c("lvl",names(auc))

  res=list(auc=auc,
           roc.mtx=roc.mtx,
           pm=pm,
           S=S,m=m,
           d=d,
           perout=perout,
           farout=farout,
           outlierType=outlierType,
           seed=seed0)

  return(res)
}

#' Generate simulated Z-matrix with injected outliers
#'
#' @description
#' INTERNAL: helper used by [sim.metaZ()]. Generates a matrix of dimension \code{m x d}
#' under multivariate normality with covariance \code{Sigma} and injects outliers according to
#' \code{outlierType}. Returns the simulated list from \code{\link[cellWise]{generateData}}
#' with attribute \code{"h"} marking rows that contain any outlying cell or are outlying rows.
#'
#' @param m Integer. Number of rows (markers).
#' @param d Integer. Number of columns (variables).
#' @param Sigma Numeric \code{d x d} covariance matrix.
#' @param perout Numeric in \[0,1\]. Proportion of outliers.
#' @param farout Numeric. Magnitude/severity of outliers.
#' @param outlierType Character. One of \code{"casewise"}, \code{"cellwiseStructured"}, or \code{"both"}.
#'
#' @return The object returned by \code{\link[cellWise]{generateData}} with an added attribute
#' \code{"h"} (integer vector of length \code{m}, values in \{0,1\}) and possibly modified \code{X}.
#'
#' @keywords internal
#' @importFrom cellWise generateData
#' @importFrom stats rnorm
#' @noRd
#'
gen.Zmtx=function(m,d,Sigma,
                  perout,farout,
                  outlierType)

{
  mu=rep(0,d)
  Zsim=generateData(m,d,
                    mu=mu,
                    Sigma=Sigma,
                    perout,
                    farout,
                    outlierType)

  h=rep(0,m)
  if (!is.null(Zsim$indrows))
  {
    h[Zsim$indrows]=1
    n.out=length(Zsim$indrows)
    Zsim$X[Zsim$indrows,]=Zsim$X[Zsim$indrows,]+rnorm(n.out*d)
  }
  if (!is.null(Zsim$indcells))
  {
    indrows=unique((Zsim$indcells-1)%%m+1)
    h[indrows]=1
    Zsim$X[indrows,]=Zsim$X[indrows,]+rnorm(length(indrows)*d)
  }
  attr(Zsim,"h")=h
  return(Zsim)
}


#' Evaluate performance on a simulated dataset
#'
#' @description
#' INTERNAL: helper used by [sim.metaZ()]. Given a result table (e.g., from [cwMCD()])
#' and a binary hypothesis indicator \code{h} for each row, computes cumulative
#' true-positive and false-positive rates as a function of decreasing criterion, along with AUCs.
#'
#' @param mZ.result \code{data.frame} with per-row statistics and columns ending in \code{".crit"}
#'   (e.g., \code{"fisher.crit"}, \code{"SSz.crit"}, \code{"MCD.crit"}).
#' @param h Integer or logical vector of length \code{nrow(mZ.result)} with 0 for null and 1 for non-null rows.
#'
#' @return A \code{data.frame} that binds \code{mZ.result} and two sets of columns:
#'   \code{"pwr_*"} and \code{"lvl_*"} (one pair per criterion). The object also carries
#'   an \code{auc} attribute: a named numeric vector of AUCs (one per criterion).
#'
#' @details
#' Rows are ranked by each \code{*.crit} in descending order. Empirical TPR (\code{F1})
#' and FPR (\code{F0}) are accumulated, de-duplicated on FPR, and trapezoid-integrated to obtain AUC.
#'
#' @keywords internal
#' @noRd

evaluate.performance=function(mZ.result, # result of metaZ
                              h)            # hypothesis indicator (0=null, 1=non-null)
{
  clms=colnames(mZ.result)
  res=cbind(mZ.result,h=h)
  colnames(res)=c(clms,"h")
  na.clms=colMeans(is.na(res))
  m=nrow(res)
  crit.clms=grep(".crit",clms,fixed=T,value=T)
  pwr=lvl=matrix(NA,m,length(crit.clms))
  auc=rep(NA,length(crit.clms))
  names(auc)=crit.clms
  colnames(pwr)=colnames(lvl)=crit.clms
  nas=na.clms[crit.clms]


  for (cr in crit.clms)
  {
    ord=rev(order(res[,cr]))
    F1=cumsum(res[ord,"h"]%in%1)
    F1=F1/F1[m]
    F0=cumsum(res[ord,"h"]%in%0)
    F0=F0/F0[m]
    u0=!duplicated(F0)
    x=c(0,F0[u0],1)
    y=c(0,F1[u0],1)
    dx=diff(x)
    hy=(y[-1]+y[-length(y)])/2
    auc[cr]=sum(dx*hy)
    pwr[ord,cr]=F1
    lvl[ord,cr]=F0
  }
  colnames(pwr)=paste0("pwr_",crit.clms)
  colnames(lvl)=paste0("lvl_",crit.clms)

  res=cbind(res,pwr,lvl)
  attr(res,"auc")=auc

  return(res)
}
