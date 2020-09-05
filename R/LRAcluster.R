#' @name LRAcluster
#' @aliases LRAcluster
#' @title integrated analysis of cancer omics data by low-rank approximation
#' @description The LRAcluster function is the main interface of this package, it gets a list of matrix as input and outputs the coordinates of the samples in the reduced space and the explained potential.
#' @param data a list of data matrix,please keep the columns are the same order of samples
#' @param types a list of data types, can be binary, gaussian, or poisson
#' @param dimension the reduced dimension
#' @param names data names
#' @return A list contains the following component:
#'
#'         \code{coordinate} A matrix of the coordinates of all the samples in the reduced space
#'
#'         \code{potential}  ratio of explained variance
#' @keywords internal
#' @author Dingming Wu, Dongfang Wang
#' @references Wu D, Wang D, Zhang MQ, Gu J (2015). Fast dimension reduction and integrative clustering of multi-omics data using low-rank approximation: application to cancer molecular classification. BMC Genomics, 16(1):1022.
LRAcluster <- function(data, types, dimension = 2, names = as.character(1:length(data)))
{
  #--------#
  # binary #
  #--------#
  epsilon.binary<-2.0
  check.binary.row<-function(arr)
  {
    if (sum(!is.na(arr))==0)
    {
      return (F)
    }
    else
    {
      idx<-!is.na(arr)
      if (sum(arr[idx])==0 || sum(arr[idx])==sum(idx))
      {
        return (F)
      }
      else
      {
        return (T)
      }
    }
  }
  check.binary<-function(mat,name)
  {
    index<-apply(mat,1,check.binary.row)
    n<-sum(!index)
    if (n>0)
    {
      w<-paste("Warning: ",name," have ",as.character(n)," invalid lines!",sep="")
      warning(w)
    }
    mat_c<-mat[index,]
    rownames(mat_c)<-rownames(mat)[index]
    colnames(mat_c)<-colnames(mat)
    return (mat_c)
  }

  base.binary.row<-function(arr)
  {
    idx<-!is.na(arr)
    n<-sum(idx)
    m<-sum(arr[idx])
    return (log(m/(n-m)))
  }

  base.binary<-function(mat)
  {
    mat_b<-matrix(0,nrow(mat),ncol(mat))
    ar_b<-apply(mat,1,base.binary.row)
    mat_b[1:nrow(mat_b),]<-ar_b
    return (mat_b)
  }

  update.binary<-function(mat,mat_b,mat_now,eps)
  {
    mat_p<-mat_b+mat_now
    mat_u<-matrix(0,nrow(mat),ncol(mat))
    idx1<-!is.na(mat) & mat==1
    idx0<-!is.na(mat) & mat==0
    index<-is.na(mat)
    arr<-exp(mat_p)
    mat_u[index]<-mat_now[index]
    mat_u[idx1]<-mat_now[idx1]+eps*epsilon.binary/(1.0+arr[idx1])
    mat_u[idx0]<-mat_now[idx0]-eps*epsilon.binary*arr[idx0]/(1.0+arr[idx0])
    return (mat_u)
  }

  stop.binary<-function(mat,mat_b,mat_now,mat_u)
  {
    index<-!is.na(mat)
    mn<-mat_b+mat_now
    mu<-mat_b+mat_u
    arn<-exp(mn)
    aru<-exp(mu)
    idx1<-!is.na(mat) & mat==1
    idx0<-!is.na(mat) & mat==0
    lgn<-sum(log(arn[idx1]/(1+arn[idx1])))+sum(log(1/(1+arn[idx0])))
    lgu<-sum(log(aru[idx1]/(1+aru[idx1])))+sum(log(1/(1+aru[idx0])))
    return (lgu-lgn)
  }

  LL.binary<-function(mat,mat_b,mat_u)
  {
    index<-!is.na(mat)
    mu<-mat_b+mat_u
    aru<-exp(mu)
    idx1<-!is.na(mat) & mat==1
    idx0<-!is.na(mat) & mat==0
    lgu<-sum(log(aru[idx1]/(1+aru[idx1])))+sum(log(1/(1+aru[idx0])))
    return (lgu)
  }

  LLmax.binary<-function(mat)
  {
    return (0)
  }

  LLmin.binary<-function(mat,mat_b)
  {
    index<-!is.na(mat)
    aru<-exp(mat_b)
    idx1<-!is.na(mat) & mat==1
    idx0<-!is.na(mat) & mat==0
    lgu<-sum(log(aru[idx1]/(1+aru[idx1])))+sum(log(1/(1+aru[idx0])))
    return (lgu)
  }

  binary_type_base <- function( data,dimension=2 ,name="test")
  {
    data<-check.binary(data,name)
    data_b<-base.binary(data)
    data_now<-matrix(0,nrow(data),ncol(data))
    data_u<-update.binary(data,data_b,data_now)
    data_u<-nuclear_approximation(data_u,dimension)
    while (T)
    {
      thr<-stop.binary(data,data_b,data_now,data_u)
      message(thr)
      if (thr<0.2)
      {
        break
      }
      data_now<-data_u
      data_u<-update.binary(data,data_b,data_now)
      data_u<-nuclear_approximation(data_u,dimension)
    }
    return (data_now)
  }

  #----------#
  # gaussian #
  #----------#

  epsilon.gaussian=0.5

  check.gaussian.row<-function(arr)
  {
    if (sum(!is.na(arr))==0)
    {
      return (F)
    }
    else
    {
      return (T)
    }
  }
  check.gaussian<-function(mat,name)
  {
    index<-array(T,nrow(mat))
    for(i in 1:nrow(mat))
    {
      if (sum(is.na(mat[i,])==ncol(mat)))
      {
        war<-paste("Warning: ",name,"'s ",as.character(i)," line is all NA. Delete this line",sep="")
        warning(war)
        index[i]<-F
      }
    }
    mat_c<-mat[index,]
    rownames(mat_c)<-rownames(mat)[index]
    colnames(mat_c)<-colnames(mat)
    return (mat_c)
  }

  base.gaussian.row<-function(arr)
  {
    idx<-!is.na(arr)
    return (mean(arr[idx]))
  }

  base.gaussian<-function(mat)
  {
    mat_b<-matrix(0,nrow(mat),ncol(mat))
    ar_b<-apply(mat,1,base.gaussian.row)
    mat_b[1:nrow(mat_b),]<-ar_b
    return (mat_b)
  }

  update.gaussian<-function(mat,mat_b,mat_now,eps)
  {
    mat_p<-mat_b+mat_now
    mat_u<-matrix(0,nrow(mat),ncol(mat))
    index<-!is.na(mat)
    mat_u[index]<-mat_now[index]+eps*epsilon.gaussian*(mat[index]-mat_p[index])
    index<-is.na(mat)
    mat_u[index]<-mat_now[index]
    return (mat_u)
  }

  stop.gaussian<-function(mat,mat_b,mat_now,mat_u)
  {
    index<-!is.na(mat)
    mn<-mat_b+mat_now
    mu<-mat_b+mat_u
    ren<-mat[index]-mn[index]
    reu<-mat[index]-mu[index]
    lgn<- -0.5*sum(ren*ren)
    lgu<- -0.5*sum(reu*reu)
    return (lgu-lgn)
  }

  LL.gaussian<-function(mat,mat_b,mat_u)
  {
    index<-!is.na(mat)
    mu<-mat_b+mat_u
    reu<-mat[index]-mu[index]
    lgu<- -0.5*sum(reu*reu)
    return (lgu)
  }

  LLmax.gaussian<-function(mat)
  {
    return (0.0)
  }

  LLmin.gaussian<-function(mat,mat_b)
  {
    index<-!is.na(mat)
    reu<-mat[index]-mat_b[index]
    lgu<- -0.5*sum(reu*reu)
    return (lgu)
  }

  gaussian_base<-function(data,dimension=2,name="test")
  {
    data<-check.gaussian(data,name)
    data_b<-base.gaussian(data)
    data_now<-matrix(0,nrow(data),ncol(data))
    data_u<-update.gaussian(data,data_b,data_now)
    data_u<-nuclear_approximation(data_u,dimension)
    while(T)
    {
      thr<-stop.gaussian(data,data_b,data_now,data_u)
      message(thr)
      if (thr<0.2)
      {
        break
      }
      data_now<-data_u
      data_u<-update.gaussian(data,data_b,data_now)
      data_u<-nuclear_approximation(data_u,dimension)
    }
    return (data_now)
  }

  #---------#
  # poisson #
  #---------#

  epsilon.poisson<-0.5

  check.poisson.row<-function(arr)
  {
    if (sum(!is.na(arr))==0)
    {
      return (F)
    }
    else
    {
      idx<-!is.na(arr)
      if (sum(arr[idx]<0)>0)
      {
        return (F)
      }
      else
      {
        return (T)
      }
    }
  }

  check.poisson<-function(mat,name)
  {
    w<-paste(name," is poisson type. Add 1 to all counts",sep="")
    message(w)
    index<-apply(mat,1,check.poisson.row)
    n<-sum(!index)
    if (n>0)
    {
      w<-paste("Warning: ",name," have ",as.character(n)," invalid lines!",sep="")
      warning(w)
    }
    mat_c<-mat[index,]+1
    rownames(mat_c)<-rownames(mat)[index]
    colnames(mat_c)<-colnames(mat)
    return (mat_c)
  }

  base.poisson.row<-function(arr)
  {
    idx<-!is.na(arr)
    m<-sum(log(arr[idx]))
    n<-sum(idx)
    return(m/n)
  }

  base.poisson<-function(mat)
  {
    mat_b<-matrix(0,nrow(mat),ncol(mat))
    ar_b<-apply(mat,1,base.poisson.row)
    mat_b[1:nrow(mat_b),]<-ar_b
    return (mat_b)
  }

  update.poisson<-function(mat,mat_b,mat_now,eps)
  {
    mat_p<-mat_b+mat_now
    mat_u<-matrix(0,nrow(mat),ncol(mat))
    index<-!is.na(mat)
    mat_u[index]<-mat_now[index]+eps*epsilon.poisson*(log(mat[index])-mat_p[index])
    index<-is.na(mat)
    mat_u[index]<-mat_now[index]
    return (mat_u)
  }

  stop.poisson<-function(mat,mat_b,mat_now,mat_u)
  {
    index<-!is.na(mat)
    mn<-mat_b+mat_now
    mu<-mat_b+mat_u
    lgn<-sum(mat[index]*mn[index]-exp(mn[index]))
    lgu<-sum(mat[index]*mu[index]-exp(mu[index]))
    return (lgu-lgn)
  }

  LL.poisson<-function(mat,mat_b,mat_u)
  {
    index<-!is.na(mat)
    mu<-mat_b+mat_u
    lgu<-sum(mat[index]*mu[index]-exp(mu[index]))
    return (lgu)
  }

  LLmax.poisson<-function(mat)
  {
    index<-!is.na(mat)
    lgu<-sum(mat[index]*log(mat[index])-mat[index])
    return (lgu)
  }

  LLmin.poisson<-function(mat,mat_b)
  {
    index<-!is.na(mat)
    lgu<-sum(mat[index]*mat_b[index]-exp(mat_b[index]))
    return (lgu)
  }

  poisson_type_base<-function(data,dimension=2,name="test")
  {
    data<-check.poisson(data,name)
    data_b<-base.poisson(data)
    data_now<-matrix(0,nrow(data),ncol(data))
    data_u<-update.poisson(data,data_b,data_now)
    data_u<-nuclear_approximation(data_u,dimension)
    while(T)
    {
      thr<-stop.poisson(data,data_b,data_now,data_u)
      message(thr)
      if (thr<0.2)
      {
        break
      }
      data_now<-data_u
      data_u<-update.poisson(data,data_b,data_now)
      data_u<-nuclear_approximation(data_u,dimension)
    }
    return (data_now)
  }

  #----#
  # na #
  #----#

  nuclear_approximation<-function(mat,dimension)
  {
    svd<-svd(mat,nu=0,nv=0)
    if (dimension<length(svd$d))
    {
      lambda<-svd$d[dimension+1]
      svd<-svd(mat,nu=dimension,nv=dimension)
      indexh<-svd$d>lambda
      indexm<-svd$d<lambda
      dia<-array(svd$d,length(svd$d))
      dia[indexh]<-dia[indexh]-lambda
      dia[indexm]<-0
      mat_low<-svd$u%*%diag(c(dia[1:dimension],0))[1:dimension,1:dimension]%*%t(svd$v)
    }
    else
    {
      mat_low<-mat
    }
    return (mat_low)
  }

  #------------#
  # LRAcluster #
  #------------#
  check.matrix.element<-function(x)
  {
    if (!is.matrix(x))
    {
      return (T)
    }
    else
    {
      return (F)
    }
  }

  ncol.element<-function(x)
  {
    return (ncol(x))
  }

  nrow.element<-function(x)
  {
    return (nrow(x))
  }

  check<-function(mat,type,name)
  {
    if (type=="binary")
    {
      return (check.binary(mat,name))
    }
    else if (type=="gaussian")
    {
      return (check.gaussian(mat,name))
    }
    else if (type=="poisson")
    {
      return (check.poisson(mat,name))
    }
    else
    {
      e<-paste("unknown type ",type,sep="")
      stop(e)
    }
  }

  eps<-0.0
  if (!is.list(data))
  {
    stop("the input data must be a list!")
  }
  c<-sapply(data,check.matrix.element)
  if (sum(c)>0)
  {
    stop("each element of input list must be a matrix!")
  }
  c<-sapply(data,ncol.element)
  if (length(levels(factor(c)))>1)
  {
    stop("each element of input list must have the same column number!")
  }
  if (length(data)!=length(types))
  {
    stop("data and types must be the same length!")
  }
  nSample<-c[1]
  loglmin<-0
  loglmax<-0
  loglu<-0.0
  nData<-length(data)
  for (i in 1:nData)
  {
    data[[i]]<-check(data[[i]],types[[i]],names[[i]])
  }
  nGeneArr<-sapply(data,nrow.element)
  nGene<-sum(nGeneArr)
  indexData<-list()
  k=1
  for(i in 1:nData)
  {
    indexData[[i]]<- (k):(k+nGeneArr[i]-1)
    k<-k+nGeneArr[i]
  }
  base<-matrix(0,nGene,nSample)
  now<-matrix(0,nGene,nSample)
  update<-matrix(0,nGene,nSample)
  thr<-array(0,nData)
  for (i in 1:nData)
  {
    if (types[[i]]=="binary")
    {
      base[indexData[[i]],]<-base.binary(data[[i]])
      loglmin<-loglmin+LLmin.binary(data[[i]],base[indexData[[i]],])
      loglmax<-loglmax+LLmax.binary(data[[i]])
    }
    else if (types[[i]]=="gaussian")
    {
      base[indexData[[i]],]<-base.gaussian(data[[i]])
      loglmin<-loglmin+LLmin.gaussian(data[[i]],base[indexData[[i]],])
      loglmax<-loglmax+LLmax.gaussian(data[[i]])
    }
    else if (types[[i]]=="poisson")
    {
      base[indexData[[i]],]<-base.poisson(data[[i]])
      loglmin<-loglmin+LLmin.poisson(data[[i]],base[indexData[[i]],])
      loglmax<-loglmax+LLmax.poisson(data[[i]])
    }
  }
  for (i in 1:nData)
  {
    if (types[[i]]=="binary")
    {
      update[indexData[[i]],]<-update.binary(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
    }
    else if (types[[i]]=="gaussian")
    {
      update[indexData[[i]],]<-update.gaussian(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
    }
    else if (types[[i]]=="poisson")
    {
      update[indexData[[i]],]<-update.poisson(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
    }
  }
  update<-nuclear_approximation(update,dimension)
  nIter<-0
  thres<-array(Inf,3)
  epsN<-array(Inf,2)
  while(T)
  {
    for (i in 1:nData)
    {
      if (types[[i]]=="binary")
      {
        thr[i]<-stop.binary(data[[i]],base[indexData[[i]],],now[indexData[[i]],],update[indexData[[i]],])
      }
      else if (types[[i]]=="gaussian")
      {
        thr[i]<-stop.gaussian(data[[i]],base[indexData[[i]],],now[indexData[[i]],],update[indexData[[i]],])
      }
      else if (types[[i]]=="poisson")
      {
        thr[i]<-stop.poisson(data[[i]],base[indexData[[i]],],now[indexData[[i]],],update[indexData[[i]],])
      }
    }
    nIter<-nIter+1
    thres[1]<-thres[2]
    thres[2]<-thres[3]
    thres[3]<-sum(thr)
    epsN[1]<-epsN[2]
    epsN[2]<-eps
    if (nIter>5)
    {
      if (runif(1)<thres[1]*thres[3]/(thres[2]*thres[2]+thres[1]*thres[3]))
      {
        eps<-epsN[1]+0.05*runif(1)-0.025
      }
      else
      {
        eps<-epsN[2]+0.05*runif(1)-0.025
      }
      if (eps< -0.7)
      {
        eps<- 0
        epsN<-c(0,0)
      }
      if (eps > 1.4)
      {
        eps<-0
        epsN<-c(0,0)
      }
    }
    if (sum(thr)<nData*0.2)
    {
      break
    }
    now<-update
    for (i in 1:nData)
    {
      if (types[[i]]=="binary")
      {
        update[indexData[[i]],]<-update.binary(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
      }
      else if (types[[i]]=="gaussian")
      {
        update[indexData[[i]],]<-update.gaussian(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
      }
      else if (types[[i]]=="poisson")
      {
        update[indexData[[i]],]<-update.poisson(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
      }
    }
    update<-nuclear_approximation(update,dimension)
  }
  for (i in 1:nData)
  {
    if (types[[i]]=="binary")
    {
      loglu<-loglu+LL.binary(data[[i]],base[indexData[[i]],],update[indexData[[i]],])
    }
    else if (types[[i]]=="gaussian")
    {
      loglu<-loglu+LL.gaussian(data[[i]],base[indexData[[i]],],update[indexData[[i]],])
    }
    else if (types[[i]]=="poisson")
    {
      loglu<-loglu+LL.poisson(data[[i]],base[indexData[[i]],],update[indexData[[i]],])
    }
  }
  sv<-svd(update,nu=0,nv=dimension)
  coordinate<-diag(c(sv$d[1:dimension],0))[1:dimension,1:dimension]%*%t(sv$v)
  colnames(coordinate)<-colnames(data[[1]])
  rownames(coordinate)<-paste("PC ",as.character(1:dimension),sep="")
  ratio<-(loglu-loglmin)/(loglmax-loglmin)
  return (list("coordinate"=coordinate,"potential"=ratio))
}
