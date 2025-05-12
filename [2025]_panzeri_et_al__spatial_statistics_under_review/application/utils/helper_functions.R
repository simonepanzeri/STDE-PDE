## COMPUTE THE KERNEL DENSITY ESTIMATE ON A LINEAR NETWORK USING THE HEAT EQUATION -----
density1Dkernel <- function(X,at=c("points","pixels"),leaveoneout=FALSE,
                            kernel=NULL,...){
  
  if(!inherits(X, "stlpp")) stop("X should an object of class stlpp")
  
  if(missing(at)) at <- "pixels"
  
  n <- npoints(X)
  
  Xs <- as.lpp.stlpp(X)
  
  Xt <- lpp(X=cbind(X$data$t,rep(0,n)), 
            L=linnet_interval(startp=X$time[1], endp=X$time[2]))
  
  bw <- bw.lppl(X=Xs)
  if(is.null(kernel)){
    IntEstL <- densityHeat(x=Xs,sigma=as.numeric(bw),...)
  } else {
    IntEstL <- densityEqualSplit(x=Xs,sigma=as.numeric(bw),kernel=kernel,...)
  } 
  
  bw <- bw.lppl(X=Xt)
  if(is.null(kernel)){
    IntEstT <- densityHeat(x=Xt,sigma=as.numeric(bw),...)
  } else {
    IntEstT <- densityEqualSplit(x=Xt,sigma=as.numeric(bw),kernel=kernel,...)
  } 
  
  tgrid <- IntEstT$xcol
  
  IntEstTv <- IntEstT$v[!is.na(IntEstT$v)]/n
  
  out <- lapply(X=(1:length(IntEstTv)), FUN=function(j){IntEstL*IntEstTv[j]}) 
  
  if(at=="points"){
    t <- X$data$t
    id <- findInterval(t,tgrid)
    out1 <- c()
    for (i in 1:n){
      out1[i] <- out[[id[i]]][Xs[i]]
    }
    out <- out1
  }
  
  if(at=="points") class(out) <- c("numeric","stlppint")
  if(at=="pixels") class(out) <- c("list","stlppint")
  
  attr(out,"lint") <- IntEstL
  attr(out,"tint") <- IntEstT$v[!is.na(IntEstT$v)]
  
  attr(out,"tgrid") <- tgrid
  attr(out,"time") <- X$data$t
  attr(out,"stlpp") <- X
  
  return(out)
  
}

linnet_interval <- function(startp=0, endp=1,...){
  
  Wt <- c(startp, endp)
  vertices <- ppp(x=c(Wt[1],Wt[2]), y=c(0,0), 
                  window=owin(Wt,c(-0.5,0.5)))
  m <- matrix(data=c(FALSE,TRUE,TRUE,FALSE), nrow=2, ncol=2, 
              byrow=TRUE)
  out <- linnet(vertices=vertices, m=m,...)
  
  return(out)
  
}