# GLOBAL VARIABLES: #

# PROCEDURE CODE: #
lr_var <- function(u,v,side){

  #if (side==1) {band=trunc(nobs^(0.2))} else{band=0}
  urprint <- 1
  
  
  tu <- nrow(u)
  p <- ncol(u)
  tv <- nrow(v)
  q <- ncol(v)
  e <- u
  te <- tu
  if (band == 0){
    eb <- as.matrix(e[1:(te-1),])
    ef <- as.matrix(e[2:te,])
    ae <- as.matrix(colSums(eb*ef)/colSums(eb^2))
    ee <- ef - eb*(matrix(1,nrow(eb),1)%*%t(ae))
    se <- as.matrix(colMeans(ee^2))
    ad <- sum((se/((1-ae)^2))^2)
    a1 <- 4*sum((ae*se/(((1-ae)^3)*(1+ae)))^2)/ad
    a2 <- 4*sum((ae*se/((1-ae)^4))^2)/ad
    if (kernel == 2){                   # Quadratic Spectral #
      bandw <- 1.3221*((a2*te)^.2)
    }   
    if (kernel == 1){                   # Parzen #
      bandw <- 2.6614*((a2*te)^.2)               
      if (bandw > (te-2)) bandw <- te-2
    }
    if (kernel == 3){                   # Bartlett #   
      bandw <- 1.1447*((a1*te)^.333)
      if (bandw > (te-2)) bandw <- te-2
    }
  }else{
    bandw <- band
  }
 
  # Estimate Covariances #
  if (kernel == 2){                             # Quadratic Spectral Kernel #
    tm <- te-1
    jb <- as.matrix(seq(1,tm,1)/bandw)
    jband <- jb*1.2*pi
    kern <- ((sin(jband)/jband - cos(jband))/(jband^2))*3
  }
  if (kernel == 1){                             # Parzen kernel #
    tm <- floor(bandw)   
    if (tm > 0){
      jb <- as.matrix(seq(1,tm,1)/bandw)
      kern <- (1 - (jb^2)*6 + (jb^3)*6)*(jb <= .5)
      kern <- kern + ((1-jb)^3)*(jb > .5)*2
    }
  }
  if (kernel == 3){                             # Bartlett kernel #
    tm <- floor(bandw)
    if (tm > 0) kern <- as.matrix(1 - seq(1,tm,1)/bandw)
  }
  
  if (kernel == 4){                             # flat-top kernel #
    tm <- floor(bandw)   
    if (tm > 0){
      jb <- as.matrix(seq(1,tm,1)/bandw)
      kern <- (1.1-abs(jb))*(jb >= .1 & jb < 1.1) + 1*(jb <= .1) + 0*(jb >= 1.1)  
     }
  }
  
  
  if (side==2 & tm>=1){
    lam <- matrix(0,p,q)
    lam2 <- matrix(0,p,q) 
    for (j in 1:(tm)){
      kj <- kern[j]
      lam <- lam + (t(as.matrix(u[1:(te-j),]))%*%as.matrix(v[(1+j):te,]))*kj/(te-j)
  #    lam2 <- lam2 + (t(as.matrix(u[(1+j):te,]))%*%as.matrix(v[1:(te-j),]))*kj/(te)
    }
    omega <- (t(u)%*%v/te + lam + t(lam))
  }
    if (side==2 & tm<=0){
     omega <- (t(u)%*%v/te)
  }
  if (side==1){
    lam <- matrix(0,p,q)
    for (j in 1:(tm)){
      kj <- kern[j]
      lam <- lam + (t(as.matrix(u[1:(te-j),]))%*%(as.matrix(v[(1+j):te,])))*kj/(te)
    }
    omega <- (t(u)%*%v/te + lam )
  }
  
  list(omega=omega,bandw=bandw)
}
