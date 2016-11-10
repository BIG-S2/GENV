ginv <- function (X, tol = sqrt(.Machine$double.eps)) 
{
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  Xsvd <- svd(X)
  if (is.complex(X)) 
    Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive)) 
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive)) 
    array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
                                               t(Xsvd$u[, Positive, drop = FALSE]))
}

auxinit <- function(M, U, MU, u, ng, n, p){
  
  tmp.MU <- eigen(MU)
  invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
  invMU2 <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / sqrt(tmp.MU$values), '*') %*% t(tmp.MU$vectors)
  
  startv <- function(a){
    out <- list(length=p)
    for (i in 1:p){
      out[[i]] <- t(a)%*% U[[i]] %*% a
    }
    sout <- Reduce("+", out)
    return(sout)
  }
  tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
  tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
  init <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]]) 
  
  obj1a <- c()
  for(i in 1:p){
    eig1 <- eigen(t(init) %*% M[[i]] %*% init)
    eig2 <- eigen(t(init) %*% invMU %*% init)
    obj1a[i] <- sum(log(eig1$values))*ng[i]/n + sum(log(eig2$values))/p
  }
  obj1 <- sum(obj1a)
  
  startv2 <- function(a){
    out <- list(length=p)
    for (i in 1:p){
      out[[i]] <- t(a)%*% invMU2 %*% tcrossprod(U[[i]], invMU2)  %*% a
    }
    sout <- Reduce("+", out)
    return(sout)
  }
  tmp2.MU <- apply(tmp.MU$vectors, 2, startv2)
  tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
  init.MU <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]])
  obj1a <- c()
  for(i in 1:p){
    e1 <- eigen(t(init.MU) %*% M[[i]] %*% init.MU)
    e2 <- eigen(t(init.MU) %*% invMU %*% init.MU)
    obj1a[i] <- sum(log(e1$values))*ng[i]/n + sum(log(e2$values))/p	
  }
  obj2 <- sum(obj1a)
  
  if (obj2 < obj1) {
    init <- init.MU
    obj1 <- obj2
  }
  
  Ma = Ua <- list(length=p)
  for (i in 1:p){
    Ma[[i]] <- M[[i]]*ng[i]/n
    Ua[[i]] <- U[[i]]/p
  }
  M1 <- Reduce("+", Ma)
  U1 <- Reduce("+", Ua)
  tmp.M <- eigen(M1)
  startv3 <- function(a) t(a) %*% U1 %*% a
  tmp2.M <- apply(tmp.M$vectors, 2, startv3)
  tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)	
  init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])
  e1 <- eigen(t(init.M) %*% M1 %*% init.M)
  e2 <- eigen(t(init.M) %*% invMU %*% init.M)
  obj3 <- sum(log(e1$values)) + sum(log(e2$values))	
  if (obj3 < obj1) {
    init <- init.M
    obj1 <- obj3
  }
  
  invM2 <- sweep(tmp.M$vectors, MARGIN = 2, 1 / sqrt(tmp.M$values), '*') %*% t(tmp.M$vectors)
  midmatrix <- invM2 %*% tcrossprod(U1, invM2) 
  startv4 <- function(a) t(a) %*% midmatrix %*% a
  tmp2.M <- apply(tmp.M$vectors, 2, startv4)
  tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
  init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])				
  e1 <- eigen(t(init.M) %*% M1 %*% init.M)
  e2 <- eigen(t(init.M) %*% invMU %*% init.M)			
  obj4 <- sum(log(e1$values)) + sum(log(e2$values))	
  if (obj4 < obj1) {
    init <- init.M
    obj1 <- obj4
  }
  
  return(list(init = init, obj1 = obj1, invMU = invMU))
  
}


GE <- function(A) {
  
  # Gaussian elimination, p must be less than or equal to n
  a <- dim(A)
  n <- a[1]
  p <- a[2]
  idx <- rep(0, p)
  res.idx <- 1:n
  
  i <- 1
  while (i <= p) {
    tmp <- max(abs(A[res.idx, i]))
    Stmp <- setdiff(which(abs(A[, i]) == tmp), idx)
    idx[i] <- Stmp[1]
    res.idx <- setdiff(res.idx, idx[i])
    for (j in 1:(n-i)) {
      A[res.idx[j], ] <- A[res.idx[j], ] - A[res.idx[j], i] / A[idx[i], i] * A[idx[i], ]
    }
    i <- i + 1			
  }
  c(idx, res.idx)
}


genv <- function(X, Y, Z, u){
  
  XX <- as.factor(X)
  Y <- as.matrix(Y)
  Z <- as.matrix(Z)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  c <- ncol(Z)
  p <- nlevels(XX)
  ncumx <- c()
  for (i in 1:p){
    ncumx[i] <- length(which(XX==i-1))
  }
  ncum <- cumsum(ncumx)
  ng <- diff(c(0,ncum))
  sortx <- sort(X, index.return=T)
  Xs <- sortx$x
  ind <- sortx$ix
  Ys <- Y[ind,]
  Zs <- Z[ind,]
  mY <- apply(Y, 2, mean)
  sigres = sigYcg = sigresc <- list(length=p)
  mYg <- matrix(rep(0, r*p), ncol=p)
  mZg <- matrix(rep(0, c*p), ncol=p)
  for (i in 1:p){
    if(i>1){
      yy <- Ys[(ncum[i-1]+1):ncum[i],]
      zz <- Zs[(ncum[i-1]+1):ncum[i],]
      sigres[[i]]<- cov(yy)*(ng[i]-1)/ng[i]
      mYg[,i] <- apply(yy, 2, mean)
      mZg[,i] <- apply(zz, 2, mean)
      Ycg <- scale(yy, center=T, scale=F)
      Zcg <- scale(zz, center=T, scale=F)
      QXC <- diag(ng[i]) - Zcg%*%solve(t(Zcg)%*%Zcg)%*%t(Zcg)
      sigYcg[[i]] <- t(Ycg)%*%Ycg
      sigresc[[i]] <- t(Ycg)%*%QXC%*%Ycg/ng[i]
    }else{
      yy <- Ys[1:ncum[i],]
      zz <- Zs[1:ncum[i],]
      sigres[[i]]<- cov(yy)*(ng[i]-1)/ng[i]
      mYg[,i] <- apply(yy, 2, mean)
      mZg[,i] <- apply(zz, 2, mean)
      Ycg <- scale(yy, center=T, scale=F)
      Zcg <- scale(zz, center=T, scale=F)
      QXC <- diag(ng[i]) - Zcg%*%solve(t(Zcg)%*%Zcg)%*%t(Zcg)
      sigYcg[[i]] <- t(Ycg)%*%Ycg
      sigresc[[i]] <- t(Ycg)%*%QXC%*%Ycg/ng[i]
    }
  }
  sigY <- cov(Y)*(n-1)/n
  eigtemY <- eigen(sigY)$values
  logDetSigY <- log(prod(eigtemY[eigtemY > 0]))
  invsigY <- solve(sigY)
  
  sigYc <- matrix(rep(0, r*r), ncol=r)
  for (i in 1:p){
    sigYc <- sigYc + sigYcg[[i]]/n
  }
  eigtemYc <- eigen(sigYc)$values
  logDetSigYc <- log(prod(eigtemYc[eigtemYc > 0]))
  invsigYc <- solve(sigYc)
  
  U = M <- list(length=p)
  for (i in 1:p){
    M[[i]] <- sigresc[[i]]
    U[[i]] <-  sigYc - M[[i]]
  }
  MU <- sigYc
  tmp <- genvMU(M, U, MU, u, n, ng, p)
  Gammahat <- tmp$Gammahat
  Gamma0hat <- tmp$Gamma0hat
  
  
  if(u==0){
    
    Yfit <- matrix(rep(0, n*r), ncol=r)
    Sigmahat <- sigYc
    etahat <- NULL
    Omegahat <- NULL
    Omega0hat <- sigYc
    betahat <- list(length=p)
    for (i in 1:p){
      betahat[[i]] <- matrix(rep(0, r*c), ncol=c)
    }
    muhat <- mYg
    bb <- rep(0, p)
    for (i in 1:p){
      if(i>1){
        yy <- Ys[(ncum[i-1]+1):ncum[i],]
        Ycg <- scale(yy, center=T, scale=F)
        bb[i] <- sum(diag(Ycg%*%solve(sigYc)%*%t(Ycg)))
      }else{
        yy <- Ys[1:ncum[i],]
        Ycg <- scale(yy, center=T, scale=F)
        bb[i] <- sum(diag(Ycg%*%solve(sigYc)%*%t(Ycg)))
      }
    }
    b <- sum(bb)
    loglik <- -n*r*log(2*pi)/2 - n*logDetSigYc/2 - b/2
    paranum <- p*r + p*(r-u)*(r-u+1)/2
    
  }else if (u==r){
    
    Sigmahat <- sigresc 
    aa <- rep(0, p)
    etahat <- list(length=p)
    for (i in 1:p){
      if(i>1){
        yy <- Ys[(ncum[i-1]+1):ncum[i],]
        zz <- Zs[(ncum[i-1]+1):ncum[i],]
        Ycg <- scale(yy, center=T, scale=F)
        Zcg <- scale(zz, center=T, scale=F)
        etahat[[i]] <- t(Ycg)%*%Zcg%*%solve(t(Zcg)%*%Zcg)
        QXC <- diag(ng[i]) - Zcg%*%solve(t(Zcg)%*%Zcg)%*%t(Zcg)
        aa[i] <- sum(diag(QXC%*%Ycg%*%solve(sigresc[[i]])%*%t(Ycg)%*%QXC))
      }else{
        yy <- Ys[1:ncum[i],]
        zz <- Zs[1:ncum[i],]
        Ycg <- scale(yy, center=T, scale=F)
        Zcg <- scale(zz, center=T, scale=F)
        etahat[[i]] <- t(Ycg)%*%Zcg%*%solve(t(Zcg)%*%Zcg)
        QXC <- diag(ng[i]) - Zcg%*%solve(t(Zcg)%*%Zcg)%*%t(Zcg)
        aa[i] <- sum(diag(QXC%*%Ycg%*%solve(sigresc[[i]])%*%t(Ycg)%*%QXC))
      }
    }
    a <- sum(aa)
    Omegahat <- sigresc
    Omega0hat <- NULL
    betahat <- etahat
    muhat <- mYg
    for (i in 1:p){
      muhat[,i]<- muhat[,i] - betahat[[i]]%*%mZg[,i]
    }
    Yfit <- matrix(rep(0, n*r), ncol=r)
    for (i in 1:p){
      if(i>1){
        Yfit[ind[(ncum[i-1]+1):ncum[i]],] <- rep(1,ng[i])%*%t(muhat[,i]) +
          Z[ind[(ncum[i-1]+1):ncum[i]],]%*% t(betahat[[i]]) 
      }else{
        Yfit[ind[1:ncum[1]],] <- rep(1,ng[1])%*%t(muhat[,1]) +
          Z[ind[1:ncum[1]],]%*% t(betahat[[1]]) 
      }
    }
    loglik <- -n*r*log(2*pi)/2 - a/2
    for (i in 1:p){
      eig <- eigen(sigresc[[i]])
      loglik <- loglik - ng[i]*sum(log(eig$values))/2
    }
    paranum <- p*r + p*u*c + u*(r-u) + p*u*(u+1)/2 + (r-u)*(r-u+1)/2
    
  }else{
    
    Omega0hat  <- t(Gamma0hat)%*%sigYc%*%Gamma0hat
    Sigmahat = Omegahat <- list(length=p)
    for(i in 1:p){
      Omegahat[[i]] <- t(Gammahat)%*%sigresc[[i]]%*%Gammahat
      Sigmahat[[i]] <- Gammahat%*%Omegahat[[i]]%*%t(Gammahat) +
        Gamma0hat%*%Omega0hat%*%t(Gamma0hat)
    }
    aa = bb <- rep(0, p)
    etahat = betahat <- list(length=p)
    for (i in 1:p){
      if(i>1){
        yy <- Ys[(ncum[i-1]+1):ncum[i],]
        zz <- Zs[(ncum[i-1]+1):ncum[i],]
        Ycg <- scale(yy, center=T, scale=F)
        Zcg <- scale(zz, center=T, scale=F)
        etahat[[i]] <- t(Gammahat)%*%t(Ycg)%*%Zcg%*%solve(t(Zcg)%*%Zcg)
        betahat[[i]] <- Gammahat%*%etahat[[i]]
        QXC <- diag(ng[i]) - Zcg%*%solve(t(Zcg)%*%Zcg)%*%t(Zcg)
        aa[i] <- sum(diag(QXC%*%Ycg%*%Gammahat%*%solve(t(Gammahat)%*%sigresc[[i]]%*%Gammahat)%*%t(Gammahat)%*%t(Ycg)%*%QXC))
        bb[i] <- sum(diag(Ycg%*%Gamma0hat%*%solve(t(Gamma0hat)%*%sigYc%*%Gamma0hat)%*%t(Gamma0hat)%*%t(Ycg)))
      }else{
        yy <- Ys[1:ncum[i],]
        zz <- Zs[1:ncum[i],]
        Ycg <- scale(yy, center=T, scale=F)
        Zcg <- scale(zz, center=T, scale=F)
        etahat[[i]] <- t(Gammahat)%*%t(Ycg)%*%Zcg%*%solve(t(Zcg)%*%Zcg)
        betahat[[i]] <- Gammahat%*%etahat[[i]]
        QXC <- diag(ng[i]) - Zcg%*%solve(t(Zcg)%*%Zcg)%*%t(Zcg)
        aa[i] <- sum(diag(QXC%*%Ycg%*%Gammahat%*%solve(t(Gammahat)%*%sigresc[[i]]%*%Gammahat)%*%t(Gammahat)%*%t(Ycg)%*%QXC))
        bb[i] <- sum(diag(Ycg%*%Gamma0hat%*%solve(t(Gamma0hat)%*%sigYc%*%Gamma0hat)%*%t(Gamma0hat)%*%t(Ycg)))
      }
    }
    muhat <- mYg
    for (i in 1:p){
      muhat[,i]<- muhat[,i] - betahat[[i]]%*%mZg[,i]
    }
    Yfit <- matrix(rep(0, n*r), ncol=r)
    for (i in 1:p){
      if(i>1){
        Yfit[ind[(ncum[i-1]+1):ncum[i]],] <- rep(1,ng[i])%*%t(muhat[,i]) +
          Z[ind[(ncum[i-1]+1):ncum[i]],]%*% t(betahat[[i]]) 
      }else{
        Yfit[ind[1:ncum[1]],] <- rep(1,ng[1])%*%t(muhat[,1]) +
          Z[ind[1:ncum[1]],]%*% t(betahat[[1]]) 
      }
    }
    a <- sum(aa)
    b <- sum(bb)
    eig2 <- eigen(t(Gamma0hat)%*%sigYc%*%Gamma0hat)
    r1 <- sum(log(eig2$values))
    loglik <- -n*r*log(2*pi)/2 - a/2 -b/2 - n*r1/2
    for (i in 1:p){
      eig <- eigen(t(Gammahat)%*%sigresc[[i]]%*%Gammahat)
      loglik <- loglik - ng[i]*sum(log(eig$values))/2
    }
    paranum <- p*r + p*u*c + u*(r-u) + p*u*(u+1)/2 + (r-u)*(r-u+1)/2
    
  }
  
  return(list(Yfit=Yfit, Gamma = Gammahat, Gamma0 = Gamma0hat,
              Sigma = Sigmahat,
              eta = etahat, Omega = Omegahat,
              Omega0 = Omega0hat, beta = betahat,
              mu = muhat, loglik = loglik,
              paranum = paranum, ng = ng))
}


genvasy <- function(X, Y, Z, u){
  
  XX <- as.factor(X)
  Y <- as.matrix(Y)
  Z <- as.matrix(Z)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  c <- ncol(Z)
  p <- nlevels(XX)
  ncumx <- c()
  for (i in 1:p){
    ncumx[i] <- length(which(XX==i-1))
  }
  ncum <- cumsum(ncumx)
  ng <- diff(c(0,ncum))
  sortx <- sort(X, index.return=T)
  Xs <- sortx$x
  ind <- sortx$ix
  Ys <- Y[ind,]
  Zs <- Z[ind,]
  
  m <- genv(X, Y, Z, u)
  Sigma <- m$Sigma
  Gamma <- m$Gamma
  Gamma0 <- m$Gamma0
  Omega <- m$Omega
  Omega0 <- m$Omega0
  eta <- m$eta
  
  xx <- list(length=p)
  for(i in 1:p){
    if(i>1){
      xx[[i]] <- cov(Zs[((ncum[i-1]+1):ncum[i]),])
    }else{
      xx[[i]] <- cov(Zs[(1:ncum[i]),])
    }
  }
  
  if(u==0){
    
    asySEmu <- NULL
    asySEbeta <- NULL
    
  }else if (u==r){
    
    asymu=asybeta <- list(length=p)
    semu <- matrix(rep(0, r*p), ncol=p)
    sebeta = asySEbeta <- list(length=p)
    for(i in 1:p){
      asymu[[i]] <- n*Sigma[[i]]/ng[[i]]
      asybeta[[i]] <- n*kronecker(solve(xx[[i]]), Omega[[i]])/ng[i]
      semu[,i] <- diag(asymu[[i]])
      sebeta[[i]] <- diag(asybeta[[i]])
      asySEbeta[[i]] <- matrix(sqrt(sebeta[[i]]), ncol=c) 
    }
    
    asySEmu <- sqrt(semu)
    
  }else{
    
    asymu=asybeta <- list(length=p)
    semu <- matrix(rep(0, r*p), ncol=p)
    sebeta = asySEbeta <- list(length=p)
    bb <- matrix(rep(0, u*(r-u)*u*(r-u)), ncol=u*(r-u))
    for(i in 1:p){
      b <- kronecker(eta[[i]]%*%xx[[i]]%*%t(eta[[i]]), solve(Omega0))+
        kronecker(Omega[[i]], solve(Omega0))+
        kronecker(solve(Omega[[i]]), Omega0)-
        2*kronecker(diag(u), diag(r-u))
      bb <- bb + ng[i]*b/n
    }
    for(i in 1:p){
      asymu[[i]] <- n*Sigma[[i]]/ng[[i]]
      aux <- kronecker(t(eta[[i]]), Gamma0)%*%solve(bb)%*%kronecker(eta[[i]], t(Gamma0))
      asybeta[[i]] <- n*kronecker(solve(xx[[i]]), Gamma%*%Omega[[i]]%*%t(Gamma))/ng[i] + aux
      semu[,i] <- diag(asymu[[i]])
      sebeta[[i]] <- diag(asybeta[[i]])
      asySEbeta[[i]] <- matrix(sqrt(sebeta[[i]]), ncol=c) 
    }
    
    asySEmu <- sqrt(semu)
    
  }
  
  return(list(asySEmu=asySEmu, asySEbeta=asySEbeta))
}


genvMU <- function(M, U, MU, u, n, ng, L){
  
  p <- L 
  dimM <- dim(M[[1]])
  dimU <- dim(U[[1]])
  r <- dimM[1]
  
  if(u==0){
    Gammahat <- NULL
    Gamma0hat <- diag(r)
  }else if (u==r){
    Gammahat <- diag(r)
    Gamma0hat <- NULL
  }else if (u==r-1){
    
    maxiter = 100
    ftol = 1e-3
    initout <- auxinit(M, U, MU, u, ng, n, p)
    init <- initout$init
    obj1 <- initout$obj1
    invMU <- initout$invMU 
    
    GEidx <- GE(init)
    Ginit <- init %*% solve(init[GEidx[1:u], ])
    j <- GEidx[r]
    
    fobj <- function(x) {
      res <- -2 * log(1 + sum(x^2)) 
      for(i in 1:p){
        res <- res + auxf1(M[[i]], U[[i]], u, n, ng[i], p, init, x, r)
      }
      return(res)
    }
    
    gobj <- function(x) {
      res <- -4 * x %*% solve(1 + sum(x^2))
      for(i in 1:p){
        res <- res + auxg1(M[[i]], U[[i]], u, n, ng[i], p, init, x, r)
      }
      return(res)
    }
    
    i <- 1
    while (i < maxiter) {
      
      res <- optim(Ginit[j,], fobj, gobj, method = "BFGS")
      Ginit[j, ] <- res$par
      a <- qr(Ginit)
      Gammahat <- qr.Q(a)
      obj5a <- c()
      for(i in 1:p){
        e1 <- eigen(t(Gammahat) %*% M[[i]] %*% Gammahat)
        e2 <- eigen(t(Gammahat) %*% invMU %*% Gammahat)		
        obj5a[i] <- sum(log(e1$values))*ng[i]/n + sum(log(e2$values))/p
      }
      obj5 <- sum(obj5a)
      
      if (abs(obj1 - obj5) < ftol * abs(obj1)) {
        break
      } else {
        obj1 <- obj5
        i <- i + 1
      }
    }
    Gamma0hat <- qr.Q(a, complete = TRUE)[, (u+1):r, drop = FALSE]
    
  }else{
    
    maxiter <- 100
    ftol <- 1e-3
    
    initout <- auxinit(M, U, MU, u, ng, n, p)
    init <- initout$init
    obj1 <- initout$obj1
    
    GEidx <- GE(init)
    Ginit <- init %*% solve(init[GEidx[1:u], ])
    
    GUG = GVG <- list(length=p)
    for (k in 1:p){
      MU <- M[[k]] + U[[k]]
      tmp.MU <- eigen(MU)
      invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
      GUG[[k]] <- crossprod(Ginit, (M[[k]] %*% Ginit))	
      GVG[[k]] <- crossprod(Ginit, (invMU %*% Ginit))	
    }
    t4 <- crossprod(Ginit[GEidx[(u+1):r],], Ginit[GEidx[(u+1):r], ]) + diag(u)
    i <- 1
    while (i < maxiter) {
      
      for (j in GEidx[(u+1):r]) {
        g <- as.matrix(Ginit[j, ])
        t4 <- t4 - tcrossprod(g, g)
        invt4 <- chol2inv(chol(t4))	
        
        t2 = t3 = GUGt2 = GVGt2 = invc1 = invc2 <- list(length=p)
        for(k in 1:p){
          Maux <- M[[k]]
          MU <- M[[k]] + U[[k]]
          tmp.MU <- eigen(MU)
          invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
          t2[[k]] <- crossprod(Ginit[-j, ], as.matrix(Maux[-j, j])) / Maux[j, j]
          t3[[k]] <- crossprod(Ginit[-j, ], as.matrix(invMU[-j, j])) / invMU[j, j]
          
          GUGt2[[k]] <- g + t2[[k]]
          GUG[[k]] <- GUG[[k]] - tcrossprod(GUGt2[[k]], GUGt2[[k]]) * Maux[j, j]
          
          GVGt2[[k]] <- g + t3[[k]]
          GVG[[k]] <- GVG[[k]] - tcrossprod(GVGt2[[k]], GVGt2[[k]]) * invMU[j, j] 
          
          invc1[[k]] <- ginv(GUG[[k]]) #chol2inv(chol(GUG[[k]]))
          invc2[[k]] <- ginv(GVG[[k]]) #chol2inv(chol(GVG[[k]]))
        }
        fobj <- function(x) {
          res <- -2 * log(1 + x %*% invt4 %*% x)
          for(k in 1:p){
            res <- res + auxf2(M[[k]], U[[k]], t2[[k]], t3[[k]], 
                               invc1[[k]], invc2[[k]], ng[k], n, p, x, j)
          }
          return(res)
        }
        
        gobj <- function(x) {
          res <- -4	* invt4 %*% x / as.numeric(1 + x %*% invt4 %*% x)
          for(k in 1:p){
            res <- res + auxg2(M[[k]], U[[k]], t2[[k]], t3[[k]], 
                               invc1[[k]], invc2[[k]], ng[k], n, p, x, j)
          }
          return(res)
        }
        
        res <- optim(Ginit[j,], fobj, gobj, method = "BFGS")
        Ginit[j, ] <- res$par
        g <- as.matrix(Ginit[j, ])
        t4 <- t4 + tcrossprod(g, g)
        
        for(k in 1:p){
          Maux <- M[[k]]
          MU <- M[[k]] + U[[k]]
          tmp.MU <- eigen(MU)
          invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
          GUGt2[[k]] <- g + t2[[k]]
          GUG[[k]] <- GUG[[k]] + tcrossprod(GUGt2[[k]], GUGt2[[k]]) * Maux[j, j]
          
          GVGt2[[k]] <- g + t3[[k]]
          GVG[[k]] <- GVG[[k]] + tcrossprod(GVGt2[[k]], GVGt2[[k]]) * invMU[j, j] 
        }
        
      }
      a <- qr(Ginit)
      Gammahat <- qr.Q(a)
      obj5a <- c()
      for(i in 1:p){
        e1 <- eigen(t(Gammahat) %*% M[[i]] %*% Gammahat)
        e2 <- eigen(t(Gammahat) %*% invMU %*% Gammahat)		
        obj5a[i] <- sum(log(e1$values))*ng[i]/n + sum(log(e2$values))/p
      }
      obj5 <- sum(obj5a)
      if (abs(obj1 - obj5) < ftol * abs(obj1)) {
        break
      } else {
        obj1 <- obj5
        i <- i + 1
      }
    }
    Gamma0hat <- qr.Q(a, complete = TRUE)[, (u+1):r]
  }
  
  return(list(Gammahat = Gammahat, Gamma0hat = Gamma0hat))
  
}


u_genv <- function(X, Y, Z) {
  
  XX <- as.factor(X)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  c <- ncol(Z)
  p <- nlevels(XX)
  
  paranum = loglik.seq <- c()
  for (u in 0:r){
    loglik.seq[u+1] <- genv(X, Y, Z, u)$loglik
    paranum[u+1] <- p*r + p*u*c + u*(r-u) + p*u*(u+1)/2 + (r-u)*(r-u+1)/2
  }
  bic.seq <- -2 * loglik.seq + log(n) * paranum
  u.bic <- which.min(bic.seq) - 1
  
  return(list(u.bic = u.bic, bic.seq = bic.seq))
}


boot_genv <- function(X, Y, Z, u, B) {
  
  XX <- as.factor(X)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  c <- ncol(Z)
  p <- nlevels(XX)
  ncumx <- c()
  for (i in 1:p){
    ncumx[i] <- length(which(XX==i-1))
  }
  ncum <- cumsum(ncumx)
  ng <- diff(c(0,ncum))
  sortx <- sort(X, index.return=T)
  Xs <- sortx$x
  ind <- sortx$ix
  
  fit <- genv(X, Y, Z, u)
  Yfit <- fit$Yfit
  res <- Y - Yfit
  
  bootgenv <- function(i) {
    out <- list(length=p)
    res.boot <- matrix(rep(0, n*r), ncol=r)
    for (j in 1:p){
      if(j>1){
        res.boot[ind[(ncum[j-1]+1):ncum[j]],] <- res[sample(ind[(ncum[j-1]+1):ncum[j]], ng[j], replace=T),]
      }else{
        res.boot[ind[1:ncum[1]],] <- res[sample(ind[1:ncum[1]], ng[1], replace=T),]
      }
    }
    Y.boot <- Yfit + res.boot
    for(k in 1:p){
      out[[k]] <- genv(X, Y.boot, Z, u)$beta[[k]]
    }
    return(out)
  }
  
  bootsebeta <- list(length=p)
  for (k in 1:p){
    out1 <- lapply(1:B, function(i) bootgenv(i)[[k]])
    bootbeta <- matrix(unlist(out1), nrow = B, byrow = TRUE)
    bootsebeta[[k]] <- matrix(apply(bootbeta, 2, sd), nrow = r)
  }
  
  return(list(bootsebeta = bootsebeta))
  
}


cv_genv <- function(X, Y, Z, u, m, nperm){
  
  XX <- as.factor(X)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  c <- ncol(Z)
  p <- nlevels(XX)
  
  prederr <- matrix(rep(0, m*nperm), ncol = nperm)
  PE <- rep(0, nperm)
  for (i in 1:nperm)	{
    id <- sample(n, n)
    Xn <- X[id]
    Yn <- Y[id, ]
    Zn <- Z[id, ]
    for (j in 1:m) {
      id.test <- (floor((j - 1) * n / m) + 1) : ceiling(j * n / m)
      id.train <- setdiff(1:n, id.test)
      X.train <- Xn[id.train]
      Y.train <- Yn[id.train, ]
      Z.train <- Zn[id.train, ]
      X.test <- Xn[id.test]
      Y.test <- Yn[id.test, ]
      Z.test <- Zn[id.test, ]
      n.test <- length(id.test)
      
      fit <- genv(X.train, Y.train, Z.train, u)
      betahat <- fit$beta
      muhat <- fit$mu
      
      testn <- length(id.test)
      traceres <- 0
      for(l in 1:testn){
        iG <- X.test[l]+1
        resi <- Y.test[l,] - t(muhat[,iG]) - Z.test[l,]%*%t(betahat[[iG]])
        traceres <- traceres + resi%*%t(resi)
      }
      
      prederr[j, i] <- traceres
    }
    PE[i] <- sqrt(sum(prederr[,i])/n)
  }
  
  out <- mean(PE)
  return(out)
  
}

auxf1 <- function(M1, U1, u, n, ng, p, init, x, r){
  M <- M1
  U <- U1
  MU <- M + U
  tmp.MU <- eigen(MU)
  invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
  
  GEidx <- GE(init)
  Ginit <- init %*% solve(init[GEidx[1:u], ])		
  
  j <- GEidx[r]
  g <- as.matrix(Ginit[j, ])
  t2 <- crossprod(Ginit[-j, ], as.matrix(M[-j, j])) / M[j, j]
  t3 <- crossprod(Ginit[-j, ], as.matrix(invMU[-j, j])) / invMU[j, j]
  
  GUGt2 <- g + t2
  GUG <- crossprod(Ginit, (M %*% Ginit)) - tcrossprod(GUGt2, GUGt2) * M[j, j]
  
  GVGt2 <- g + t3
  GVG <- crossprod(Ginit, (invMU %*% Ginit)) - tcrossprod(GVGt2, GVGt2) * invMU[j, j] 
  
  invc1 <- chol2inv(chol(GUG))
  invc2 <- chol2inv(chol(GVG))
  
  tmp2 <- x + t2
  tmp3 <- x + t3
  T2 <- invc1 %*% tmp2	
  T3 <- invc2 %*% tmp3
  out <- ng*log(1 + M[j, j] * crossprod(tmp2, T2))/n + 
    log(1 + invMU[j, j] * crossprod(tmp3, T3))/p
  return(out)
}


auxf2 <- function(M1, U1, t2, t3, invc1, invc2, ng, n, p, x, j){
  M <- M1
  U <- U1
  MU <- M + U
  tmp.MU <- eigen(MU)
  invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
  
  tmp2 <- x + t2
  tmp3 <- x + t3
  
  T2 <- invc1 %*% tmp2	
  T3 <- invc2 %*% tmp3
  out <- ng*log(1 + M[j, j] * crossprod(tmp2, T2))/n + 
    log(1 + invMU[j, j] * crossprod(tmp3, T3))/p
  return(out)
  
}


auxg1 <- function(M1, U1, u, n, ng, p, init, x, r){
  M <- M1
  U <- U1
  MU <- M + U
  tmp.MU <- eigen(MU)
  invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
  
  GEidx <- GE(init)
  Ginit <- init %*% solve(init[GEidx[1:u], ])		
  
  j <- GEidx[r]
  g <- as.matrix(Ginit[j, ])
  t2 <- crossprod(Ginit[-j, ], as.matrix(M[-j, j])) / M[j, j]
  t3 <- crossprod(Ginit[-j, ], as.matrix(invMU[-j, j])) / invMU[j, j]
  
  GUGt2 <- g + t2
  GUG <- crossprod(Ginit, (M %*% Ginit)) - tcrossprod(GUGt2, GUGt2) * M[j, j]
  
  GVGt2 <- g + t3
  GVG <- crossprod(Ginit, (invMU %*% Ginit)) - tcrossprod(GVGt2, GVGt2) * invMU[j, j] 
  
  invC1 <- chol2inv(chol(GUG))
  invC2 <- chol2inv(chol(GVG))
  
  tmp2 <- x + t2
  tmp3 <- x + t3
  T2 <- invC1 %*% tmp2	
  T3 <- invC2 %*% tmp3
  out <-  2 * ng * T2 / (n*
                           as.numeric(1 / M[j, j] + crossprod(tmp2, T2))) + 2 * T3 /(p* 
                                                                                       as.numeric(1 / invMU[j, j] + crossprod(tmp3, T3)))
  return(out)
}


auxg2 <- function(M1, U1, t2, t3, invc1, invc2, ng, n, p, x, j){
  M <- M1
  U <- U1
  MU <- M + U
  tmp.MU <- eigen(MU)
  invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
  
  tmp2 <- x + t2
  tmp3 <- x + t3
  
  T2 <- invc1 %*% tmp2	
  T3 <- invc2 %*% tmp3
  out <-  2 * T2 * ng / (n*as.numeric(1 / M[j, j] + crossprod(tmp2, T2))) +
    2 * T3 / (p*as.numeric(1 / invMU[j, j] + crossprod(tmp3, T3)))	
  return(out)
  
}
