#' Function simulateCentersCors
#' 
#' Simulate a set of centers drawn from a bivariate distribution, with a given average
#' correlation and a given standard deviation on that correlation (on the Fisher scale).
#'
#' @param cc The mean correlation
#' @param rho Standard deviation on the correlation
#' @param Np Number of patients for each center. Number of centers drawn is given by the length of Np.
#' @return A matrix with 2 columns with the bivariate normal data
#' @export
simulateCentersCors = function(cc, rho, Np)
{ ctr = tanh(atanh(cc)+rnorm(length(Np))*rho);
  do.call(rbind, lapply(seq_along(ctr), function(i) mvrnorm(Np[i], mu=c(0,0), Sigma=matrix(c(1,ctr[i],ctr[i],1), nrow=2))));
} 

#' Function calcPars
#' 
#' Calculate correlations, margin s.d. on each center in a dataset.
#'
#' @param x A matrix with two features on which to calculate the correlations (etc)
#' @param ctr The center to which each observation belong
#' @return A list of a list with the following components
#'   \item{cc}{Correlations value per center}
#'   \item{N}{Number of non-NA entries in the dataset per center}
#'   \item{sd}{Marginal sds per center (one per line)}
#' @export
#' @seealso \code{\link{testCorFisher}}, \code{\link{testCorMargin}}
calcPars = function(x, ctr)
{ x = scale(x[complete.cases(x),]);
  res = tapply(1:nrow(x), ctr, function(i)
  {  d2 = x[i, ];
     if (nrow(d2) < 4) { return(NULL); }
     cv = cov(d2);
     if (any(cv==0)) { return(NULL); }
     cr = cov2cor(cv);
     return(list(x=cr[1,2], cv=cv, N=nrow(d2)));
  } );
  res = res[!sapply(res, is.null)];
  
  cc = sapply(res, function(i) i$x);
  sd = do.call(rbind, lapply(res, function(i) diag(i$cv)));
  N = sapply(res, function(i) i$N);
  return(list(cc=cc, N=N, sd=sd));
}

#' Function testCorFisher
#' 
#' Apply the Fisher test on a set of correlations.
#'
#' @param corrs A vector with the correlations per center
#' @param Ni Number of data points used to calculate the correlations, per center
#' @param retAll if FALSE, return only the p-values per center
#' @return If retAll is TRUE, A list of a list with the following components
#'   \item{p}{P-values per center. P-values signs indicate whether centers are higher or lower than expected. }
#'   \item{mu}{Mean correlation}
#'   \item{raneff}{Standard deviation on the correlation)}
#' @export
#' @seealso \code{\link{calcPars}}, \code{\link{testCorMargin}}
testCorFisher = function(corrs, Ni, retAll=FALSE)
{ likel = function(pars, cc, NiVar)
  { mu = atanh(pars[1]); raneff = exp(pars[2]); 
    -sum(dnorm(cc, mean=mu, sd=sqrt(raneff+NiVar), log = TRUE));
  }
  
  fishercorr<-atanh(corrs);
  prod2<-(Ni-3)*fishercorr;
  #mu = sum(prod2)/sum(Ni-3);
  mu = median(corrs);
  par = optim(c(tanh(mu), log(.05)), likel, cc=fishercorr, NiVar = 1/(Ni-3))$par;
  mu = atanh(par[1]); raneff = exp(par[2]);
  zscore = (fishercorr-mu)/(sqrt(raneff+1/(Ni-3)));
  p1 = pnorm(zscore); p2 = pnorm(zscore, lower.tail=FALSE);
  p = 2*pmin(p1,p2)*ifelse(p1>p2, +1, -1);
  if (!retAll) { return(p); }
  return(list(p=p, mu=mu, raneff=raneff));
}

#' Function testCorMargin
#' 
#' Apply the Margin restricted test on a set of correlations.
#'
#' @param corrs A vector with the correlations per center
#' @param Ni Number of data points used to calculate the correlations, per center
#' @param sds marginal standard deviations per center (one column for each margin)
#' @param retAll if FALSE, return only the p-values per center
#' @return If retAll is TRUE, A list of a list with the following components
#'   \item{p}{P-values per center. P-values signs indicate whether centers are higher or lower than expected. }
#'   \item{mu}{Mean correlation}
#'   \item{raneff}{Standard deviation on the correlation)}
#' @export
#' @seealso \code{\link{calcPars}}, \code{\link{testCorFisher}}
testCorMargin = function(corrs, Ni, sds, retAll=FALSE)
{ mrho = median(atanh(corrs));
  msig = log(mad(atanh(corrs)));
  bp = optim(c(mrho,msig), toOptv3, cc=corrs, Ni=Ni, sds=sds, control=list(reltol=0.00001));
  
  bp$par[2] = exp(bp$par[2]);
  
  p = do.call(rbind, lapply(seq_along(corrs), function(i) 
  { c(ctrPv3(corrs[i], Ni[i], sds[i,], bp$par[1], bp$par[2], lower.tail = TRUE),
      ctrPv3(corrs[i], Ni[i], sds[i,], bp$par[1], bp$par[2], lower.tail = FALSE))
  }))
  p = 2*pmin(p[,1],p[,2])*ifelse(p[,1]>p[,2], +1, -1);
  if (!retAll) { return(p); }
  return(list(p=p, mu=bp[1], raneff=bp[2]));
} 

toOptv3 = function(p, cc, Ni, sds)
{ mu = p[1]; sigma=exp(p[2]);
  res = rep(NA, length(cc));
  for (i in 1:length(cc)) {
    res[i] = suppressWarnings(log(ctrLikelv3(cc[i], Ni[i], sds[i,], mu, sigma)));
  }
  mr = median(res);
  res[res < mr-20] = mr-20;
  res = sum(res, na.rm=TRUE);
  return(-res);
}

ctrLikelv3 = function(x, N, cv, mu, sigma)
{ 
    if(x>1-1e-10){x=1-1e-10}
    if(x<(-1+1e-10)){x=-1+1e-10}
    
    baseM = sqrt(cv[2]/cv[1]);
    baseS = sqrt(N-1) / (N * sqrt(cv[1]));
    toI = function(r)
        {r2 = tanh(r);
         theta  = r2 / sqrt(1 - r2^2) * sqrt(cv[2]) * sqrt(N)
         xx = sqrt((N - 2) * x ^ 2 / (1 - x ^ 2));
         if(x<0){xx=-xx}
         if(N<100){ dnorm(r, mean=mu, sd=sigma) * suppressWarnings(dt(xx,N-2,theta));}
         else{dnorm(r, mean=mu, sd=sigma) * dnorm(xx,mean=theta,sd=1)}
     }

    if (baseS < sigma) { b = atanh(x)/baseM + 12*c(-baseS, baseS)/baseM; }
    else {  b = mu + 6*c(-sigma, +sigma); }
    #b1[b1>1] = 1; b1[b1< -1] =-1;
    #b1 = atanh(b1);
    #b2 = mu + 6*c(-sigma, +sigma);
    #b = c(max(b1[1],b2[1]), min(b1[2], b2[2]));
    
    while (toI(b[1]) > 1e-20) { b[1] = b[1] - .1; }
    while (toI(b[2]) > 1e-20) { b[2] = b[2] + .1; }
  
    tmp = integrate(toI, b[1], b[2], subdivisions=1e2, stop.on.error=FALSE);
    if (tmp$message != "OK" || tmp$value<1e-30)
    { i = seq(b[1], b[2], len=1e3);
      j = c(0,toI(i),0);
      i = c(b[1], i, b[2]);
      if (j[2]==0) { b[1] = i[which(j>0)[1]-1]; };
      if (is.na(b[1])) { return(1e-38); }
      if (j[length(j)-1]==0) { b[2] = i[which(j==0 & i>b[1])[1]]; }
      tmp = integrate(toI, b[1], b[2], subdivisions=1e2);
    }
    return(tmp$value);
}

ctrPv3 = function(x, N, cv, mu, sigma, lower.tail = TRUE)
{   baseM = sqrt(cv[2]/cv[1]);
    baseS = sqrt(N-1) / (N * sqrt(cv[1]));
    if(x>1-1e-10){x=1-1e-10}
    if(x<(-1+1e-10)){x=-1+1e-10}
    
    toI = function(r,lower.tail) {
        r2 = tanh(r);
        
        r2[r2==1]=1-1e-10;
        r2[r2==-1]=-1+1e-10;
        
        theta  = r2 / sqrt(1 - r2 ^ 2) * sqrt(cv[2]) * sqrt(N);
        xx = sqrt((N - 2) * x ^ 2 / (1 - x ^ 2));
        if (x < 0) {
            xx = -xx
        }
        
        pt= suppressWarnings(pt(xx,N - 2,theta,lower.tail = lower.tail));
        p= dnorm(r, mean = mu, sd = sigma) *pt;
    }
    b = mu + 6 * c(-sigma,+sigma);
    if (!lower.tail) {
        tmp = x - 12 * baseS / baseM;
        if (tmp > -1) {
            b[1] = max(b[1], atanh(tmp));
        }
    } else {
        tmp = x + 12 * baseS / baseM;
        if (tmp < 1) {
            b[2] = min(b[2], atanh(x + 12 * baseS / baseM));
        }
    }
    while (toI(b[1], lower.tail) > 1e-20) {
        b[1] = b[1] - .05;
    }
    while (toI(b[2], lower.tail) > 1e-20) {
        b[2] = b[2] + .05;
    }
    
    tmp = integrate(toI, b[1], b[2], lower.tail=lower.tail, subdivisions=1e2, stop.on.error=FALSE);
    if (tmp$message != "OK")
    { i = seq(b[1], b[2], len=1e3);
      j = c(0,toI(i, lower.tail),0);
      i = c(b[1], i, b[2]);
      if (j[2]==0) { b[1] = i[which(j>0)[1]-1]; };
      if (is.na(b[1])) { return(0); }
      if (j[length(j)-1]==0) { b[2] = i[which(j==0 & i>b[1])[1]]; }
      tmp = integrate(toI, b[1], b[2], lower.tail=lower.tail, subdivisions=1e3);
    }
    return(tmp$value);
}