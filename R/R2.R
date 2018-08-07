R2 <- function(poi,V,poiType,pval = T) {
  #
  #   R2(poi,V,poiType)
  #
  # Args:
  # - V is a p*k matrix, where the rows corresponds to the samples
  # - poi is a matrix p*l, representing the phenotypes of interest
  # - poiType (1*l) is the types of poi: 'continuous' (then a linear
  # regression is used) or 'categorical' (then the mean by class is used)
  #
  # Outputs:
  # - R2(l), higher R^2 value between a column of V and poi(l)
  # - idxCorr(l), index of the column of V giving the higher R^2 value (if many,
  # takes the first one)
  # - allR2(k,l),  R2 value for column k of V with poi l
  #
  #    IF pval =TRUE, return also:   #
  # - pv(l) smaller p-value association between a column of V and poi(l)
  # - idxcorr2(l) index of the column of V giving the smaller p-value (if many,
  #                                                                          # takes the first one)
  # - allpv(k,l),  p-value for column k of V with poi l
  #
  # if missing information in poi, remove the corresponding samples in the R2 computation


  if (is.vector(V) ){ V = matrix(V,ncol=1)}
  if (is.vector(poi)){poi = matrix(poi,nrow =length(poi))}

  p = nrow(V)   # number of samples
  k = ncol(V)    # number of components
  l = length(poiType)  # number of cf/poi to test
  if (is.null(l)){stop("POI type(s) neeeded")}
  p2 = nrow(poi)
  l2 = ncol(poi)

  if( l2 != l){ # checking poi and poiType dimensions compatiblity
    if (p2 == l){  # if poi is transposed (l*p)
      poi = t(poi)
      warning("Transposing poi to match poiType dimension")
      p2 = nrow(poi)
    } else {
      print(poi)
      print(poiType)
      stop("poi dimensions doesn't match poiType dimension")
    }
  }


  if (p != p2){ # checking poi and V dimensions compatiblity
    if (p2 == k){
      warnings("Transposing V to match poi dimension")
      V =t(V)
      k = p
      p = p2
    } else {
      stop("poi and V dimensions incompatible")
    }
  }






  R2 = rep(-1,l)
  names(R2) = colnames(poi)
  idxcorr = R2
  R2_tmp <- matrix(rep(-1,k*l),k,l,dimnames=list(colnames(V),colnames(poi)))    # r2_tmp(k,l) hold the R2 value for column k of V with poi l

  if (pval){
    pv = R2
    idxcorr2 = R2
    pv_tmp <- R2_tmp   # r2_tmp(k,l) hold the R2 value for column k of V with poi l
  }

  for (cmpt in 1:k){    # for each column of V
    cmpt2an <- V[,cmpt]
    for (ipoi in 1:l){
      idx_finite = is.finite(as.factor(poi[,ipoi]))
      poi2an = poi[idx_finite,ipoi]
      cmpt2an_finite=cmpt2an[idx_finite]
      if (poiType[ipoi] == "continuous") {  # estimation by linear regression
        coefs <- coef(lm(cmpt2an_finite~as.numeric(poi2an)))
        cmpt2an_est <- coefs[2]*as.numeric(poi2an)+coefs[1]
        nc <- 2;
      } else if (poiType[ipoi]=="categorical"){  # estimation by classe mean
        classes <- unique(poi2an)
        nc <- length(classes)
        cmpt2an_est <- rep(NA,length(cmpt2an_finite))
        for (icl in 1:length(classes) ){
          idxClasse <- which(poi2an==classes[icl])
          cmpt2an_est[idxClasse] <- mean(cmpt2an_finite[idxClasse])
        }
      } else {
        stop("Incorrect poiType. Select 'continuous' or 'categorical'. ")
      }
      sse <- sum((cmpt2an_finite-cmpt2an_est)^2)
      sst <- sum((cmpt2an_finite-mean(cmpt2an_finite))^2)
      R2_tmp[cmpt,ipoi] <-  1 - sse/sst
      if (pval){
        F <- ((sst-sse)/(nc-1))/(sse/(p-nc))
        pv_tmp[cmpt,ipoi] = 1-pf(F,nc-1,p-nc);
        if (!is.finite(pv_tmp[cmpt,ipoi])) {
          warning(paste("Non finite p-value for component ",cmpt," (pv=",pv_tmp[cmpt,ipoi],", F=",F,"), assigning NA", sep=""))
          pv_tmp[cmpt,ipoi] <- NA
        }
      }


    }
  }

  for (ipoi in 1:l){
    if (pval){
      pv[ipoi] <- min(pv_tmp[,ipoi])
      idxcorr2[ipoi] <- which(pv_tmp[,ipoi] == pv[ipoi])[1]   # if more than one component gives the best R2, takes the first one
    }
    R2[ipoi] <- max(R2_tmp[,ipoi])
    idxcorr[ipoi] <- which(R2_tmp[,ipoi] == R2[ipoi])[1]   # if more than one component gives the best R2, takes the first one
  }

  if (pval){
    return(list(R2=R2,idxcorr=idxcorr,allR2 = R2_tmp, pv=pv,idxcorr2=idxcorr2,allpv = pv_tmp))
  } else {
    return(list(R2=R2,idxcorr=idxcorr,allR2 = R2_tmp))
  }

}

