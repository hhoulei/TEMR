#'
#' @title TEMR
#'
#' @description Improving statistical power of causal effects in the one or
#' multiple target populations using a auxiliary population. Details refer to
#' TEMR: Trans-ethnic Mendelian Randomization Method using Large-scale GWAS
#' Summary Datasets by Lei Hou, Sijia Wu, Zhongshang Yuan, Hongkai Li, Fuzhong Xue
#'
#' @param betaXG nSNP*nRace matrix of G-X association (beta coefficient),
#' the last column is auxiliary population.
#'
#' @param betaYG nSNP*nRace matrix of G-Y association (beta coefficient),
#' the last column is auxiliary population.
#'
#' @param sebetaXG nSNP*nRace matrix of G-X association (standard error),
#' the last column is auxiliary population.
#'
#' @param sebetaYG nSNP*nRace matrix of G-Y association (standard error),
#' the last column is auxiliary population.
#'
#' @param rho nRace*nRace trans-ethnics genetic correlation matrix.
#'
#' @param meth Optimization method to use, including "Nelder-Mead", "BFGS",
#' "CG", "L-BFGS-B", "SANN" and "Brent".
#'
#' @importFrom MASS ginv
#' @importFrom stats4 mle
#'
#' @return A list object,
#'     resultall are the results of causal effect estimation in the
#'     target populations using MRTP;
#'     betaXY_IVW is the causal effect estimation of all populations using IVW;
#'     betaXY_Egger is the causal effect estimation of all populations using MR-Egger;
#'     betaXY_Final is the causal effect estimation of all populations
#'     removing the impact of pleiotropy.
#'
#' @examples
#' betaXG <- NULL
#' sebetaXG <- NULL
#' sebetaYG <- NULL
#' for(i in 1:5){
#'   betaXG <- cbind(betaXG,
#'                   runif(100,0.05,0.2)
#'   )
#'   sebetaXG <- cbind(sebetaXG,
#'                     runif(100,0.01,0.1)
#'   )
#'   sebetaYG <- cbind(sebetaYG,
#'                     runif(100,0.01,0.1)
#'   )
#' }
#' betaYG <- betaXG*0.2 + rnorm(100)
#' race_name <- c('target1','target2','target3','target4','auxiliary')
#' colnames(betaXG) <- race_name
#' colnames(betaYG) <- race_name
#' colnames(sebetaXG) <- race_name
#' colnames(sebetaYG) <- race_name
#' rho <- NULL
#' meth <- 'CG'
#' TEMR(betaXG,betaYG,sebetaXG,sebetaYG,rho,meth)
#'
#' @export
#'


TEMR <- function(betaXG,betaYG,sebetaXG,sebetaYG,rho,meth){

  race_name <- colnames(betaXG)

  race <- ncol(betaXG)
  gb <- nrow(betaXG)

  P_pleio <- NULL
  alpha <- NULL
  betaXY_Egger <- NULL
  for(i in 1:race){
    Fit <- lm(betaYG[,i]~betaXG[,i],weights = (sebetaYG[,i])^(-2))
    P_pleio <- c(P_pleio,summary(Fit)$coef[1,4])
    betaXY_Egger <- rbind(betaXY_Egger,
                          summary(Fit)$coef[2,])
    alphai <- betaYG[,i]-betaXG[,i]*summary(Fit)$coef[1,1]
    alpha <- cbind(alpha,alphai)
  }

  Wald_Ratio <- function(betaXG,betaYG,sebetaXG,sebetaYG){
    WR <- betaYG/betaXG
    varWR_2 <- (betaYG^2)*(sebetaXG^2)/(betaXG^4)+
      (sebetaYG^2)/(betaXG^2)

    return(data.frame(WR=WR,
                      varWR_2=varWR_2))
  }

  BetaWR <- NULL
  seBetaWR <- NULL
  betaXY_IVW <- NULL
  betaXY_Final <- NULL
  for(i in 1:race){
    Fit <- lm(betaYG[,i]~betaXG[,i]-1,weights = (sebetaYG[,i])^(-2))
    betaXY_IVW <- rbind(betaXY_IVW,summary(Fit)$coef[1,])

    if(P_pleio[i]<0.05){
      once <- Wald_Ratio(betaXG[,i],betaYG[,i]-alpha[,i],sebetaXG[,i],sebetaYG[,i])
      BetaWR <- cbind(BetaWR,once$WR)
      seBetaWR <- cbind(seBetaWR,once$varWR_2)
      betaXY_Final <- rbind(betaXY_Final,
                            betaXY_Egger[i,])
    }else{
      once <- Wald_Ratio(betaXG[,i],betaYG[,i],sebetaXG[,i],sebetaYG[,i])
      BetaWR <- cbind(BetaWR,once$WR)
      seBetaWR <- cbind(seBetaWR,once$varWR_2)
      betaXY_Final <- rbind(betaXY_Final,
                            summary(Fit)$coef[1,])
    }
  }

  if(is.null(rho)){
    rho <- matrix(NA,ncol=race,nrow=race)
    for(k1 in 1:race){
      for(k2 in 1:race){
        z1 <- BetaWR[,k1]/seBetaWR[,k1]
        z2 <- BetaWR[,k2]/seBetaWR[,k2]
        rho[k1,k2] <- cor(z1,z2)
      }
    }
  }

  sigma_b_all <- list()
  for(i in 1:gb){
    once <- matrix(NA, nrow=race,ncol=race)
    for(k1 in 1:race){
      for(k2 in 1:race){
        if(k1==k2) once[k1,k2] <- seBetaWR[i,k1]^2
        else once[k1,k2] <- seBetaWR[i,k1]*rho[k1,k2]*seBetaWR[i,k2]
      }
    }
    sigma_b_all[[i]] <- once
  }


  resultall <- NULL
  for(k in 1:(race-1)){

    nll <- function(bbe1) {

      bbe <- betaXY_Final[-k,1]

      sigma_bb <- list()
      pp <- NULL
      det_sigma <- NULL

      for (i in 1:gb) {

        sigma_a0 <- as.numeric(sigma_b_all[[i]][k, k])
        sigma_a1 <- sigma_b_all[[i]][-k, k]
        sigma_aa <- sigma_b_all[[i]][-k, -k]
        inv_sigma <- ginv(sigma_aa)
        sigma_bb <- t(sigma_a1) %*% inv_sigma %*% sigma_a1
        det_sigma[i] <- as.numeric(sigma_a0-sigma_bb)

        diff <- BetaWR[i,k] - bbe1 - sigma_a1 %*% inv_sigma %*% c(BetaWR[i,-k] - bbe)
        pp[i] <- as.numeric(diff * diff / det_sigma[i] )

      }

      sum(gb * log(2 * pi) + 0.5 * log(abs(det_sigma)) + 0.5 * pp)

    }

    fit <- mle(minuslog=nll,
               start=list(bbe1=0),
               method = meth)
    res <- fit@coef

    fit1 <- mle(minuslog=nll,
                start=list(bbe1=0),
                fixed = list(bbe1=0),
                method = meth)

    stat_beta1=  2 * (fit@min - fit1@min)
    pvalue_beta1 = pchisq(-stat_beta1,1,lower.tail=F)

    resultall[[k]] <- list(beta=res,
                           stat=stat_beta1,
                           pvalue=pvalue_beta1)

  }

  rownames(betaXY_IVW) <- race_name
  rownames(betaXY_Egger) <- race_name
  rownames(betaXY_Final) <- race_name
  names(resultall) <- race_name[1:(race-1)]

  return(list(resultall=resultall,
              betaXY_IVW=betaXY_IVW,
              betaXY_Egger=betaXY_Egger,
              betaXY_Final=betaXY_Final
  ))

}
