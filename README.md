# TEMR

***Installation***  
`devtools::install_github("hhoulei/TEMR")`  

***Toy Example***  
`betaXG <- NULL`  
`sebetaXG <- NULL`  
`sebetaYG <- NULL`  
`for(i in 1:5){`  
`  betaXG <- cbind(betaXG,runif(100,0.05,0.2))`  
`  sebetaXG <- cbind(sebetaXG,runif(100,0.01,0.1))`  
`  sebetaYG <- cbind(sebetaYG,runif(100,0.01,0.1))`  
`}`  
`betaYG <- betaXG*0.2 + rnorm(100)`  
`race_name <- c('target1','target2','target3','target4','auxiliary')`  
`colnames(betaXG) <- race_name`  
`colnames(betaYG) <- race_name`  
`colnames(sebetaXG) <- race_name`  
`colnames(sebetaYG) <- race_name`  
`rho <- NULL`  
`meth <- 'CG'`  
`RESULT <- TEMR(betaXG,betaYG,sebetaXG,sebetaYG,rho,meth)`  
`RESULT`  


***Citation***:  
TEMR: Trans-ethnic Mendelian Randomization Method using Large-scale GWAS Summary Datasets
Lei Hou, Sijia Wu, Zhongshang Yuan, Fuzhong Xue, Hongkai Li

Please contact houlei@pku.edu.cn for any questions. We will continue to update this R package and reduce the problems that may be encountered during its installation.
