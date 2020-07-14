# keepSampEn
 A R function for sample entropy calculation when signal contains missing values
## Background
  In common R packages, when a missing value appears in a sequence, the sample entropy cannot be calculated directly. We improved the original algorithm to calculate the signal with missing value. Compared with other traditional missing value processing, our method: KeepSampEn can minimize the error from missing value. It has been proved that KeepSampEn can more accurately calculate sample entropy regardless of physiological signal type, data size and generating mechanism when the signal with missing values. 
 
 ## Install package
   Before installing this package, users need to insrall the Rtools by themselves. If users have installed the Rcpp package in R language, they need to remove it first. If “remove.packages(Rcpp)” cannot completely remove the Rcpp package, users can use "Sys.geyenv(R_LIBS_USER)" to find the address of Rcpp package and delete it. 
   
        install(devtools)
        library(devtools)
        devtools::install_github("AzureChan17/keepSampEn")
        library(keepSampEn)
 ## Data
 This algorithm can calculate the sample entropy as original algorithm when the signal doesn't contain the missing value. When the signal contains missing value, missing values need to be represented by NA. Here we conduct a write noise signal with missing value to be an example.
 
        set.seed(30)
        y <- ts(rnorm(1000, 2,0.5))  #conduct a write noise signal with 1000 points.
        missing <- c()
        missing <- sample(1:length(y), 100)  #generate 100 positions for missing values
        y[missing] <- NA  #vector y is a write noise which contain 100 missing values
 
 ## Usage 
        
        KeepNaSampEn(y,2,0.15)
        #   m       r SampEn
        # m 2 0.15*SD  2.463
        KeepNaSampEn(y,2,0.15)$SampEn
        # [1] 2.463   ## The calculation result is not unique
        
 ## Reference
   The algorithm can be found from the reference: Dong, X.; Chen, C.; Geng, Q.; Cao, Z.; Chen, X.; Lin, J.; Jin, Y.; Zhang, Z.; Shi, Y.; Zhang, X.D. An Improved Method of Handling Missing Values in the Analysis of Sample Entropy for Continuous Monitoring of Physiological Signals. Entropy 2019, 21, 274. https://www.mdpi.com/1099-4300/21/3/274 This paper also provide an C code and guideline for users.
   
 
