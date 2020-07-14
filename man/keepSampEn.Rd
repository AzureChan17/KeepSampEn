\name{KeepNaSampEn}
\alias{KeepNaSampEn}
\title{An improved sample entropy algorithm for signals with missing values}
\usage{
KeepNaSampEn(x, m = 2, r = 0.15)
}
\arguments{
  \item{x}{A numeric vector, with data for a regularly spaced time series. Missing value must be shown as NA.}
  \item{m}{Vector of positive integers giving the window size for the entropy calculations in the second step of the algorithm:  the number of consecutive _bins_ over which similarity between subsequences is of interest. Typical values are 1, 2, or 3.}
  \item{r}{Vector of coefficients for similarity thresholds, such as r=0.15, r*sd(x) must be in the same units as 'x'. Averages in two bins are defined to be similar if they differ by 'r*sd(x)' or less. NOTE: Currently only a single threshold is allowed per run; i.e., 'r' must be a scalar.}
}
\value{
A data frame with with one row for each combination of 'm' and 'r'.
Columns are  "m", "r", and "SampEn" (the calculated sample entropy with KeepSampEn).
The data frame will also have an attribute "SD", the standard deviation
of 'x'. r = r*SD(x).
}
\description{
To perform the sample entropy for the signals which contain the missing values.
Singals without missing value also can be calculated by this algorithm.
Return the results as an R data frame.
reference: Dong X#, Chen C#, Geng Q, Cao Z, Chen X, Lin J, Yu J, Zhang Z, Shi Y, Zhang XD*. 2019 An improved method of handling missing values in the analysis of sample entropy for continuous monitoring of physiological signals. Entropy 21(3), 274 (doi: 3390/e21030274)
}
\examples{
# create a write noise data with missing value 
# or load data by yourself
set.seed(30)
y <- ts(rnorm(1000, 2, 0.15))
missing <- c()
missing <- sample(1:length(y), 100)
y[missing] <- NA
############################################################
## calculate the sample entropy of the rawdata
############################################################
keepNaSampEn(y, m=2, r=0.15)
keepNaSampEn$SampEn
}
