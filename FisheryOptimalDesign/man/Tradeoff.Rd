\name{Tradeoff}
\alias{Tradeoff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Trade-off between the Dockside Intercept Sample Size and Reporting Rate
}
\description{
Visualize all possible combinations of the dockside intercept sample size and reporting rate for certain PSE
}
\usage{
Tradeoff(CVy, MeanL_dockside, target_PSE, target_n2 = NULL, p1_obs = NULL,
              Mean_report = NULL, CVy_report = NULL,
              EstMean_report = NULL, EstCVy_report = NULL,
              deff = NULL, R = NULL, type = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CVy}{
CV of the dockside intercept sample
}
  \item{Mean_dockside}{
Mean fish catch of the dockside intercept sample
}
  \item{target_PSE}{
Desired estimation precision
}
  \item{target_n2}{
Desired dockside intercept sample size, optional variable
}
  \item{p1_obs}{
Current reporting rate, optional variable from the self-reported sample
}
  \item{Mean_report}{
Mean fish catch of the self-reported sample, optional variable if the information is available for the self-reported sample
}
  \item{CVy_report}{
CV of the self-reported sample, optional variable if the information is available for the self-reported sample
}
  \item{EstMean_report}{
Estimated actual mean fish catch of the self-reported sample, optional variable if the information is available for the self-reported sample
}
  \item{EstCVy_report}{
Estimated CV of the self-reported sample, optional variable if the information is available for the self-reported sample
}
\item{deff}{
Design effect of the dockside intercept sample. If set to NULL, will be assumed as 2.5.
}
  \item{R}{
Correlation coefficient of the actual fish catch and its self-reported value among the self-reported sample. If set to NULL, will be assumed as 1
}
  \item{type}{
Optional variable of different type of measurement error model assumed on the self-reported.sample "CME" refers to the classical measurement model, "Berkson" refers to the Berkson model. If set to NULL, the self reports are assumed to be accurate in terms of CV
}
}
\details{
The function could be used either with or without information about the self-reported sample. However, when such information is available, the resulting scenario could avoid extra assumptions. To describe the self-reported sample, 'p1_obs' is required and at least one pair of arguments:'Mean_report' and 'CVy_report' or 'EstMean_report' and 'EstCVy_report' is required. When either one pair of arguments are missing, 'type' is used to calculate the missing arguments from the other pair.
}
\value{
None
}
\references{
Liu, B., Stokes, L., Topping, T., & Stunz, G. (2017). Estimation of a total from a population of unknown size and application to estimating recreational red snapper catch in Texas. Journal of Survey Statistics and Methodology, 5(3), 350-371.
}
\author{
Zhaoce Liu and Lynne Stokes
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
Tradeoff(CVy = 0.68, Mean_dockside = 9, target_PSE = 0.15)
Tradeoff(CVy = 0.68, Mean_dockside = 9,target_PSE = 0.15,target_n2 = 421)
Tradeoff(CVy = 0.68, Mean_dockside = 9, target_PSE = 0.15,
         target_n2 = 421,Mean_report= 10,CVy_report = 0.76,
         p1_obs = 0.4, R = 0.66,deff  =1.4,type = "CME")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
