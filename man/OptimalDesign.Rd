\name{OptimalDesign}
\alias{OptimalDesign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimal Sampling Design
}
\description{
Visualize the effect of different combinations of the dockside intercept sample size and reporting rate on PSE under certain budget

}
\usage{
OptimalDesign(CVy, Mean_dockside, cost_ratio, RelBudget,
              p1_obs = NULL, n_obs = NULL,
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
  \item{cost_ratio}{
Relative cost for increasing one percent reporting rate in terms of increasing the number of dockside intercept trips
}
  \item{RelBudget}{
Relative budget in terms of the dockside intercept sample size
}
  \item{p1_obs}{
Current reporting rate, optional variable from the self-reported sample
}
  \item{n_obs}{
Current dockside intercept sample size, optional variable
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
Optional variable, different type of measurement error model assumed on the self-reported sample. "CME" refers to the classical measurement model, "Berkson" refers to the Berkson model. If set to NULL, the self reports are assumed to be accurate in terms of CV
}
}
\details{
The function serves two purposes: 1) Pilot study, which requires set 'n2_obs' as NULL. Under this case, the arguments about the self-repoted sample are optional. 2) Optimal sampling design for improving current estimation precision, which requires the arguments of 'p1_obs', 'n2_obs' and at least one pair of arguments:'Mean_report' and 'CVy_report' or 'EstMean_report' and 'EstCVy_report'. When either one pair of these arguments are missing, 'type' is used to calculate the missing arguments from the other pair.
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
OptimalDesign(CVy = 0.68, Mean_dockside = 9, cost_ratio = 50, RelBudget = 1500)
OptimalDesign(CVy = 0.68, Mean_dockside = 9,cost_ratio = 50, RelBudget = 800,
              Mean_report= 10,CVy_report = 0.76,p1_obs = 0.15, n_obs = 400)
OptimalDesign(CVy = 0.68, Mean_dockside = 9,cost_ratio = 50, RelBudget = 800,
              Mean_report= 10,CVy_report = 0.76,p1_obs = 0.02, n_obs = 400,
              R = 0.66,deff  =1.4)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
