\name{CV_est}
\alias{CV_est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample Statistics Calculation
}
\description{
Calculate the mean, variance and coefficient of variation (CV) of the fishery sample, could be used for both the dockside intercept sample and self-reported sample.
}
\usage{
CV_est(Landings_pct, Landings_mean, Landings_var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Landings_pct}{
Percentage of fishing trips with non-zero landings
}
  \item{Landings_mean}{
Mean of catch among fishing trips with non-zero landings
}
  \item{Landings_var}{
Variance of catch among fishing trips with non-zero landings
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{Mean }{Mean of the sample}
\item{Var }{Variance of the sample}
\item{CV }{CV fo the sample}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Zhaoce Liu and Lynne Stokes

}
\note{
The arguments are only for the trips with non-zero landings. Users could omit this function if the information about the full sample or population is available from other source.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
CV_est(0.254,2.13,0.040889)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
