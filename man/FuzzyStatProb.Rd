% File FuzzyStatProb/man/FuzzyStatProb.Rd

\name{FuzzyStatProb}
\alias{fuzzyStationaryProb}
\alias{print.FuzzyStatObj}
\alias{summary.FuzzyStatObj}
\title{Fuzzy stationary probabilities of Markov chains from observations}
\description{
Computation of LR fuzzy numbers representing fuzzy stationary probabilities of an unknown Markov chain from which a sequence of observations has been drawn.
The fuzzy Markov chain considered during the processing follows the approach proposed by J. Buckley (see the reference section).
}
\author{Pablo J. Villacorta Iglesias, Department of Computer Science and Artificial Intelligence, University of Granada (Spain).

\email{pjvi@decsai.ugr.es} - \url{http://decsai.ugr.es/~pjvi}
}
\usage{
fuzzyStationaryProb(data,options,step=0.05)
\method{print}{FuzzyStatObj}(x, ...)
\method{summary}{FuzzyStatObj}(object, ...)
}
\arguments{
\item{data}{ 
An array of natural numbers \eqn{\ge 1} representing the observed states of the chain at consecutive time points. The function first coerces the elements to integer. 
The minimum value should equal 1; otherwise all elements are substracted the difference between 1 and the minimum value to adapt the range and force 1 to be the minimum. 
After such transformation, the maximum of the vector is taken as the number of states of the chain \emph{n}.
}
\item{options}{ A tagged list containing the following parameters:
\itemize{
\item \code{verbose}: boolean, set to TRUE if progress information should be printed during the process. It is set to FALSE if this option is not specified. 
\item \code{states}: an array of naturals \eqn{\ge 1} indicating the states for which the stationary distribution should be computed. The values should match those specified in the
\code{data} argument, and will be coerced to integer, checked and transformed if necessary in the same way as with \code{data}. 
If this option is not specified, the fuzzy stationary probabilities are computed for every state of the chain.
\item \code{regression}: a string with the type of the regression to be applied at the end of the algorithm for fitting the membership functions of the fuzzy stationary probabilities. 
Possible values are 'linear', 'quadratic', 'cubic', 'gaussian', 'spline' and 'piecewise' (piecewise linear interpolation). 
In all cases (including the gaussian), a different curve is fitted for each side of tue fuzzy number.
The \code{gaussian} option fits curves of the form \eqn{\mu(x)} = exp \eqn{( -1/2 |(x-c)/s|^m)}.
The \code{spline} option performs interpolation by a monotone cubic spline according to the Hyman method (see \code{splinefun} documentation) while \code{piecewise} computes a
piecewise linear membership function by connecting consecutive points of the \eqn{\alpha}-cuts with straight lines, using the built-in \code{PiecewiseLinearFuzzyNumber} subclass of the 
\pkg{FuzzyNumbers} package. If this option is not specified, quadratic regression is carried out by default.
\item \code{acutsonly}: boolean, set to TRUE if no regression should be done after computing the \eqn{\alpha}-cuts. This option is set to FALSE if not specified.
\item \code{parallel}: boolean, set to TRUE in order to make use of all available physical cores and reduce computation time using the built-in \code{parLapply} function of the \pkg{parallel} package.
Each processor takes care of all the computations involving one of the values of \eqn{\alpha} that have to be sampled. This option is set to FALSE if not specified.
\item \code{cores}: integer representing the maximum number of cores that can be used when running in parallel. If \code{parallel} is set to \code{FALSE}, 
this option is ignored if specified. If \code{parallel} is set to \code{TRUE}, then this value must be strictly positive. If this option is not specified when 
\code{parallel} has been set to \code{TRUE}, then all physical cores are used.
}
}
\item{step}{ 
Step size for sampling \eqn{\alpha} when computing the \eqn{\alpha}-cuts. The smallest \eqn{alpha} is always present and equals 0.001, 
and the rest of values are calculated as \eqn{\alpha = k } \code{step} for \eqn{k \ge 1}. 
The greatest sampled value that is always present as well is \eqn{\alpha = 0.999}.
It is set to 0.05 when this option is not specified.
}
\item{x}{
An object returned by the fuzzyStationaryProb function
}
\item{object}{
An object returned by the fuzzyStationaryProb function
}
\item{...}{
Further arguments passed to the summary and prints functions.
}
}
\details{
Given a sequence of consecutive observations of the state of the chain, a fuzzy transition matrix is constructed according to the approach proposed in
J. Buckley's \emph{Fuzzy Probabilities} book. Fuzzy transition probabilities are constructed as the superposition of intervals (\eqn{\alpha}-cuts), 
which in this case represent simultaneous confidence intervals for multinomial proportions, and are computed using the input sequence of observations 
drawn from the chain. For each value of \eqn{\alpha}, the \eqn{\alpha}-cuts of such fuzzy transition probabilities define a matrix space
where we seek for the the matrices producing respectively the minimum and maximum possible stationary probability for each state of the chain,
using heuristic optimization tools (Differential Evolution). 
Both points define a closed real interval that is indeed an \eqn{\alpha} cut of the output fuzzy number representing the fuzzy stationary probability for that state.
Solving these problems for different \eqn{\alpha} allows to reconstruct the fuzzy stationary probabilities from their \eqn{\alpha}-cuts, 
applying the decomposition theorem. Regression is applied at the final stage to compute the membership functions of the stationary probabilities.
}
\value{
The value returned by \code{fuzzyStationaryProb} is an object of the new S3 class \code{FuzzyStatObj}, which is a tagged list with the following components:

\item{fuzzyStatProb}{A list of \code{FuzzyNumber} objects (see the \pkg{FuzzyNumbers} package cited before). 
The length of the list equals that of the \code{states} tag of the \code{options} argument. The object at a given position \code{i} corresponds to the 
fuzzy stationary probability of the state indicated at position i of the \code{states} vector.
If any of the states indicated in the \code{states} option is not found in the \code{data} input vector, the corresponding position in \code{fuzzyStatProb} 
will be NA. If the function was called with \code{acutsonly} set to TRUE, then the returned object will not have a \code{fuzzyStatProb} tag.
}
\item{acuts}{
A list of data frame objects containing the \eqn{\alpha}-cuts of every output stationary probability, represented as bidimensional points
(lowerBound,\eqn{\alpha}) and (upperBound,\eqn{\alpha}) where \eqn{\pi(\alpha) = [lowerBound, upperBound]} is an \eqn{\alpha}-cut of the fuzzy stationary number
\eqn{\pi}. The length of the list also equals that of the \code{states} tag of the \code{options} argument. Again, object at position \code{i} corresponds
to \eqn{\alpha}-cuts of the state indicated at position \code{i} of the \code{states} vector of the option list.
If any of the states indicated in the \code{states} option is not found in the \code{data} input vector, the corresponding position in \code{acuts} will be NA.
}
Class \code{FuzzyStatObj} has only two specific methods, 
\code{print} and \code{summary}, and both print exactly the same: a brief summary of the processing that has been done, including the number of states, 
number of observations, regression type, step size and time elapsed, and information
about the names of the two tags described above that should be queried to retrieve the results. 
}
\references{
Buckley, J.J. Fuzzy Probabilities: New Approach and Applications, 2nd edition, volume
115 of Studies in Fuzziness and Soft Computing. Springer, 2005.

Glaz, J. and C.P. Sison. Simultaneous confidence intervals for multinomial proportions. 
Journal of Statistical Planning and Inference 82:251-262 (1999). 

May, W.L. and W.D. Johnson. Constructing two-sided simultaneous confidence intervals for 
multinomial proportions for small counts in a large number of cells. Journal of Statistical Software 5(6) (2000).
Paper and code available at \url{http://www.jstatsoft.org/v05/i06}.

Gagolewski M. FuzzyNumbers Package: Tools to deal with fuzzy numbers in R (2012).
Tutorial available at \url{http://www.ibspan.waw.pl/~gagolews/FuzzyNumbers/doc/FuzzyNumbers-Tutorial.pdf}
}

\examples{
# Simulate 200 observations of a 10-state Markov chain, 
# and compute fuzzy stationary probability of state 1
simula <- function(trans,N,nstates) {
 states = 1:nstates;
 sim <- rep(NA,N)
 sim[1] <- sample(states,1)
 for (i in 2:N) {
  sim[i] <- sample(states,1,prob=trans[sim[i-1],])
 }
 sim
}
nstates = 10;
trans <- matrix(data=
c(0,0.05,0,0,0.95,0,0,0,0,0,0.41,
0,0.59,0,0,0,0,0,0,0,
0,0.23,0,0.44,0,0.33,0,0,0,0,
0,0,1,0,0,0,0,0,0,0,
0.65,0,0,0,0,0,0.35,0,0,0,
0,0,0.55,0,0,0,0,0,0.45,0,
0,0,0,0,0.53,0,0,0.47,0,0,
0,0,0,0,0,0,0.48,0,0.52,0,
0,0,0,0,0,0.35,0,0.35,0,0.3,
0,0,0,0,0,0,0,0,1,0),nrow=nstates,ncol=nstates,byrow=TRUE)
set.seed(666);
data = simula(trans,200,nstates);
# Simplest case: compute only alpha-cuts for alpha=0.001 and alpha=0.999
linear = fuzzyStationaryProb(data,list(verbose=TRUE, parallel=FALSE,
states=1, regression="piecewise"),step=1)
summary(linear)
linear$fuzzyStatProb[[1]]
plot(linear$fuzzyStatProb[[1]])
points(linear$acuts[[1]])
# A more accurate approximation, using membership steps of 0.1 (takes much longer!)
\dontrun{
quadratic = fuzzyStationaryProb(data,list(verbose=TRUE, parallel=TRUE, cores = 2,
states=1, regression="quadratic"),step=0.1)
summary(quadratic)
quadratic$fuzzyStatProb[[1]]
plot(quadratic$fuzzyStatProb[[1]])
points(quadratic$acuts[[1]])
}
}
\keyword{fuzzy}
\keyword{Markov chain}
\keyword{stationary probability}
\keyword{probabilities}