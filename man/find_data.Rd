% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMR-find_data.R
\name{find_data}
\alias{find_data}
\alias{find_data.default}
\alias{find_data.lm}
\alias{find_data.glm}
\alias{find_data.coxph}
\alias{find_data.aalen}
\alias{find_data.timecox}
\alias{find_data.cox.aalen}
\alias{find_data.aareg}
\title{Internal Method for Extracting Data for IMR}
\usage{
find_data(m)

\method{find_data}{default}(m)

\method{find_data}{lm}(m)

\method{find_data}{glm}(m)

\method{find_data}{coxph}(m)

\method{find_data}{aalen}(m)

\method{find_data}{timecox}(m)

\method{find_data}{cox.aalen}(m)

\method{find_data}{aareg}(m)
}
\arguments{
\item{m}{a model object}
}
\value{
a list with the following three entries
  \describe{
    \item{X}{the fitted values of the model. In hazard models this is the
      time points corresponding to observed events.}
    \item{Y}{a list with the raw and standerdized residuals. In hazard models
      this is the standerdized residuals found when fitting a line through the
      estimated cumulative hazard.}
    \item{name}{a ("short") string describing the model. This value will be
      used as the model name in the final output.}
  }
}
\description{
An internal generic method for extracting fitted response and standerdize
residuals.
}
\details{
This function finds the relevant model data and returns it in a suitable
manner to the \code{\link{IMR}} function.
}
\keyword{internal}
