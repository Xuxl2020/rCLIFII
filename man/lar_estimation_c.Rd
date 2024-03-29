% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_c.R
\name{lar_estimation_c}
\alias{lar_estimation_c}
\title{Estimation of variance inflation factor in QAIC under LAR models}
\usage{
lar_estimation_c(g_m, g_n, Aij, Ai, tauij, mtau = 1000)
}
\arguments{
\item{g_m}{\eqn{\sum_{i,j} \{A_{t_i,t_j}|\tau_{ij}=\tau\}}, with \eqn{\tau_{ij}=|t_i-t_j|}}

\item{g_n}{\eqn{\sum_{i,j} \{A_i|\tau_{ij}=\tau\}}, with \eqn{\tau_{ij}=|t_i-t_j|}}

\item{Aij}{Number of observed associated individuals (in pairs) observed at both time t_i and t_j}

\item{Ai}{Number of observed associated individuals (in pairs) at time t_i}

\item{tauij}{Time lag between time t_i and t_j}

\item{mtau}{Maximum allowable lag time}
}
\description{
Estimation of variance inflation factor in QAIC under LAR models
}
