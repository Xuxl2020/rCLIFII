% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models_internal.R
\name{models_internal}
\alias{models_internal}
\alias{Model1}
\alias{Model2}
\alias{Model3}
\alias{Model4}
\alias{Model5}
\alias{Model6}
\alias{lir.model.res}
\alias{lar.model.res}
\title{Build the models}
\usage{
Model1(mij, nij, tauij, mtau = 1000)

Model2(mij, nij, tauij, mtau = 1000)

Model3(mij, nij, tauij, mtau = 1000)

Model4(Aij, Ai, tauij, mtau = 1000)

Model5(Aij, Ai, tauij, mtau = 1000)

Model6(Aij, Ai, tauij, mtau = 1000)

lir.model.res(
  model,
  mij,
  nij,
  tauij,
  mtau,
  model_cl_fun = NULL,
  cl.H = NULL,
  model.K = NULL
)

lar.model.res(
  model,
  Aij,
  Ai,
  tauij,
  mtau,
  model_cl_fun = NULL,
  cl.H = NULL,
  model.K = NULL
)
}
\arguments{
\item{mij}{The number of individuals identified at both time t_i and t_j}

\item{nij}{The number of individuals identified at time t_i}

\item{tauij}{Time lag between time t_i and t_j}

\item{mtau}{The maximum allowable lag time.}

\item{Aij}{The number of observed associated individuals (in pairs) observed at both time t_i and t_j}

\item{Ai}{The number of observed associated individuals (in pairs) at time t_i}

\item{model}{Models of LIR and LAR, or model is formulated by yourself}

\item{model_cl_fun}{If you formulate your model, please input function to calculate the composite likelihood about your model}

\item{cl.H}{If you formulate your model, please input the sensitivity matrix with respect to parameters in your model}

\item{model.K}{If you formulate your model, please input the number of parameters in your model}
}
\description{
Build the models
}
