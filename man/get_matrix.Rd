% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extra_functions.R
\encoding{UTF-8}
\name{get_matrix}
\alias{get_matrix}
\title{Converts a vector into an adjacency matrix}
\usage{
get_matrix(beta, type = c("undirected", "directed"))
}
\arguments{
\item{beta}{Vectorized adjacency matrix. If the network is undirected,
the vector is assumed to represent the upper triangular part of the adjacency 
matrix in column major order. For undirected network, the vector contains the
entries ordered by column and excluding the diagonal elements.}

\item{type}{Specifies whether the vector represents an \code{undirected}
or \code{directed} network. Default is \code{undirected}.}
}
\value{
Adjacency matrix.
}
\description{
Given a vector that encodes an adjacency matrix, returns the matrix representation.
}
\examples{

# Obtain the adjacency matrix of a COBRE data subject
data(COBRE.data)
A <- get_matrix(COBRE.data$X.cobre[1,], type = "undirected")
plot_adjmatrix(A)

}
