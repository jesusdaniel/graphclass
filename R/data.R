#' fMRI brain networks of the COBRE dataset
#'
#' The COBRE dataset contains the fMRI brain networks of 124  subjects
#' and class labels of schizophrenia brain disease status.
#'
#' @format \code{COBRE.data} is a list with two elements
#' \describe{
#'   \item{X.cobre}{A matrix containing the upper triangular part of networks. Each row represents a subject,
#'   and the columns represent edges. The edge weights are calculated based on the correlation networks of the
#'   fMRI time series, after which a rank transformation among all weights in the network and a normalization
#'   between subjects is computed. The networks are composed of 263 nodes and 34453 different edges.
#'   For a description of the preprocessing steps to obtain the network edge weights, see \insertCite{relion2017network;textual}{graphclass}.
#'   }
#'   \item{Y.cobre}{Class labels of the subjects in the dataset. \code{Y=1} represents schizophrenia status.}
#' }
#' @source \url{http://fcon_1000.projects.nitrc.org/indi/retro/cobre.html}
"COBRE.data"





#' fMRI brain networks of the UMich dataset
#'
#' The UMich dataset contains the fMRI brain networks of 78  subjects
#' and class labels of schizophrenia brain disease status.
#'
#' @format \code{UMich.data} is a list with two elements
#' \describe{
#'   \item{X.cobre}{A matrix containing the upper triangular part of networks. Each row represents a subject,
#'   and the columns represent edges. The edge weights are calculated based on the correlation networks of the
#'   fMRI time series, after which a rank transformation among all weights in the network and a normalization
#'   between subjects is computed. Each network has 264 labeled nodes and 34716 different edges.
#'   For a description of the preprocessing steps to obtain the network edge weights, see \insertCite{relion2017network;textual}{graphclass}.
#'   }
#'   \item{Y.cobre}{Class labels of the subjects in the dataset.  \code{Y = 1} represents schizophrenia status.}
#' }
"UMich.data"


#' Node assignments of the Power parcellation
#'
#' Table containing the node assignments to brain systems of the 264 regions of interest from \insertCite{power2011functional;textual}{graphclass}. 
#' These nodes were used to construct COBRE.data and UMich.data. Note that node 75 is missing in the COBRE data.
#'
#' @format \code{power.parcellation} is a data frame in which the nodes represent nodes. The data frame has three columns:
#' \describe{
#'   \item{Master.Assignment}{Brain system number}
#'   \item{Color}{Color representing the system.}
#'   \item{Suggested.System}{Suggested brain system.}
#'   
#' }
#' @source \url{http://www.nil.wustl.edu/labs/petersen/Resources_files/Consensus264.xls}
"power.parcellation"

