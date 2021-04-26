#' fMRI brain networks of the COBRE dataset
#'
#' This dataset contains functional connectivity brain networks of 124 subjects
#' and corresponding class labels indicating schizophrenia status (70 healthy controls and 54
#' schizophrenic subjects) obtained from the COBRE dataset.
#' 
#' @docType data
#' 
#' @format \code{COBRE.data} is a list with three elements
#' \describe{
#'   \item{X.cobre}{A matrix containing the upper triangular part of the networks in column major order and without the
#'   diagonal entries. The networks are composed of 263 nodes obtained from the Power parcellation \insertCite{power2011functional}{graphclass}
#'   and 34453 different edges. Each row corresponds to a subject in the data,
#'   and the columns correspond to edges in the networks. The edge weights are the Fisher-transformed correlation 
#'   between the fMRI time series of the nodes. Nuisance covariates like age, gender, motion (meanFD and meanFDquad) and handedness
#'   have been regressed out. For a description of the preprocessing steps to obtain the network edge weights, 
#'   see \insertCite{relion2017network;textual}{graphclass}.
#'   }.
#'   \item{Y.cobre}{Class labels of the subjects in the dataset. \code{Y=1} represents schizophrenia status.}
#'   \item{subject.label}{Subject ID number in the COBRE dataset.}
#' }
#' 
#' @usage
#' data(COBRE.data)
#' 
#' @examples 
#' data(COBRE.data)
#' plot_adjmatrix(COBRE.data$X.cobre[1,])
#' 
#' @keywords datasets
#' 
#' @source \url{http://fcon_1000.projects.nitrc.org/indi/retro/cobre.html}
#' @references 
#' \insertRef{Aine2017}{graphclass}
#' 
#' \insertRef{relion2017network}{graphclass}
#' 
#' \insertRef{power2011functional}{graphclass}
"COBRE.data"





#' fMRI brain networks of the UMich dataset
#'
#' This dataset contains functional connectivity brain networks of 79 subjects
#' and corresponding class labels indicating schizophrenia status (40 healthy controls and 39
#' schizophrenic subjects). The fMRI data used to obtain these networks were collected by Professor
#' Stephen F. Taylor's lab at the University of Michigan.
#'
#' @docType data
#' 
#' @format \code{UMich.data} is a list with two elements
#' \describe{
#'   \item{X.cobre}{A matrix containing the upper triangular part of the networks in column major order and without the
#'   diagonal entries. The networks are composed of 264 nodes obtained from the Power parcellation \insertCite{power2011functional}{graphclass} 
#'   and 34716 different edges. Each row corresponds to a subject in the data,
#'   and the columns correspond to edges in the networks. The edge weights are the Fisher-transformed correlations 
#'   between the fMRI time series of the nodes. Nuisance covariates like age, gender, motion (meanFD and meanFDquad) and handedness
#'   have been regressed out. For a description of the preprocessing steps to obtain the network edge weights, 
#'   see \insertCite{relion2017network;textual}{graphclass}.
#'   }
#'   \item{Y.cobre}{Class labels of the subjects in the dataset.  \code{Y = 1} represents schizophrenia status.}
#' }
#' 
#' @keywords datasets
#' 
#' @usage 
#' data(UMich.data)
#' 
#' @examples 
#' data(UMich.data)
#' plot_adjmatrix(UMich.data$X.umich[1,])
#' 
#' @references 
#' \insertRef{relion2017network}{graphclass}
#' 
#' \insertRef{power2011functional}{graphclass}
"UMich.data"


#' Node assignments of the Power parcellation
#'
#' Table containing the node assignments to brain systems of the 264 regions of interest from \insertCite{power2011functional;textual}{graphclass}. 
#' These nodes were used to construct COBRE.data and UMich.data. Note that node 75 is missing in the COBRE data.
#'
#' @docType data
#' 
#' @format \code{power.parcellation} is a data frame in which the nodes represent nodes. The data frame has three columns:
#' \describe{
#'   \item{Master.Assignment}{Brain system number}
#'   \item{Color}{Color representing the system.}
#'   \item{Suggested.System}{Suggested brain system.}
#'   
#' }
#' 
#' @examples 
#' data(power.parcellation)
#' 
#' @keywords datasets
#' 
#' @source \url{http://www.nil.wustl.edu/labs/petersen/Resources_files/Consensus264.xls}
#' @references 
#' \insertRef{power2011functional}{graphclass}
#' 
"power.parcellation"

