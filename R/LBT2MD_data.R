# 135 characters #####################################################################################################################
#' @title MrDAG dataset: Lifestyle and behavioural exposures that might impact mental health phenotypes
#' @description The dataset contains lifestyle and behavioural traits that are considered exposures for the risk of mental health 
#' phenotypes. As outcomes, seven mental health phenotypes are considered, including (in alphabetical order): attention deficit 
#' hyperactivity disorder (ADHD), anorexia nervosa (AN), autism spectrum disorder (ASD), bipolar disorder (BD), cognition (COG), 
#' major depressive disorder (MDD), and schizophrenia (SCZ). As exposures, six lifestyle and behavioural traits that have 
#' previously been investigated for their protective/risk effects on mental health are considered, including (in alphabetical order): 
#' alcohol consumption (ALC), education (in years) (EDU), leisure screen time (LST), physical activity (PA), lifetime smoking index 
#' (SM), and sleep duration (SP)
#'
#' @docType data
#'
#' @format A list \code{list(data, gencov)}:
#' \describe{
#'   \item{\code{data}}{ A matrix of 708 independent Instrumental Variables (IVs) (rows) by trait (columns, outcomes first then
#'         exposures), selected to be associated at genome-wide significance with the exposures after pruning or clumping. For
#'         details, see \insertCite{Zuber2025;textual}{MrDAGSuite} }
#'   \item{\code{gencov}}{ A list of 708 13-by-13 per-IV sample-overlap covariance matrices: the off-diagonal entries are the LDSC cross-trait intercepts times the standardised standard errors, and the diagonal is each trait's squared standardised standard error (stored, but not used by \code{MrDAG}, which corrects only the off-diagonal sample-overlap terms), one per IV in
#'         \code{data}, in the same trait order, ready to pass to \code{\link{MrDAG}} }
#' }
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @examples
#' \dontrun{
#' # Example
#' data(LBT2MD_data)
#' head(LBT2MD_data$data)
#' }

"LBT2MD_data"