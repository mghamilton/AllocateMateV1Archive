#' allocate.mate.tilapia.ped
#' 
#' @description 
#' This function generates a mating list for a set of parents.  
#' The mating list can be generated i) to minimise the average inbreeding coefficient (F) of families generated or ii) according to assortative mating principles.
#' Inputs include a list of parents and a 3-column pedigree file specifying the ancestry of these candidates.
#' @param ped is a 3-column dataframe with the following columns (class in parentheses):
#' \itemize{
#'  \item{'ID' is the individual identifier of parents and their ancestors (character).} 
#'  \item{'DAM' is the identifier of the individual's dam (NA if unknown) (character).} 
#'  \item{'SIRE' is the identifier of the individual's sire (NA if unknown) (character).} 
#' }
#' @param parents data frame with the following columns (class in parentheses):
#' \itemize{
#'  \item{'ID' is the individual identifier (character).} 
#'  \item{'SEX' is the sex of the individual - 'M' or 'F', for male and female respectively (character).} 
#'  \item{'EBV' Estimated breeding value (EBV) of the individual - can not be NA if applying assortative mating as the method (numeric).} 
#' }
#' @param max_F is the maxiumum inbreeding coefficient allowed in the offspring of parents (numeric between 0 and 1)
#' @param method either 'min_F' (to minimise the average inbreeding in offspring) or 'assortative' (to apply assortative mating) (character)
#' @param n_fam_crosses is the maxiumum number of crosses allowed between parents from the same families (integer)
#' @return 'summary' is a data frame containing a summary of all possible families generated from matings between parents: 
#' \itemize{
#'  \item{'SELECTED' N represents unselected families, Y represents selected families (i.e. mate allocations to be generated) and All represents all possible families.}           
#'  \item{'COUNT_FAMS' count of families.}
#'  \item{MEAN_EBV' mean of family EBVs.}
#'  \item{'SD_EBV' standard deviation of family EBVs.}
#'  \item{'MIN_EBV' minimum of family EBVs.} 
#'  \item{'MIN_EBV' maximum of family EBVs.}  
#'  \item{MEAN_EBV' mean of family inbreeding coefficients (F).}
#'  \item{'SD_EBV' standard deviation of family inbreeding coefficients (F).}
#'  \item{'MIN_EBV' minimum of family inbreeding coefficients (F).} 
#'  \item{'MIN_EBV' maximum of family inbreeding coefficients (F).}  
#' } 
#' @return 'all_families' is a data frame containing details of all possible families able to be generated from matings between parents: 
#' \itemize{
#'  \item{'SIRE' male parent.}           
#'  \item{'DAM' female parent.} 
#'  \item{'F' inbreeding coefficient of family members (i.e. the 'family F').}       
#'  \item{'EBV' mean of parent EBVs (i.e. the 'family EBV').}
#'  \item{'SELECTED' N represents an unselected family (i.e. family is not to be generated), Y represents a selected family (i.e. mate allocated - family to be generated).}           
#' } 
#' @return 'optimal_families' is a data frame containing details of mate allocations (i.e.families to be generated from matings between parents): 
#' \itemize{
#'  \item{'SIRE' male parent.}           
#'  \item{'DAM' female parent.} 
#'  \item{'F' inbreeding coefficient of family members (i.e. the 'family F').}       
#'  \item{'EBV' mean of parent EBVs (i.e. the 'family EBV').}          
#' } 
#' @examples
#' #Retrieve example data
#' data(ped)
#' data(parents)
#' 
#' output <- allocate.mate.ped(ped, parents, max_F = 0.08, method = "min_F", n_fam_crosses = 1) 
#' output$summary
#' head(output$all_families)
#' head(output$optimal_families)


allocate.mate.ped <- function(ped, parents, max_F = 1, method = "min_F", n_fam_crosses = 1) {
  #mhamilton@cgiar.org
  #Feb 2021
  
  if("nadiv" %in% installed.packages()[, "Package"] == F) {install.packages("nadiv")}   
  library(nadiv)
  
  if(sum(colnames(parents) == "N_AS_PARENT") > 0) {
    if(sum(parents[,"N_AS_PARENT"] != 1) > 0) {
      stop("The allocate.mate.tilapia.ped function assumes N_AS_PARENT is equal to one for all parents  Rename or remove this column from \'parents\'")
    }
  }
  
  parents[,"N_AS_PARENT"] <- 1
  
  check.ped(ped)
  check.parents(parents)
  check.n_fam_crosses(n_fam_crosses)
  check.max_F(max_F)
  check.method(method)
  
  ped <- reduce.ped(ped = ped, parents = parents)
  ped[ped$DAM  == 0 & !is.na(ped$DAM), "DAM"]  <- NA
  ped[ped$SIRE == 0 & !is.na(ped$SIRE),"SIRE"] <- NA
  ped <- prepPed(ped)
  
  H <- makeA(ped)  
  H <- H[rownames(H) %in% parents$ID, colnames(H) %in% parents$ID]
  H <- as.matrix(H)
  
  families <- generate.fams(H = H, ped = ped, parents = parents, max_F = max_F) 
  
  if(method == "assortative") {
    output <- assortative(families = families, parents = parents, n_fam_crosses = n_fam_crosses)
  }
  
  if(method == "min_F") {
    output <- min.F(families = families, parents = parents, n_fam_crosses = n_fam_crosses)
  }
  
  return(output)
  
}