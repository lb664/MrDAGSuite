# 135 characters #####################################################################################################################
#' @title simMrDAG: Data simulator of causal graphical models for joint Mendelian randomization analysis of multiple exposures and 
#' outcomes
#' @description Data simulator to generate multiple outcomes and multiple exposures with an underlying DAG structure with and without 
#' reverse causation in one- or two-sample, and individual- or summary-level designs
#'
#' @param nRep Number of simulated replicates. Default value is \code{1}
#' @param q Number of simulated outcomes. When \code{graph_Y = "dir"}, at least three outcomes are needed to simulate a directed graph with more than one edge
#' @param p Number of simulated exposures
#' @param nSubject_X Number of individuals (default \code{50,000}) in the exposures sample
#' @param ratioSubject_XY Ratio (default \code{1}) between the number of individuals in the exposures and outcomes samples. This 
#' corresponds to a two-sample design. If \code{ratioSubject_XY == NULL}, a one-sample design is considered with \code{nSubject_X} 
#' individuals
#' @param nIV Number of genetic variants (default \code{100}) used as instrumental variables
#' @param MAF Minor allele frequency used to simulate the genotypes in each sample's individual-level data. Default value is set at 
#' \code{0.05} for all genetic variants. If two values are specified, a different MAF is randomly drawn for each genetic variant 
#' centred in the first specified value of MAF, while the second value corresponds to the dispersion
#' @param LDmax LD structure between the genetic variants. By default \code{NULL}, genetic variants are simulated independently. If a 
#' positive value is specified, it corresponds to the max pairwise LD between two consecutive genetic variants, while other variants 
#' have decreasing levels of LD depending on how far apart they are
#' @param thetaRange Sought range of simulated direct effects. Default value is set at \code{(-2,2)}
#' @param thetaSparsness Proportion of simulated direct effects between exposures and outcomes. Default value is set at \code{0.2}
#' @param graph_Y Type of graph simulated within the outcomes. Directed acyclic graphical (DAG) model \code{"dir"} is set as default,
#' but an undirected graphical model \code{"undir"} with dependencies between consecutive exposures can also be designed. If
#' \code{graph_X = "reverse"}, this parameter is not used
#' @param pDAG_Y If \code{graph_Y = "dir"}, probability of connecting a node to another node with higher topological ordering within
#' the outcomes. Default is \code{0.1}. If \code{graph_X = "reverse"}, this parameter is not used
#' @param methodDAG_Y String, specifying the method used for generating the random graph within the outcomes. For details, see
#' function \code{randDAG} in \insertCite{Kalisch2012;textual}{MrDAGSuite}. If \code{graph_X = "reverse"}, this parameter is not used
#' @param dDAG_Y Positive number, corresponding to the expected number of neighbours per node within the outcomes, more precisely
#' the expected sum of the in- and out-degree. For details, see function \code{randDAG} in
#' \insertCite{Kalisch2012;textual}{MrDAGSuite}. If \code{graph_X = "reverse"}, this parameter is not used
#' @param par1DAG_Y Optional additional arguments, dependent on \code{methodDAG_Y}. If \code{NULL}, the default parametrisation is
#' used, for details see function \code{randDAG} in \insertCite{Kalisch2012;textual}{MrDAGSuite}. If \code{graph_X = "reverse"},
#' this parameter is not used
#' @param par2DAG_Y Optional additional arguments, dependent on \code{methodDAG_Y}. If \code{NULL}, the default parametrisation is
#' used, for details see function \code{randDAG} in \insertCite{Kalisch2012;textual}{MrDAGSuite}. If \code{graph_X = "reverse"},
#' this parameter is not used
#' @param graph_X Type of graph simulated within the exposures. Directed acyclic graphical (DAG) model \code{"dir"} is set as default,
#' but an undirected graphical model \code{"undir"} with dependencies between consecutive exposures can also be designed. If
#' \code{"reverse"} is specified, a DAG model between all exposures and outcomes is simulated, allowing for reverse causation. The
#' parametrization for the simulation is encoded only by the \code{graph_X} parametrization (and not by \code{graph_Y})
#' @param pDAG_X If \code{graph_X = "dir"} or \code{graph_X = "reverse"}, probability of connecting a node to another node with
#' higher topological ordering within the exposures. Default is \code{0.1}
#' @param methodDAG_X String, specifying the method used for generating the random graph within the exposures. For details, see
#' function \code{randDAG} in \insertCite{Kalisch2012;textual}{MrDAGSuite}
#' @param dDAG_X Positive number, corresponding to the expected number of neighbours per node within the exposures, more precisely
#' the expected sum of the in- and out-degree. For details, see function \code{randDAG} in \insertCite{Kalisch2012;textual}{MrDAGSuite}
#' @param par1DAG_X Optional additional arguments, dependent on \code{methodDAG_X}. If \code{NULL}, the default parametrisation is
#' used, for details see function \code{randDAG} in \insertCite{Kalisch2012;textual}{MrDAGSuite}
#' @param par2DAG_X Optional additional arguments, dependent on \code{methodDAG_X}. If \code{NULL}, the default parametrisation is
#' used, for details see function \code{randDAG} in \insertCite{Kalisch2012;textual}{MrDAGSuite}
#' @param herit_Y Heritability of the outcomes. Default value is set at \code{0.25}
#' @param betaHatRange_X Sought range of simulated genetic effects on the exposures. Default value is set at \code{(-2,2)}
#' @param herit_X Heritability of the exposures. Default value is set at \code{0.1}
#' @param gamma_Y Average mediation between dependent outcomes or correlation between consecutive ones if \code{graph_Y = "undir"}.
#' Default value is set at \code{1}. For details, see \insertCite{Zuber2025;textual}{MrDAGSuite}
#' @param rho_X Average mediation between dependent exposures or correlation between consecutive ones if \code{graph_X = "undir"}.
#' Default value is set at \code{0.8}. For details, see \insertCite{Zuber2025;textual}{MrDAGSuite}
#' @param psi_Y Impact of the confounder on the outcomes. Default value is set at \code{1}
#' @param psi_X Impact of the confounder on the exposures. Default value is set at \code{2}
#' @param saveAll If logical (default) \code{FALSE}, only summary-level data are saved, otherwise also simulated individual-level data 
#' are saved
#' @param fileName Name of the file for the \code{simMrDAG} output object (default \code{"Sim_MrDAG_object"}). If \code{NULL}, the 
#' output is not saved
#' @param filePath Directory where the \code{simMrDAG} object is saved. If the directory does not exist, \code{simMrDAG} will create 
#' it by using the current working directory as the root directory. If the directory is not specified (default \code{NULL}), the 
#' \code{simMrDAG} object is saved in the current working directory
#' @param progress If logical (default) \code{TRUE} and \code{nRep > 1}, the simulations' progress is displayed on the screen 
#' @param seed Seed to be used in the initialization of the simulator (default \code{3112021})
#'
#' @details
#' For details regarding the simulation structure, see \insertCite{Zuber2025;textual}{MrDAGSuite}
#'
#' @export
#'
#' @return For each replicate, depending on whether \code{saveAll = TRUE}, the value returned is a list object
#' \code{list(Y, X, G, betaHat_Y, seBetaHat_Y, pvalueBetaHat_Y, betaHat_X, seBetaHat_X, pvalueBetaHat_X,}
#' \code{Theta, Gamma_Y, Rho_X, DAG_Y, }
#' \code{DAG_X, A_Y, A_X, parSim, optSim)} or if \code{saveAll = FALSE},
#' \code{list(betaHat_Y, seBetaHat_Y, pvalueBetaHat_Y, betaHat_X, seBetaHat_X, pvalueBetaHat_X, Theta, Gamma_Y, Rho_X, DAG_Y, DAG_X, }
#' \code{A_Y, A_X, parSim, optSim)}:
#' \describe{
#'   \item{\code{Y}}{ Simulated individual-level outcomes }
#'   \item{\code{X}}{ Simulated individual-level exposures }
#'   \item{\code{G}}{ Simulated genotypes split into G_Y for the outcomes and G_X for the exposures, unless a one-sample design is 
#'                    selected }
#'   \item{\code{simMAF}}{ Simulated MAF for each IV split into simMAF_Y for the outcomes and simMAF_X for the exposures, unless a 
#'                         one-sample design is selected }
#'   \item{\code{betaHat_Y}}{ Simulated summary-level outcomes }
#'   \item{\code{seBetaHat_Y}}{ Standard error of simulated summary-level outcomes }
#'   \item{\code{pvalueBetaHat_Y}}{ P-value of simulated summary-level outcomes }
#'   \item{\code{betaHat_X}}{ Simulated summary-level exposures }
#'   \item{\code{seBetaHat_X}}{ Standard error of simulated summary-level exposures }
#'   \item{\code{pvalueBetaHat_X}}{ P-value of simulated summary-level exposures }
#'   \item{\code{Theta}}{ Matrix (exposures times outcomes) of simulated direct effects }
#'   \item{\code{Gamma_Y}}{ Matrix (outcomes times replicates) of simulated mediation effects between dependent outcomes }
#'   \item{\code{Rho_X}}{ Matrix (exposures times exposures) of simulated mediation effects between dependent exposures }
#'   \item{\code{DAG_Y}}{ graphNEL object that describes the dependency relations within the outcomes }
#'   \item{\code{DAG_X}}{ graphNEL object that describes the dependency relations within the exposures }
#'   \item{\code{A_Y}}{ Adjacency matrix (outcomes times outcomes) of dependency relations within the outcomes }
#'   \item{\code{A_X}}{ Adjacency matrix (exposures times exposures) of dependency relations within the exposures }
#'   \item{\code{parSim}}{ List of parameters used in the simulation \code{list(nRep, q, p, nSubject_Y, nSubject_X, ratioSubject_XY, }
#'                         \code{nIV, MAF, LDmax, thetaRange, thetaSparsness, graph_Y, pDAG_Y, dDAG_Y, methodDAG_Y, par1DAG_Y, }
#'                         \code{par2DAG_Y, graph_X, pDAG_X, dDAG_X, methodDAG_X, par1DAG_X, par2DAG_X, herit_Y, herit_X, }
#'                         \code{betaHatRange_X, gamma_Y, rho_X, psi_Y, psi_X)} }
#'   \item{\code{optSim}}{ List of options used in \code{simMrDAG} \code{list(saveAll, seed)} } }
#' If \code{graph_X = "reverse"} is specified, a DAG model between all exposures and outcomes is simulated allowing for reverse 
#' causation. For each replicate, depending on whether \code{saveAll = TRUE}, the value returned is a list object: 
#' \code{list(Y, X, G, betaHat_Y, seBetaHat_Y, pvalueBetaHat_Y, betaHat_X, seBetaHat_X, pvalueBetaHat_X, Rho_YX, DAG_YX, A_YX, }
#' \code{parSim, optSim)} or if \code{saveAll = FALSE}, \code{list(betaHat_Y, seBetaHat_Y, pvalueBetaHat_Y, betaHat_X, seBetaHat_X, }
#' \code{pvalueBetaHat_X, Rho_YX, DAG_YX, A_YX, parSim, optSim)}. They are similar to the list above, with two exceptions:
#' \describe{
#'   \item{\code{Rho_YX}}{ Matrix (traits times traits) of simulated mediations, direct and reverse effects between all traits 
#'                         with outcomes first followed by the exposures }
#'   \item{\code{DAG_YX}}{ graphNEL object that describes the dependency relations between all traits with outcomes first followed 
#'                         by the exposures} } 
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @examples
#' \dontrun{
#' # Example
#' # In this first example, two-sample summary-level design data are simulated with five outcomes 
#' # and 15 exposures with 100 IVs on 50,000 individuals in both samples. The graphical structure 
#' # within the outcomes and exposures is a topologically ordered DAG with a default probability 
#' # of connecting a node to another node (pDAG_Y = 0.1 and pDAG_X = 0.1) as well as a default 
#' # proportion of simulated direct effects between exposures and outcomes 
#' # (thetaSparsness = 0.2)
#'
#' q <- 5
#' p <- 15
#'
#' simData <- simMrDAG(nRep = 1, q = q, p = p, nSubject_X = 50000, 
#'                     graph_Y = "dir", graph_X = "dir", 
#'                     progress = TRUE, fileName = NULL, filePath = NULL, 
#'                     seed = 28061971)
#'
#' # The indices of the traits that define the outcomes and exposures are provided in a list 
#' # object 
#'
#' MRconstraint <- NULL
#' MRconstraint$Y_idx <- 1 : q
#' MRconstraint$X_idx <- q + (1 : p)
#'
#' # The MrDAG algorithm is run to generate 2,000 posterior samples of all unknowns
#'
#' output <- MrDAG(data = cbind(simData$betaHat_Y, simData$betaHat_X), 
#'                 niter = 150000, burnin = 50000, thin = 50, pp = 0.1, tempMax = 20, 
#'                 MRconstraint = MRconstraint, fileName = NULL, 
#'                 seed = 28061971)
#'
#' # Finally, PEIPs and direct effects are estimated
#'
#' PEIP <- get_postEdgeProb(output)
#' genPredEffects <- get_genPredEffects(output)
#'
#' # compared with the simulated ones
#'
#' round(sum((PEIP$postEdgeProb[q + (1 : p), q + (1 : p)] - simData$A_X) ^2, na.rm = TRUE), 3)
#' round(sum((PEIP$postEdgeProb[1 : q, 1 : q] - simData$A_Y) ^2, na.rm = TRUE), 3)
#' round(sum((genPredEffects$genPredEffects - simData$Theta) ^2), 3)
#'
#' # Example
#' # In this second example, two-sample summary-level design data are simulated with five 
#' # outcomes and 15 exposures with 100 IVs, MAF = 0.2 on 50,000 individuals in both samples. 
#' # The graphical structure within the outcomes and exposures is an Erdős-Rényi graph ("er") 
#' # with an expected number of neighbours per node within the outcomes dDAG_Y = 2 and the 
#' # exposures dDAG_X = 1, respectively. The default proportion of simulated direct effects 
#' # between exposures and outcomes is used (thetaSparsness = 0.2). Individual-level data are 
#' # also saved
#'
#' q <- 5
#' p <- 15
#'
#' simData <- simMrDAG(nRep = 1, q = q, p = p, nSubject_X = 50000, MAF = 0.2, 
#'                     graph_Y = "dir", methodDAG_Y = "er", dDAG_Y = 2, 
#'                     graph_X = "dir", methodDAG_X = "er", dDAG_X = 1, 
#'                     saveAll = TRUE, progress = FALSE, fileName = NULL, 
#'                     seed = 28061971)
#'
#' # The indices of the traits that define the outcomes and exposures are provided in a list 
#' # object 
#'
#' MRconstraint <- NULL
#' MRconstraint$Y_idx <- 1 : q
#' MRconstraint$X_idx <- q + (1 : p)
#'
#' # The MrDAG algorithm is run to generate 2,000 posterior samples of all unknowns
#'
#' output <- MrDAG(data = cbind(simData$betaHat_Y, simData$betaHat_X), 
#'                 niter = 150000, burnin = 50000, thin = 50, pp = 0.1, tempMax = 20, 
#'                 MRconstraint = MRconstraint, fileName = NULL, 
#'                 seed = 28061971)
#'
#' # Finally, PEIPs and direct effects are estimated
#'
#' PEIP <- get_postEdgeProb(output)
#' genPredEffects <- get_genPredEffects(output)
#'
#' # compared with the simulated ones
#'
#' round(sum((PEIP$postEdgeProb[q + (1 : p), q + (1 : p)] - simData$A_X) ^2, na.rm = TRUE), 3)
#' round(sum((PEIP$postEdgeProb[1 : q, 1 : q] - simData$A_Y) ^2, na.rm = TRUE), 3)
#' round(sum((genPredEffects$genPredEffects - simData$Theta) ^2), 3)
#'
#' # Example
#' # In this last example, reverse causation in a two-sample summary-level design data is 
#' # simulated with five outcomes and 15 exposures on 50,000 individuals for the exposures and 
#' # 100,000 for the outcomes. The graphical structure between all traits is an Erdős-Rényi 
#' # graph ("er") with an expected number of neighbours per node between all traits dDAG_X = 3
#'
#' q <- 5
#' p <- 15
#'
#' simData <- simMrDAG(nRep = 1, q = q, p = p, nSubject_X = 50000, ratioSubject_XY = 0.5, 
#'                     graph_X = "reverse", methodDAG_X = "er", dDAG_X = 3, 
#'                     progress = TRUE, fileName = NULL, 
#'                     seed = 28061971)
#'
#'  MRconstraint <- NULL
#'  MRconstraint$Y_idx <- 1 : q
#'  MRconstraint$X_idx <- q + (1 : p)
#'
#'  output <- MrDAG(data = cbind(simData$betaHat_Y, simData$betaHat_X), 
#'                  niter = 150000, burnin = 50000, thin = 50, pp = 0.1, tempMax = 20, 
#'                  MRconstraint = MRconstraint, fileName = NULL, 
#'                  seed = 28061971)
#'
#' # Finally, PEIPs and direct effects are estimated
#'
#' PEIP <- get_postEdgeProb(output)
#' genPredEffects <- get_genPredEffects(output)
#'
#' # compared with the simulated ones
#'
#' round(sum((PEIP$postEdgeProb[q + (1 : p), q + (1 : p)] - simData$A_XY[1 : p, 1 : p]) ^2, 
#'           na.rm = TRUE), 3)
#' round(sum((PEIP$postEdgeProb[1 : q, 1 : q] - simData$A_XY[p + (1 : q), p + (1 : q)]) ^2, 
#'           na.rm = TRUE), 3)
#' round(sum((genPredEffects$genPredEffects - simData$Rho_XY[1 : p, p + (1 : q)]) ^2), 3)
#' 
#' # Results seem not to be influenced to a large extent by the simulated reverse causation
#' 
#' round(sum(simData$A_XY[p + (1 : q), 1 : p]) / (q * p), 3)
#' }

########## simMrDAG ##########

simMrDAG<- function(nRep = 1, q, p, nSubject_X = 50000, ratioSubject_XY = 1, nIV = 100, MAF = 0.05, LDmax = NULL, 
                    thetaRange = c(-2, 2), thetaSparsness = 0.2, 
                    graph_Y = "dir", pDAG_Y = 0.1, methodDAG_Y = NULL, dDAG_Y = NULL, par1DAG_Y = NULL, par2DAG_Y = NULL, 
                    graph_X = "dir", pDAG_X = 0.1, methodDAG_X = NULL, dDAG_X = NULL, par1DAG_X = NULL, par2DAG_X = NULL, 
                    herit_Y = 0.25, betaHatRange_X = c(-2, 2), herit_X = 0.1, 
                    gamma_Y = 1, rho_X = 1, 
                    psi_Y = 1, psi_X = 2, 
                    saveAll = FALSE, fileName = "Sim_MrDAG_object", filePath = NULL, progress = TRUE, 
                    seed = 31122021)
{
  if (graph_Y == "dir" && q < 3)
    stop("When graph_Y = dir, at least three outcomes are needed (q >= 3)")
  
  cat("##\n")
  cat("## Simulator for MrDAGSuite\n")
  cat("## Version: 0.1.0\n")
  cat("##\n")
  cat("## Copyright (C) 2026 L. Bottolo and V. Zuber\n")
  cat("##\n\n")
  
  cat("##\n")
  cat("## ## The parameters of the simulation are:\n")
  cat("## ## Number of replicates:", nRep, "\n")
  cat("## ## Number of outcomes:", q, "\n")
  cat("## ## Number of exposures:", p, "\n")
  
  if (is.null(ratioSubject_XY))
  {
    cat("## ## Number of individuals one-sample design:", nSubject_X, "\n")
  } else {
    cat("## ## Number of individuals outcomes Y in two-sample design:", round(nSubject_X / ratioSubject_XY), "\n")
    cat("## ## Number of individuals exposures X in two-sample design:", nSubject_X, "\n")
  }
  
  cat("## ## Number of instrumental variables (IVs):", nIV, "\n")
  
  if (length(MAF) == 1)
  {
    cat("## ## Minor Allele Frequency:", MAF, "\n")
  } else {
    cat("## ## Minor Allele Frequency average:", MAF[1], "\n")
    cat("## ## Minor Allele Frequency standard deviation:", MAF[2], "\n")
  }
  
  if (!is.null(LDmax))
  {
    cat("## ## Max pairwise LD:", LDmax, "\n")
  } else {
    cat("## ## Independent genetic variants simulated\n")
  }
  
  if (graph_X != "reverse")
  {
    
    if (graph_Y == "dir")
    {
      
      if (is.null(methodDAG_Y))
      {
        cat("## ## Type of graph simulated within outcomes Y: Directed\n")
        cat("## ## Type of random graph model within outcomes Y: Topologically ordered\n")
        cat("## ## Probability of connecting a node to another node within outcomes Y:", round(pDAG_Y, digits = 2), "\n")
      } else {
        
        if (methodDAG_Y == "bipartite")
        {
          
          if (is.null(par1DAG_Y))
          {
            par1DAG_Y <- 0.5
          }
        
        } else if (methodDAG_Y == "barabasi") {
         
         if (is.null(par1DAG_Y))
         {
           par1DAG_Y <- 1
         }
         
        } else if (methodDAG_Y == "geometric") {
          
          if (is.null(par1DAG_Y))
          {
            par1DAG_Y <- 2
          }
           
        } else if (methodDAG_Y == "interEr") {
          
          if (is.null(par1DAG_Y))
          {
            par1DAG_Y <- 2
            par1DAG_Y <- 0.25
          }
        
        }
        
        cat("## ## Type of graph simulated within outcomes Y: Directed\n")
        cat("## ## Type of random graph model within outcomes Y:", methodDAG_Y, "\n")
        cat("## ## Expected number of neighbours per node within the outcomes Y:", dDAG_Y, "\n")
        cat("## ## Parameter 1 of random graph model within outcomes Y:", par1DAG_Y, "\n")
        cat("## ## Parameter 2 of random graph model within outcomes Y:", par2DAG_Y, "\n")
      }
     
    } else {
      cat("## ## Type of graph simulated within outcomes Y: Indirected\n")
    }
    
    if (graph_X == "dir")
    {
      
      if (is.null(methodDAG_X))
      {
        cat("## ## Type of graph simulated within exposures X: Directed\n")
        cat("## ## Type of random graph model within outcomes X: Topologically ordered\n")
        cat("## ## Probability of connecting a node to another node within exposures X:", round(pDAG_X, digits = 2), "\n")
      } else {
        
        if (methodDAG_X == "bipartite")
        {
          
          if (is.null(par1DAG_X))
          {
            par1DAG_X <- 0.5
          }
        
        } else if (methodDAG_X == "barabasi") {
         
         if (is.null(par1DAG_X))
         {
           par1DAG_X <- 1
         }
         
        } else if (methodDAG_X == "geometric") {
          
          if (is.null(par1DAG_X))
          {
            par1DAG_X <- 2
          }
           
        } else if (methodDAG_X == "interEr") {
          
          if (is.null(par1DAG_X))
          {
            par1DAG_X <- 2
            par1DAG_X <- 0.25
          }
        
        }
        
        cat("## ## Type of graph simulated within exposures X: Directed\n")
        cat("## ## Type of random graph model within exposures X:", methodDAG_X, "\n")
        cat("## ## Expected number of neighbours per node within the exposures X:", dDAG_X, "\n")
        cat("## ## Parameter 1 of random graph model within exposures X:", par1DAG_X, "\n")
        cat("## ## Parameter 2 of random graph model within exposures X:", par2DAG_X, "\n")
      }
    
    } else {
      cat("## ## Type of graph simulated within exposures X: Indirected\n")
    }
    
    cat("## ## Range of direct effects:", thetaRange, "\n")
    cat("## ## Level of sparsity of the direct effects:", thetaSparsness, "\n")
    cat("## ## Proportion of variance explained for each response Y:", herit_Y, "\n")
    cat("## ## Proportion of variance explained for each exposure X:", herit_X, "\n")
    cat("## ## Range of beta X:", betaHatRange_X, "\n")
    cat("## ## Average level of mediation effects within outcomes Y:", gamma_Y, "\n")
    cat("## ## Average level of mediation effects within exposure X:", rho_X, "\n")
    cat("## ## Level of confounding effects on outcomes Y:", psi_Y, "\n")
    cat("## ## Level of confounding effects on exposures X:", psi_X, "\n")
    cat("##\n\n")
  } else {
    
    if (is.null(methodDAG_X))
    {
      cat("## ## Type of graph simulated between all traits: Directed with reverse causation\n")
      cat("## ## Type of random graph model between all traits: Topologically ordered\n")
      cat("## ## Probability of connecting a node to another node between all traits:", round(pDAG_X, digits = 2), "\n")
      cat("## ## Range of direct effects:", thetaRange, "\n")
    } else {
      cat("## ## Type of graph simulated between all trait: Directed with reverse causation\n")
      cat("## ## Type of random graph model between all traits:", methodDAG_X, "\n")
      cat("## ## Expected number of neighbours per node between all traits:", dDAG_X, "\n")
      cat("## ## Parameter 1 of random graph model between all traits:", par1DAG_X, "\n")
      cat("## ## Parameter 2 of random graph model between all traits:", par2DAG_X, "\n")
      cat("## ## Range of direct and reverse effects:", thetaRange, "\n")
    }
    
    cat("## ## Proportion of variance explained for each response Y:", herit_Y, "\n")
    cat("## ## Proportion of variance explained for each exposure X:", herit_X, "\n")
    cat("## ## Average level of mediation effects within outcomes Y:", gamma_Y, "\n")
    cat("## ## Average level of mediation effects within exposure X:", rho_X, "\n")
    cat("## ## Level of confounding effects on outcomes Y:", psi_Y, "\n")
    cat("## ## Level of confounding effects on exposures X:", psi_X, "\n")
    cat("##\n\n")
  }
  
  digits <- 3
  startTime <- Sys.time()
  
  if (is.null(ratioSubject_XY))
  {
    nSubject_Y <- nSubject_X
    nSubject <- nSubject_X
    idx_subject_X <- 1 : nSubject_X
    idx_subject_Y <- 1 : nSubject_X
  } else {
    nSubject_Y <- round(ratioSubject_XY * nSubject_X)
    nSubject <- nSubject_X + nSubject_Y
    idx_subject_X <- 1 : nSubject_X
    idx_subject_Y <- (nSubject_X + 1) : nSubject
  }
  
  replist <- seq(1, nRep)
  
  ngamma_Y <- length(gamma_Y)
  gamma_Ylist <- gamma_Y
  nrho_X <- length(rho_X)
  rho_Xlist <- rho_X
  
  npsi_Y <- length(psi_Y)
  psi_Ylist <- psi_Y
  npsi_X <- length(psi_X)
  psi_Xlist <- psi_X
  
  idx_gamma_Y <- 1
  
  for (gamma_Y in gamma_Ylist)
  {
    idx_rho_X <- 1
    
    for (rho_X in rho_Xlist)
    {
      idx_psi_Y <- 1
      
      for (psi_Y in psi_Ylist)
      {
        idx_psi_X <- 1
        
        for (psi_X in psi_Xlist)
        {
          
          idx_rep <- 1
          
          for (rep in replist)
          {
            
            if (progress == TRUE)
            {
              cat("## Replicate:", rep, "\n")
            } else {
              cat("## Replicate:", rep, "\n\n")
            }
            
            ########## Storing ##########
            
            Y <- matrix(NA, nSubject_Y, q, dimnames = list(paste0("subject", seq(1, nSubject_Y)), 
                                                      paste0("Y", seq(1, q))))
            X <- matrix(NA, nSubject_X, p, dimnames = list(paste0("subject", seq(1, nSubject_X)), 
                                                      paste0("X", seq(1, p))))
            
            betaHat_Y <- matrix(NA, nIV, q, dimnames = list(paste0("IV", seq(1, nIV)), 
                                                            paste0("betaHat_Y", seq(1, q))))
            seBetaHat_Y <- matrix(NA, nIV, q, dimnames = list(paste0("IV", seq(1, nIV)), 
                                                              paste0("seBetaHat_Y", seq(1, q))))
            pvalueBetaHat_Y <- matrix(NA, nIV, q, dimnames = list(paste0("IV", seq(1, nIV)), 
                                                                  paste0("pvalueBetaHat_Y", seq(1, q))))
            betaHat_X <- matrix(NA, nIV, p, dimnames = list(paste0("IV", seq(1, nIV)), 
                                                            paste0("betaHat_X", seq(1, p))))
            seBetaHat_X <- matrix(NA, nIV, p, dimnames = list(paste0("IV", seq(1, nIV)), 
                                                              paste0("seBetaHat_X", seq(1, p))))
            pvalueBetaHat_X <- matrix(NA, nIV, p, dimnames = list(paste0("IV", seq(1, nIV)), 
                                                                  paste0("pvalueBetaHat_X", seq(1, p))))
            
            set.seed(seed + rep + round(rho_X + gamma_Y + psi_X + psi_Y))
            
            if (length(MAF) == 1 | !is.null(LDmax))
            {
              MAF <- rep(MAF, nIV)
            } else {
              MAF <- abs(stats::rnorm(n = nIV, mean = MAF[1], sd = MAF[2]))
              MAF[MAF >= 0.5] <- (1 - MAF[MAF >= 0.5])
              MAF[MAF < 0] <- abs(MAF[MAF < 0])
            }
            
            if (is.null(LDmax))
            {
              G <- matrix(NA, nrow = nSubject, ncol = nIV)
              
              for (IV in 1 : nIV)
              {
                G[, IV] <- stats::rbinom(nSubject, prob = MAF[IV], size = 2)
              }
            
            } else {
              resultMAF <- 1
              cost <- 0.5 
              betaCoeff <- 1
              xdata <- data.frame(rep(stats::rnorm(nSubject), each = nIV))
              colnames(xdata) <- "x"
              latentCorrelationMatrix <- kronecker(LDmax ^ (abs(outer(1 : nIV, 1 : nIV, "-"))), diag(3))
              
              while (resultMAF > MAF[1])
              {
                betaIntercept <- c(2, 1) + cost
                sim <- SimCorMultRes::rmult.acl(clsize = nIV, intercepts = betaIntercept, 
                                                betas = betaCoeff, xformula = ~ x, xdata = xdata, cor.matrix = latentCorrelationMatrix)
                table <- table(sim$Ysim - 1)
                resultMAF <- (table[2] + 2 * table[3]) / (nSubject * nIV)
                cost <- cost + 0.1
                
                cat("Simulating matrix genotypes G with LD structure - Average MAF:", round(resultMAF, digits = 2), "\n")
              }
              
              G <- sim$Ysim - 1
              rm(xdata)
              rm(latentCorrelationMatrix)
              rm(sim)
              cat("\n")
            }
            
            U <- stats::rnorm(n = nSubject, mean = 0, sd = 1)
            
            if (graph_X != "reverse")
            {
              
              if (graph_X == "undir")
              {
                pDAG_X <- NULL
                dDAG_X = NULL
                methodDAG_X = NULL
                par1DAG_X = NULL
                par2DAG_X = NULL
                
                if (rho_X >= 1)
                {
                  stop("When graph_X = undir, rho_X < 1")
                }
                
                Rho_X <- solve(rho_X ^ (abs(outer(1 : p, 1 : p, "-"))))
                A_X <- Rho_X
                diag(A_X) <- 0
                A_X[abs(A_X) < 10 ^ (-5)] <- 0
                A_X[A_X != 0] <- 1
                idx_ts_X <- seq(1, p)
                rank_ts_X <- seq(1, p)
                rownames(A_X)  <- paste0("expsosureX", seq(1, p))
                colnames(A_X)  <- paste0("expsosureX", seq(1, p))
                DAG_X <- methods::as(A_X, "graphNEL")
                A_X <- A_X[idx_ts_X, idx_ts_X]
                
                SigmaErrorBeta_X <- (betaHatRange_X[2] - betaHatRange_X[1]) / 4 * diag(p)
                SigmaErrorBeta_X <- chol(Rho_X) %*% SigmaErrorBeta_X %*% t(chol(Rho_X))
              } else {
                                                 
                if (is.null(methodDAG_X))
                {
                  DAG_X <- pcalg::randomDAG(p, prob = pDAG_X)
                  idx_ts_X <- seq(1, p)
                  rank_ts_X <- seq(1, p)
                } else {
                  DAG_X <- pcalg::randDAG(p, d = dDAG_X, method = methodDAG_X, par1 = par1DAG_X, par2 = par2DAG_X, weighted = FALSE)
                  idx_ts_X <- as.numeric(RBGL::tsort(DAG_X))
                  rank_ts_X <- sort(idx_ts_X, index.return = TRUE)$ix
                }
                
                A_X <- methods::as(DAG_X, "matrix")
                A_X[A_X != 0] <- 1
                rownames(A_X)  <- paste0("expsosureX", seq(1, p))
                colnames(A_X)  <- paste0("expsosureX", seq(1, p))
                DAG_X <- methods::as(A_X, "graphNEL")
                A_X <- A_X[idx_ts_X, idx_ts_X]
                
                Rho_X <- A_X * matrix(stats::runif(p * p, rho_X * 0.9, rho_X  * 1.1), p, p) * sample(c(-1, 1), size = p * p, replace = TRUE)
                diag(Rho_X) <- 1
                SigmaErrorBeta_X <- (betaHatRange_X[2] - betaHatRange_X[1]) / 4 * diag(p)
                SigmaErrorBeta_X <- solve(t(Rho_X)) %*% SigmaErrorBeta_X %*% solve(Rho_X)
                Rho_X <- -Rho_X   # ***
                diag(Rho_X) <- 1   # ***
              }
              
              beta_X <- MASS::mvrnorm(n = nIV, mu = rep(0, times = p), Sigma = SigmaErrorBeta_X)
              
              SigmaError_X <- diag(diag(stats::var(G %*% beta_X + matrix(psi_X * U, nrow = nSubject, ncol = p)))) * (1 - herit_X) / herit_X
              error_X <- MASS::mvrnorm(n = nSubject, mu = rep(0, times = p), Sigma = SigmaError_X)  
              X <- G %*% beta_X + matrix(psi_X * U, nrow = nSubject, ncol = p) + error_X
              
              Theta <- MASS::mvrnorm(n = p, mu = rep(0, times = q), Sigma = (thetaRange[2] - thetaRange[1]) / 4 * diag(q))
              idx_Theta <- matrix(sample(c(rep(1, floor(p * q * thetaSparsness)), 
                                  rep(0, times = p * q - floor(p * q * thetaSparsness)))), nrow = p, ncol = q)
              Theta <- Theta * idx_Theta
              
              if (graph_Y == "undir")
              {
                pDAG_Y <- NULL
                dDAG_Y = NULL
                methodDAG_Y = NULL
                par1DAG_Y = NULL
                par2DAG_Y = NULL
                                       
                if (gamma_Y >= 1)
                {
                  stop("When graph_Y = undir, gamma_Y < 1")
                }
                
                Gamma_Y <- solve(gamma_Y ^ (abs(outer(1 : q, 1 : q, "-"))))
                A_Y <- Gamma_Y
                diag(A_Y) <- 0
                A_Y[abs(A_Y) < 10 ^ (-5)] <- 0
                A_Y[A_Y != 0] <- 1
                idx_ts_Y <- seq(1, q)
                rank_ts_Y <- seq(1, q)
                rownames(A_Y)  <- paste0("outcomesY", seq(1, q))
                colnames(A_Y)  <- paste0("outcomesY", seq(1, q))
                DAG_Y <- methods::as(A_Y, "graphNEL")
                A_Y <- A_Y[idx_ts_Y, idx_ts_Y]
                
                SigmaError_Y <- diag(diag(stats::var(X %*% Theta + matrix(psi_Y * U, nrow = nSubject, ncol = q)))) * (1 - herit_Y) / herit_Y
                SigmaError_Y <- chol(Gamma_Y) %*% SigmaError_Y %*% t(chol(Gamma_Y))
                error_Y <- MASS::mvrnorm(n = nSubject, mu = rep(0, times = q), Sigma = SigmaError_Y)
                
                Y <- X %*% Theta + matrix(psi_Y * U, nrow = nSubject, ncol = q) + error_Y
              } else {
                A_Y <- NULL
                
                while (sum(A_Y) <= 1)
                {
                  
                  if (is.null(methodDAG_Y))
                  {
                    DAG_Y <- pcalg::randomDAG(q, prob = pDAG_Y)
                    idx_ts_Y <- seq(1, q)
                    rank_ts_Y <- seq(1, q)
                  } else {
                    DAG_Y <- pcalg::randDAG(q, d = dDAG_Y, method = methodDAG_Y, par1 = par1DAG_Y, par2 = par2DAG_Y, weighted = FALSE)
                    idx_ts_Y <- as.numeric(RBGL::tsort(DAG_Y))
                    rank_ts_Y <- sort(idx_ts_Y, index.return = TRUE)$ix
                  }
                  
                  A_Y <- methods::as(DAG_Y, "matrix")
                  A_Y[A_Y != 0] <- 1
                }
                
                rownames(A_Y)  <- paste0("outcomeY", seq(1, q))
                colnames(A_Y)  <- paste0("outcomeY", seq(1, q))
                DAG_Y <- methods::as(A_Y, "graphNEL")
                A_Y <- A_Y[idx_ts_Y, idx_ts_Y]
                
                Gamma_Y <- A_Y * matrix(stats::runif(q * q, gamma_Y * 0.9, gamma_Y  * 1.1), q, q) * sample(c(-1, 1), size = q * q, replace = TRUE)
                diag(Gamma_Y) <- 1
                SigmaError_Y <- diag(diag(stats::var(X %*% Theta + matrix(psi_Y * U, nrow = nSubject, ncol = q)))) * (1 - herit_Y) / herit_Y
                error_Y <- MASS::mvrnorm(n = nSubject, mu = rep(0, times = q), Sigma = SigmaError_Y)
                
                Y <- X %*% Theta + matrix(psi_Y * U, nrow = nSubject, ncol = q) + error_Y
                
                for (k in 2 : q)
                {
                    Y[, k] <- Y[, k] + Y[, 1 : (k - 1), drop = FALSE] %*% Gamma_Y[1 : (k - 1), k]
                }
              
              }
              
              if (progress == TRUE)
              {
                cat("## No. simulated mediation effects X:", sum(A_X), "\n")
                cat("## No. simulated mediation effects Y:", sum(A_Y), "\n")
                cat("## No. simulated direct effects:", sum((abs(Theta) > 0) * 1), "\n")
                cat("## Prop. simulated mediation effects X:", round(sum(A_X) / (p * (p - 1)), digits = 2), "\n")
                cat("## Prop. simulated mediation effects Y:", round(sum(A_Y) / (q * (q - 1)), digits = 2), "\n")
                cat("## Prop. simulated direct effects:", round(sum((abs(Theta) > 0) * 1) / (p * q), digits = 2), "\n")
              }
              
              Y <- Y[, rank_ts_Y]
              X <- X[, rank_ts_X]
              A_Y <- A_Y[rank_ts_Y, rank_ts_Y]
              A_X <- A_X[rank_ts_X, rank_ts_X]
              Gamma_Y <- Gamma_Y[rank_ts_Y, rank_ts_Y]
              Rho_X <- Rho_X[rank_ts_X, rank_ts_X]
              Theta <- Theta[rank_ts_X, rank_ts_Y]
            } else {
              graph_Y <- NULL
              pDAG_Y = NULL
              dDAG_Y = NULL
              methodDAG_Y <- NULL
              par1DAG_Y <- NULL
              par2DAG_Y <- NULL 
              
              if (is.null(methodDAG_X))
              {
                DAG_XY <- pcalg::randomDAG(p + q, prob = pDAG_X, V = c(paste0("exposureX", seq(1, p)), paste0("outcomeY", seq(1, q))))
                idx_ts_XY <- seq(1, p + q)
                rank_ts_XY <- seq(1, p + q)
              } else {
                DAG_XY <- pcalg::randDAG(p + q, d = dDAG_X, method = methodDAG_X, par1 = par1DAG_X, par2 = par2DAG_X, weighted = FALSE)
                idx_ts_XY <- as.numeric(RBGL::tsort(DAG_XY))
                rank_ts_XY <- sort(idx_ts_XY, index.return = TRUE)$ix
              }
              
              A_XY <- methods::as(DAG_XY, "matrix")
              rownames(A_XY) <- c(paste0("exposureX", seq(1, p)), paste0("outcomeY", seq(1, q)))
              colnames(A_XY) <- c(paste0("exposureX", seq(1, p)), paste0("outcomeY", seq(1, q)))
              DAG_XY <- methods::as(A_XY, "graphNEL")
              A_XY <- A_XY[idx_ts_XY, idx_ts_XY]
              A_XY[A_XY != 0] <- 1
              
              Rho_XY <- matrix(0, (p + q), (p + q))
              Rho_XY[1 : p, 1 : p] <- A_XY[1 : p, 1 : p] * matrix(stats::runif(p * p, rho_X * 0.9, rho_X  * 1.1), p, p) * sample(c(-1, 1), size = p * p, replace = TRUE)
              Rho_XY[p + (1 : q), p + (1 : q)] <- A_XY[p + (1 : q), p + (1 : q)] * matrix(stats::runif(q * q, gamma_Y * 0.9, gamma_Y  * 1.1), q, q) * sample(c(-1, 1), size = q * q, replace = TRUE)
              Rho_XY[1 : p, p + (1 : q)] <- A_XY[1 : p, p + (1 : q)] * matrix(MASS::mvrnorm(n = p, mu = rep(0, times = q), Sigma = (thetaRange[2] - thetaRange[1]) / 4 * diag(q)), p, q) * sample(c(-1, 1), size = p * q, replace = TRUE)
              Rho_XY[p + (1 : q), 1 : p] <- A_XY[p + (1 : q), 1 : p] * matrix(MASS::mvrnorm(n = q, mu = rep(0, times = p), Sigma = (thetaRange[2] - thetaRange[1]) / 4 * diag(p)), q, p) * sample(c(-1, 1), size = q * p, replace = TRUE)
              diag(Rho_XY) <- 1
              
              SigmaErrorBeta_XY <- (betaHatRange_X[2] - betaHatRange_X[1]) / 4 * diag(p + q)
              SigmaErrorBeta_XY <- solve(t(Rho_XY)) %*% SigmaErrorBeta_XY %*% solve(Rho_XY)
              beta_XY <- MASS::mvrnorm(n = nIV, mu = rep(0, times = p + q), Sigma = SigmaErrorBeta_XY)
              
              SigmaError_XY <- diag(diag(stats::var(G %*% beta_XY + matrix(psi_X * U, nrow = nSubject, ncol = p + q))))
              SigmaError_XY[1 : p, 1 : p] <- SigmaError_XY[1 : p, 1 : p] * (1 - herit_X) / herit_X
              SigmaError_XY[p + (1 : q), p + (1 : q)] <- SigmaError_XY[p + (1 : q), p + (1 : q)] * (1 - herit_Y) / herit_Y
              
              error_XY <- MASS::mvrnorm(n = nSubject, mu = rep(0, times = p + q), Sigma = SigmaError_XY)
              XY <- G %*% beta_XY + cbind(matrix(psi_X * U, nrow = nSubject, ncol = p),
                                          matrix(psi_Y * U, nrow = nSubject, ncol = q)) + error_XY
              XY <- XY[, rank_ts_XY]
              X <- XY[, 1 : p]
              Y <- XY[, p + (1 : q)]
              A_XY <- A_XY[rank_ts_XY, rank_ts_XY]
              A_YX <- A_XY[(p + q) : 1, (p + q) : 1]
              Rho_XY <- Rho_XY[rank_ts_XY, rank_ts_XY]
              Rho_XY <- -Rho_XY   # ***
              diag(Rho_XY) <- 1   # ***
              Rho_YX <- Rho_XY[(p + q) : 1, (p + q) : 1]
              
              if (progress == TRUE)
              {
                cat("## No. simulated mediation effects X:", sum(A_XY[1 : p, 1 : p]), "\n")
                cat("## No. simulated mediation effects Y:", sum(A_XY[p + (1 : q), p + (1 : q)]), "\n")
                cat("## No. simulated direct effects:", sum(A_XY[1 : p, p + (1 : q)]), "\n")
                cat("## No. simulated reverse effects:", sum(A_XY[p + (1 : q), 1 : p]), "\n")
                cat("## Prop. simulated mediation effects X:", round(sum(A_XY[1 : p, 1 : p]) / (p * (p - 1)), digits = 2), "\n")
                cat("## Prop. simulated mediation effects Y:", round(sum(A_XY[p + (1 : q), p + (1 : q)]) / (q * (q - 1)), digits = 2), "\n")
                cat("## Prop. simulated direct effects:", round(sum(A_XY[1 : p, p + (1 : q)]) / (p * q), digits = 2), "\n")
                cat("## Prop. simulated reverse effects:", round(sum(A_XY[p + (1 : q), 1 : p]) / (q * p), digits = 2), "\n")
              }
            
            }
            
            for (j in 1 : p)
            {
              
              for (i in 1 : nIV)
              {
                summary_X <- summary(stats::lm(X[idx_subject_X, j] ~ G[idx_subject_X, i]))
                betaHat_X[i, j] <- summary_X$coefficients[2, 1]
                seBetaHat_X[i, j] <- summary_X$coefficients[2, 2]
                pvalueBetaHat_X[i, j] <- summary_X$coefficients[2, 4]
              }
              
              if (progress == TRUE)
              {
                cat("## Simulating exposure X:", j, "\n")
              }
              
            }
            
            for (k in 1 : q)
            {
              
              for (i in 1 : nIV)
              {
                summary_Y <- summary(stats::lm(Y[idx_subject_Y, k] ~ G[idx_subject_Y, i]))
                betaHat_Y[i, k] <- summary_Y$coefficients[2, 1]
                seBetaHat_Y[i, k] <- summary_Y$coefficients[2, 2]
                pvalueBetaHat_Y[i, k] <- summary_Y$coefficients[2, 4]
              }
              
              if (progress == TRUE & k < q)
              {
                cat("## Simulating outcome Y:", k, "\n")
              } else if (progress == TRUE & k == q) {
                cat("## Simulating outcome Y:", k, "\n\n")
              }
            
            }
            
            Y <- Y[idx_subject_Y, ]
            X <- X[idx_subject_X, ]
            
            if (is.null(ratioSubject_XY))
            {
              simMAF <- matrix(NA, 3, nIV)
              
              for (IV in 1 : nIV)
              {
                simMAF[, IV] <- table(factor(G[, IV], levels = 0 : 2))
              }
              
              simMAF <- (simMAF[2, ] + simMAF[3, ] * 2) / (nSubject_X * 2)
              rownames(G) <- paste0("subject", seq(1, nSubject_Y))
              colnames(G) <- paste0("outcomeY", seq(1, nIV))
              rownames(Y) <- paste0("subject", seq(1, nSubject_Y))
              colnames(Y) <- paste0("outcomeY", seq(1, q))
              rownames(X) <- paste0("subject", seq(1, nSubject_Y))
              colnames(X) <- paste0("exposureX", seq(1, p))
            } else {
              G_Y <- G[idx_subject_Y, ]
              G_X <- G[idx_subject_X, ]
              simMAF_Y <- matrix(NA, 3, nIV)
              simMAF_X <- matrix(NA, 3, nIV)
              
              for (IV in 1 : nIV)
              {
                simMAF_Y[, IV] <- table(factor(G_Y[, IV], levels = 0 : 2))
                simMAF_X[, IV] <- table(factor(G_X[, IV], levels = 0 : 2))
              }
              
              simMAF_Y <- (simMAF_Y[2, ] + simMAF_Y[3, ] * 2) / (nSubject_Y * 2)
              simMAF_X <- (simMAF_X[2, ] + simMAF_X[3, ] * 2) / (nSubject_X * 2)
              rownames(G_Y) <- paste0("subjectY", seq(1, nSubject_Y))
              colnames(G_Y) <- paste0("genotype", seq(1, nIV))
              rownames(G_X) <- paste0("subjectX", seq(1, nSubject_X))
              colnames(G_X) <- paste0("genotype", seq(1, nIV))
              rownames(Y) <- paste0("subjectY", seq(1, nSubject_Y))
              colnames(Y) <- paste0("outcomeY", seq(1, q))
              rownames(X) <- paste0("subjectX", seq(1, nSubject_X))
              colnames(X) <- paste0("exposureX", seq(1, p))
            }
            
            if (graph_X != "reverse")
            {
              rownames(Theta) <- paste0("exposureX", seq(1, p))
              colnames(Theta) <- paste0("outcomeY", seq(1, q))
              
              rownames(A_Y) <- paste0("outcomeY", seq(1, q))
              colnames(A_Y) <- paste0("outcomeY", seq(1, q))
              rownames(Gamma_Y) <- paste0("outcomeY", seq(1, q))
              colnames(Gamma_Y) <- paste0("outcomeY", seq(1, q))
              
              rownames(A_X) <- paste0("exposureX", seq(1, p))
              colnames(A_X) <- paste0("exposureX", seq(1, p))
              rownames(Rho_X) <- paste0("exposureX", seq(1, p))
              colnames(Rho_X) <- paste0("exposureX", seq(1, p))
              
              parSim <- list(nRep = nRep, q = q, p = p, nSubject_Y = nSubject_Y, nSubject_X = nSubject_X, ratioSubject_XY = ratioSubject_XY, nIV = nIV, 
                             MAF = round(unique(MAF), digits = 2), LDmax = LDmax, 
                             thetaRange = thetaRange, thetaSparsness = thetaSparsness, 
                             graph_Y = graph_Y, pDAG_Y = pDAG_Y, dDAG_Y = dDAG_Y, methodDAG_Y = methodDAG_Y, par1DAG_Y = par1DAG_Y, par2DAG_Y = par2DAG_Y, 
                             graph_X = graph_X, pDAG_X = pDAG_X, dDAG_X = dDAG_X, methodDAG_X = methodDAG_X, par1DAG_X = par1DAG_X, par2DAG_X = par2DAG_X, 
                             herit_Y = herit_Y, herit_X = herit_X, betaHatRange_X = betaHatRange_X, 
                             gamma_Y = gamma_Y, rho_X = rho_X, 
                             psi_Y = psi_Y, psi_X = psi_X)
              optSim <- list(saveAll = saveAll, seed = seed)
              
              if (is.null(ratioSubject_XY))
              {
                
                if (saveAll == TRUE)
                {
                  output <- list(Y = Y, X = X, G = G, simMAF = simMAF, 
                                 betaHat_Y = betaHat_Y, seBetaHat_Y = seBetaHat_Y, pvalueBetaHat_Y = pvalueBetaHat_Y, 
                                 betaHat_X = betaHat_X, seBetaHat_X = seBetaHat_X, pvalueBetaHat_X = pvalueBetaHat_X, 
                                 Theta = Theta, Gamma_Y = Gamma_Y, Rho_X = Rho_X, 
                                 DAG_Y = DAG_Y, DAG_X = DAG_X, 
                                 A_Y = A_Y, A_X = A_X, 
                                 parSim = parSim, optSim = optSim)
                } else {
                  output <- list(betaHat_Y = betaHat_Y, seBetaHat_Y = seBetaHat_Y, pvalueBetaHat_Y = pvalueBetaHat_Y, 
                                 betaHat_X = betaHat_X, seBetaHat_X = seBetaHat_X, pvalueBetaHat_X = pvalueBetaHat_X, 
                                 Theta = Theta, Gamma_Y = Gamma_Y, Rho_X = Rho_X, 
                                 DAG_Y = DAG_Y, DAG_X = DAG_X, 
                                 A_Y = A_Y, A_X = A_X, 
                                 parSim = parSim, optSim = optSim)
                }
              
              } else {
                
                if (saveAll == TRUE)
                {
                  output <- list(Y = Y, X = X, G_Y = G_Y, G_X = G_X, simMAF_Y = simMAF_Y, simMAF_X = simMAF_X, 
                                 betaHat_Y = betaHat_Y, seBetaHat_Y = seBetaHat_Y, pvalueBetaHat_Y = pvalueBetaHat_Y, 
                                 betaHat_X = betaHat_X, seBetaHat_X = seBetaHat_X, pvalueBetaHat_X = pvalueBetaHat_X, 
                                 Theta = Theta, Gamma_Y = Gamma_Y, Rho_X = Rho_X, 
                                 DAG_Y = DAG_Y, DAG_X = DAG_X, 
                                 A_Y = A_Y, A_X = A_X, 
                                 parSim = parSim, optSim = optSim)
                } else {
                  output <- list(betaHat_Y = betaHat_Y, seBetaHat_Y = seBetaHat_Y, pvalueBetaHat_Y = pvalueBetaHat_Y, 
                                 betaHat_X = betaHat_X, seBetaHat_X = seBetaHat_X, pvalueBetaHat_X = pvalueBetaHat_X, 
                                 Theta = Theta, Gamma_Y = Gamma_Y, Rho_X = Rho_X, 
                                 DAG_Y = DAG_Y, DAG_X = DAG_X, 
                                 A_Y = A_Y, A_X = A_X, 
                                 parSim = parSim, optSim = optSim)
                }
              
              }
              
              fileNameTMP <- paste(fileName, "rep", rep, "graph_Y", graph_Y, "graph_X", graph_X, 
                                                         "thetaRange", paste(thetaRange[1], thetaRange[2], sep = ","), "gamma_Y", gamma_Y, "rho_X", rho_X, 
                                                         "psi_Y", psi_Y, "psi_X", psi_X, sep = "_")
            
            } else {
              rownames(A_XY) <- c(paste0("exposureX", seq(1, p)), paste0("outcomeY", seq(1, q)))
              colnames(A_XY) <- c(paste0("exposureX", seq(1, p)), paste0("outcomeY", seq(1, q)))
              rownames(Rho_XY) <- c(paste0("exposureX", seq(1, p)), paste0("outcomeY", seq(1, q)))
              colnames(Rho_XY) <- c(paste0("exposureX", seq(1, p)), paste0("outcomeY", seq(1, q)))
              
              parSim <- list(nRep = nRep, q = q, p = p, nSubject_Y = nSubject_Y, nSubject_X = nSubject_X, ratioSubject_XY = ratioSubject_XY, nIV = nIV, 
                             MAF = unique(MAF), LDmax = LDmax, 
                             graph_XY = graph_X, pDAG_YX = pDAG_X, dDAG_YX = dDAG_X, methodDAG_YX = methodDAG_X, par1DAG_YX = par1DAG_X, par2DAG_YX = par2DAG_X, 
                             herit_Y = herit_Y, herit_X = herit_X, 
                             gamma_Y = gamma_Y, rho_X = rho_X, 
                             psi_Y = psi_Y, psi_X = psi_X)
              optSim <- list(saveAll = saveAll, seed = seed)
              
              if (is.null(ratioSubject_XY))
              {
                
                if (saveAll == TRUE)
                {
                  output <- list(Y = Y, X = X, G = G, 
                                 betaHat_Y = betaHat_Y, seBetaHat_Y = seBetaHat_Y, pvalueBetaHat_Y = pvalueBetaHat_Y, 
                                 betaHat_X = betaHat_X, seBetaHat_X = seBetaHat_X, pvalueBetaHat_X = pvalueBetaHat_X, 
                                 Rho_XY = Rho_XY, 
                                 DAG_XY = DAG_XY, 
                                 A_XY = A_XY, 
                                 parSim = parSim, optSim = optSim)
                } else {
                  output <- list(betaHat_Y = betaHat_Y, seBetaHat_Y = seBetaHat_Y, pvalueBetaHat_Y = pvalueBetaHat_Y, 
                                 betaHat_X = betaHat_X, seBetaHat_X = seBetaHat_X, pvalueBetaHat_X = pvalueBetaHat_X, 
                                 Rho_XY = Rho_XY, 
                                 DAG_XY = DAG_XY, 
                                 A_XY = A_XY, 
                                 parSim = parSim, optSim = optSim)
                }
              
              } else {
                
                if (saveAll == TRUE)
                {
                  output <- list(Y = Y, X = X, G_Y = G_Y, G_X = G_X, 
                                 betaHat_Y = betaHat_Y, seBetaHat_Y = seBetaHat_Y, pvalueBetaHat_Y = pvalueBetaHat_Y, 
                                 betaHat_X = betaHat_X, seBetaHat_X = seBetaHat_X, pvalueBetaHat_X = pvalueBetaHat_X, 
                                 Rho_XY = Rho_XY, 
                                 DAG_XY = DAG_XY, 
                                 A_XY = A_XY, 
                                 parSim = parSim, optSim = optSim)
                } else {
                  output <- list(betaHat_Y = betaHat_Y, seBetaHat_Y = seBetaHat_Y, pvalueBetaHat_Y = pvalueBetaHat_Y, 
                                 betaHat_Y = betaHat_Y, seBetaHat_Y = seBetaHat_Y, pvalueBetaHat_Y = pvalueBetaHat_Y, 
                                 betaHat_X = betaHat_X, seBetaHat_X = seBetaHat_X, pvalueBetaHat_X = pvalueBetaHat_X, 
                                 Rho_XY = Rho_XY, 
                                 DAG_XY = DAG_XY, 
                                 A_XY = A_XY, 
                                 parSim = parSim, optSim = optSim)
                }
              
              }
              
              fileNameTMP <- paste(fileName, "rep", rep, "graph_XY", graph_X,  
                                                         "thetaRange", paste(thetaRange[1], thetaRange[2], sep = ","), "gamma_Y", gamma_Y, "rho_X", rho_X, 
                                                         "psi_Y", psi_Y, "psi_X", psi_X, sep = "_")
            }
            
            if (!is.null(fileName))
            {
              saveas <- paste(fileNameTMP, "RData", sep = ".")
              
              if (is.null(filePath))
              {
                filePath <- getwd()
              } else {
              
                if (!dir.exists(filePath))
                {
                  dir.create(filePath, recursive = TRUE)
                  cat("\nThe directory", paste0("\"", filePath, "\""), "does not exist. It is now created.\n")
                }
              
              }
              
              if (.Platform$OS.type == "windows")
              {
                filePath <- paste0(normalizePath(filePath), "\\")
              } else {
                filePath <- paste0(normalizePath(filePath), "/")
              }
              
              save(output, file = paste0(filePath, saveas))
            }
            
            idx_rep <- idx_rep + 1
          }
          
          idx_psi_X <- idx_psi_X + 1
        }
        
        idx_psi_Y <- idx_psi_Y + 1
      }
      
      idx_rho_X <- idx_rho_X + 1
    }
    
    idx_gamma_Y <- idx_gamma_Y + 1
  }
  
  endTime <- Sys.time()
  timeSimMrDAG <- as.numeric(difftime(time1 = endTime, time2 = startTime, units = "min"))
  timeSimMrDAG <- paste0(round(timeSimMrDAG, digits = 3), "m")
  cat("## Overall simulation time:", timeSimMrDAG, "\n\n")
  
  return(output) 
}