#' This function trains a model given a fixed training data set.
#'
#' @description Given features for a class of sRNA, this function will learn robust
#'              features for the given class of sRNA and return the coefficients
#'              for all the features. Non-robust features will have coefficients
#'              equal to zero.
#'
#' @usage TrainModel(Y, X, bootstrap_iterations = 1000, Kfold = 10)
#'
#' @param Y Class labels of the rows of \emph{X}. There can only be 2 class labels
#' @param X The feature matrix.
#' @param bootstrap_iterations The number of iterations for bootstrapping.
#' @param Kfold An integer specifying \emph{K} fold cross-validation.
#'
#' @return return the coefficients for the learned model. Non-robust features
#'         will be set to zero.
#'
#'
#' @author Carl Tony Fakhry , Kourosh Zarringhalam, Rahul Kulkarni and Ping Chen.
#'
TrainModel <- function(Y, X, bootstrap_iterations = 1000, Kfold = 10){

  size_of_sample = floor(nrow(X)/2)

  mat = cbind(Y,X)
  YONES = mat[mat[,1] == 1, 1, drop = FALSE]
  XONES = mat[mat[,1] == 1, 2:ncol(mat), drop = FALSE]
  YZEROES = mat[mat[,1] == 0, 1, drop = FALSE]
  XZEROES = mat[mat[,1] == 0, 2:ncol(mat), drop = FALSE]

  numG1 = nrow(YONES)
  numG2 = nrow(YZEROES)
  pool = min(numG1,numG2)
  size_of_sample = min(size_of_sample, pool)

  # fit initial model
  fit1 = cv.glmnet(X, Y, family = "binomial", type.measure="deviance", nfolds = Kfold, grouped = F)
  lambda = fit1$lambda
  lambda.min = fit1$lambda.min

  df = as.data.frame(matrix(0, ncol = 514, nrow = bootstrap_iterations))
  rownames(df) = NULL

  for(i in 1:bootstrap_iterations){

    # Random sample and balance data for cross validation
    ones = unique(sort(sample(1:pool, size_of_sample, replace = T)))
    Y1 = YONES[ones,]
    X1 = XONES[ones,]

    zeroes = unique(sort(sample(1:pool, size_of_sample, replace = T)))
    Y2 = YZEROES[zeroes,]
    X2 = XZEROES[zeroes,]

    YY = as.matrix(c(Y1,Y2))
    XX = rbind(X1,X2)

    fit <- glmnet(XX, YY, family= "binomial", lambda = lambda.min)
    nonzeroes = as.vector(coef(fit))
    nonzeroes = nonzeroes[2:515]
    nonzeroes[nonzeroes != 0] = 1
    df[i,] = nonzeroes

  }

  # Get the robust features
  counts <- colSums(df)
  goodvars = which(counts != 0)

  X2 = X[,goodvars,drop=FALSE]

  fit = glmnet(X2, Y, lambda = lambda.min, family= "binomial")
  Coefficients = as.vector(coef(fit))
  Coefficients2 = rep(0, 515)
  Coefficients2[(goodvars + 1)] = Coefficients[2:length(Coefficients)]
  Coefficients2[1] = Coefficients[1] # Set the value of the constant of the regression

  names(Coefficients2) = c("constant", getFeatureNames())

  return(Coefficients2)

}



#' This function trains a class of models given a FASTA file with seed sequences.
#'
#' @description Given a FASTA file of seed sequences of sRNA, we learn a set of models
#'              for prediction on the class of sRNAs which the seed sequences belong
#'              to.
#'
#' @usage TrainModels(seed_sequences_fasta, n_models = 100, bootstrap_iterations = 1000)
#'
#' @param seed_sequences_fasta The FASTA file containing the seed sequences.
#' @param n_models The number of models to be learned. \emph{n_models} must be
#'                 greater than zero.
#' @param bootstrap_iterations The number of iterations for bootstrapping.
#'
#' @return an object of class InvenireSRNA which has a data frame \emph{coefficients_4_features}
#'         containing the coefficients of the Boltzmann Triplet features for
#'         \emph{bootstrap_iterations} number of models trained for prediction.
#'
#' @return An object with S3 class \emph{InvenireSRNA} which contains the following elements:
#' \itemize{
#'   \item{\strong{coefficients_4_features} }{a data frame containing the coefficients of the Boltzmann Triplet features for
#'         \emph{bootstrap_iterations} number of models trained for prediction.}
#' }
#'
#'
#' @author Carl Tony Fakhry , Kourosh Zarringhalam, Rahul Kulkarni and Ping Chen.
#'
#' @export
TrainModels <- function(seed_sequences_fasta, n_models = 100, bootstrap_iterations = 1000){

  # Get the features of the positive sequences and the negative sequences
  positive_features = getBoltzmanTripletFeatures(seed_sequences_fasta)
  minimum_free_energies = getAllMinimumEnergies(seed_sequences_fasta)
  coefs_df = data.frame()

  if(n_models <= 0 || floor(n_models) != n_models) stop("n_models must be a positive integer!")
  if(bootstrap_iterations <= 0 || floor(bootstrap_iterations) != bootstrap_iterations) stop("boostrap_iterations must be a positive integer!")

  Kfold = 10
  for(i in 1:n_models){

    if((i %% n_models) == 10){
      paste(i, "models learned so far!", sep = " ")
    }

    negative_features = getNegativeSequencesFeatures(seed_sequences_fasta, min(minimum_free_energies), max(minimum_free_energies))

    X = rbind(positive_features, negative_features)
    Y = c(rep(1, nrow(positive_features)), rep(0, nrow(negative_features)))
    Coefficients2 = TrainModel(X, Y, bootstrap_iterations = bootstrap_iterations, Kfold = Kfold)
    coefs_df = rbind(coefs_df, Coefficients2)

  }

  names(coefs_df) = c("constant", getFeatureNames())
  lst = list()
  lst[["coefficients_4_features"]] = coefs_df
  class(lst) = "InvenireSRNA"

  return(lst)

}



#' Given a model for prediction, compute the probabilities of the new sequences.
#'
#' @description Given a model for prediction, compute the probabilities of the
#'              new sequences.
#'
#' @usage predict_fasta(fasta_file, InvenireSRNA_model = NULL)
#'
#' @param fasta_file A path to a FASTA file containing the RNA sequences constructed with
#'                   nucleotides A,C,G,T and U.
#' @param InvenireSRNA_model An object with S3 class \emph{InvenireSRNA} to be used for prediction. If \emph{InvenireSRNA_model = NULL}
#'              then the pretrained model to predict CsrA regulating sRNAs is used.
#'
#' @return A data frame containing the mean and standard deviation of the prediction
#'               probabilities using the given \emph{InvenireSRNA_model} object.
#'
#' @author Carl Tony Fakhry , Kourosh Zarringhalam, Rahul Kulkarni and Ping Chen.
#'
#' @export
predict_fasta <- function(fasta_file, InvenireSRNA_model = NULL){

  mean_probabilities = c()
  sd_probabilities = c()
  regions = c()

  # Read in the FASTA file
  fasta = readBStringSet(fasta_file, use.names = T)

  if(is.null(InvenireSRNA_model)){
    ff = paste(system.file(package="InvenireSRNA"), "extdata/CsrA_models.tab", sep="/")
    coefficients_4_features = read.table(ff, header = T)
  }else{
    if(class(InvenireSRNA_model) != "InvenireSRNA") stop("Please provide a model of class InvenireSRNA!")
    coefficients_4_features = InvenireSRNA_model$coefficients_4_features
  }

  for(i in 1:length(fasta)){

    rna = as.character(fasta[[i]])
    rna = tolower(rna)
    rna = gsub("t", "u", rna)

    split_seq <- strsplit(rna, split = "")[[1]]
    if(!all(split_seq %in% c("A","a","U","u","G","g","T","t","C","c"))){
      warning(paste("Prediction skipped for RNA sequence number ", i, " because it cotains letters other than a, c, g, t and u!"))
      next
    }

    RNAsubopt = 'RNAsubopt -p 1000'
    structs = system(RNAsubopt, intern = T,wait = TRUE, input = rna)
    features = get_features(rna, structs)
    features = c(1, features) # constant has feature 1

    prods = as.vector(apply(coefficients_4_features, 1, function(x) sum(x * features)))
    probabilities = sapply(prods, function(x) 1/(1+exp(-x))) # apply the sigmoid for prediction

    mean_probabilities = c(mean_probabilities, mean(probabilities))
    sd_probabilities = c(sd_probabilities, sd(probabilities))
    regions = c(regions, names(fasta)[i])

  }

  df = data.frame(regions = regions, Mean_Probability = mean_probabilities, Standard_Deviation = sd_probabilities)

  return(df)

}



#' Computes the probability of a sequence beloning to a certain class of sRNA.
#'
#' @description This function computes the probability of a given sequence belonging
#'              to a certain class of sRNA.
#'
#' @usage sRNA_prob(rna_sequence, InvenireSRNA_model = NULL)
#'
#' @param rna_sequence RNA Sequence constructed with nucleotides A,C,G,T and U.
#' @param InvenireSRNA_model object of class InvenireSRNA to be used for prediction.
#'                           If \emph{InvenireSRNA_model = NULL} then the pretrained
#'                           model to predict CsrA regulating sRNAs is used.
#'
#' @return The functions returns a list containing the following elements:
#' \itemize{
#'   \item{\strong{probability} }{probability of the given sequence belonging to a given class of sRNA.}
#'   \item{\strong{std_deviation} }{the standard deviation of the prediction.}
#' }
#'
#' @export
sRNA_prob <- function(rna_sequence, InvenireSRNA_model = NULL){

  # compute features
  rna = as.character(rna_sequence)
  rna = tolower(rna)
  rna = gsub("t", "u", rna)
  RNAsubopt = 'RNAsubopt -p 1000'
  structs = system(RNAsubopt, intern = T,wait = TRUE, input = rna)
  features = get_features(rna, structs)
  features = c(1, features)

  if(is.null(InvenireSRNA_model)){
    ff = paste(system.file(package="InvenireSRNA"), "extdata/CsrA_models.tab", sep="/")
    coefficients_4_features = read.table(ff, header = T)
    InvenireSRNA_model = list(coefficients_4_features = coefficients_4_features)
    class(InvenireSRNA_model) = "InvenireSRNA"
  }

  if(class(InvenireSRNA_model) != "InvenireSRNA") stop("Please provide a model of class InvenireSRNA!")

  # compute probability
  prods = as.vector(apply(InvenireSRNA_model$coefficients_4_features, 1, function(x) sum(x * features)))
  probabilities = sapply(prods, function(x) 1/(1+exp(-x))) # apply the sigmoid for prediction
  probability = mean(probabilities)
  std_deviation = sd(probabilities)

  lst = list(probability = probability, std_deviation = std_deviation)

  return(lst)

}
