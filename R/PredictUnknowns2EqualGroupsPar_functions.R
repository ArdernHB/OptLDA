


#' Predict unknown specimens using LDA applied to PCA with resampling to equal sample size using parallel processing
#'
#' This function is designed for predicting unknown specimens to group by
#' carrying out an LDA on reduced dimensions by PCA. The function resamples
#' reference groups (aka training groups) to equal sample size and carries out
#' PCA on the reference groups. The unknown specimens (aka test specimens) are
#' then projected into the same PCA space. This ensures there is no data leakage
#' into the identification model from either the unequal reference group sample
#' sizes or from the unknown specimens. In other words, if the reference groups
#' are unequal then the larger groups will dispropotionately influence the PCA
#' upon which later predictive models will act on, this is known as data leakage
#' and is avoided using this function.
#'
#' The function can also be applied to shape data in which case it expects an
#' array of data with landmarks as rows, dimensions as columns and specimens as
#' sliced (for shape data set ShapeGPA=TRUE). In the case of shape data, both a GPA
#' and a PCA are carried out for each resampling excercise to equal sample sizes,
#' unknown specimens are then projected into the shape space defined by that equal
#' sample size and then those rotated specimens a projected into the PCA space.
#'
#' This function is for application to a set number of PCs once the optimal number
#' is established via correct cross-validation. To establish the optimal number
#' of PCs required see the function \code{LDACVManyStepwisePar} for stepwise
#' testing of each set of consecutive PCs.
#'
#' @param TrainingData matrix where specimens are rows. If shape data (ShapeGPA=TRUE) then the data is expected in an array format where rows are landmarks, columns are landmark dimensions and slices are specimens. The order of the specimens should match the GroupMembership vector.
#' @param UnknownData as for Training data but with the unknown data.
#' @inheritParams LDACVStepwisePar
#' @return Returns a matrix of the leave-one-out classifications for all the specimens along with their known classification.
#'
#'
#' @author Ardern Hulme-Beaman
#' @import MASS
#' @import Morpho
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import abind
#' @export


PredictUnknownsEqualPar <- function(TrainingData, UnknownData, GroupMembership, EqualIter=100, SampleSize=NA, ShapeGPA=FALSE, Sliding=NULL, SlidingLMindex=NULL, SizeShape=FALSE, PClim=10){
  #DiscriminationData=Scores
  #DiscriminationData=Shape
  #GroupMembership=Groups
  #EqualIter=100
  #PClim=3
  #ShapeGPA=TRUE; Verbose=FALSE; Sliding=NA;
  #SizeShape=FALSE; PCA=FALSE
  #KFold=4

  #TrainingData = Shape[,,-c(1:2)]#BlackRatGPA$rotated
  #UnknownData = Shape[,,1:2]#BlackRatGPA$rotated[,,1:4]
  #GroupMembership=Groups[-c(1:2)]
  #EqualIter=100
  #PClim=3; ShapeGPA = TRUE; KFold = 5


  #TrainingData = CompDatasort$ShapeData[,,ArchRem]; UnknownData = CompDatasort$ShapeData[,,-ArchRem]; GroupMembership= CompDatasort$info$species[ArchRem]
  #EqualIter=100; SampleSize=NA; ShapeGPA=TRUE; Sliding=NULL; SizeShape=FALSE; PClim=16

  if (ShapeGPA==TRUE & length(dim(TrainingData))==2){
    stop('TrainingData is a matrix, but ShapeGPA set to TRUE, if the data is shape data then it must be in array format, with LMs as rows, dimensions as columns and specimens as slices')
  }

  if (ShapeGPA==TRUE){
    if (dim(TrainingData)[3]!=length(GroupMembership)){
      stop('Number of specimens in TrainingData does not appear to match the number of speicmens listed in GroupMembership')
    }
  } else {
    if (dim(TrainingData)[1]!=length(GroupMembership)){
      stop('Number of specimens in TrainingData does not appear to match the number of speicmens listed in GroupMembership')
    }
  }



  chr2nu <- function(X){
    as.numeric(as.character(X))
  }

  ParOutput <- function(PreviousResults, ResultList){

    NewResults <- PreviousResults
    for (i in 1:length(ResultList)){
      if (is.vector(ResultList[[i]])){
        NewResults[[i]] <- cbind(PreviousResults[[i]], ResultList[[i]])
      } else if (length(dim(ResultList[[i]]))>=2){
        NewResults[[i]] <- abind::abind(PreviousResults[[i]], ResultList[[i]], along = 3)
      }

    }

    return(NewResults)
  }


  BalancedGrps <- function(GroupMembership, GroupSize=SampleSize){
    #DiscriminationData=PairwiseShapePCAmat
    #GroupMembership=chr(Groups[GrpPos])

    if (is.na(GroupSize)){
      minSampSize <- min(table(as.character(GroupMembership)))
    } else {
      minSampSize <- GroupSize
    }

    sampleindex <- 1:length(GroupMembership)
    RandSamp <- stats::aggregate(x = sampleindex, by=list(factor(GroupMembership)), sample, size=minSampSize)
    Index <- c(t(RandSamp$x))
    GroupMem <- c(t(stats::aggregate(RandSamp$Group.1, list(RandSamp$Group.1), rep, minSampSize)$x))
    return(list(IndexedLocations=Index, Newfactors=GroupMem))

  }

  Array2Mat <- function(Array){
    Matrix <- matrix(NA, nrow = dim(Array)[3], ncol = length(c(t(Array[,,1]))))
    for (i in 1:dim(Array)[3]){
      #i <- 1
      Matrix[i,] <- c(t(Array[,,i]))
    }
    return(Matrix)
  }

  ParEqualIterPredict <- function(TrainingData, UnknownData, GrpMem, ShapeGPA, Sliding, SlidingLMindex, SizeShape, PClim, SampleSize){
    #DiscriminationData=TrainingData; GrpMem=GroupMembership; ParTieBreaker='Report'; ParVerbose=FALSE
    #GrpMem=Groups; PClim=3; SampleSize=NA
    #DiscriminationData; GrpMem=GroupMembership; ShapeGPA=ShapeGPA; Sliding=Sliding;  PClim=PClim; SizeShape = SizeShape
    BalancingGrps <- BalancedGrps(GrpMem, SampleSize)

    #As we're cbinding the new factors with the folding, the folding will be
    #duplicated equally until it matches the length of the factors
    #these then ensures that the now balanced sample sizes are subjected to the
    #K-fold proceedure equally
    NewFactFoldIndex <- BalancingGrps$Newfactors

    if (ShapeGPA==TRUE){
      BalDataShape <- TrainingData[,,BalancingGrps$IndexedLocations]
      invisible(utils::capture.output(BalData <- Morpho::procSym(BalDataShape[,,], sizeshape = SizeShape, outlines = Sliding, SMvector = SlidingLMindex)))
      BalPCA <- stats::prcomp(Array2Mat(BalData$orpdata))

      BalTest <- Morpho::align2procSym(BalData, UnknownData)

      if (!is.null(Sliding)){
        BalTest <- UnknownData
        for (kspec in 1:dim(UnknownData)[3]){
          invisible(utils::capture.output(BalTest[,,kspec] <- Morpho::relaxLM(reference=BalData$mshape, lm=UnknownData[,,kspec], outlines = Sliding, SMvector = SlidingLMindex)))
        }

      } else {
        BalTest <- suppressMessages(Morpho::align2procSym(BalData, UnknownData))

      }



      BalTestPCA <- stats::predict(BalPCA, newdata = Array2Mat(BalTest))
    } else {
      BalData <- TrainingData[BalancingGrps$IndexedLocations,]
      BalPCA <- stats::prcomp(x = BalData)

      BalTestPCA <- stats::predict(BalPCA, newdata = UnknownData)
    }

    LDAres <- MASS::lda(x = BalPCA$x[,1:PClim], grouping=BalancingGrps$Newfactors, CV=FALSE)
    LDApredict <- stats::predict(LDAres, newdata = BalTestPCA[,1:PClim])


    return(list(UnknownCalls=as.character(LDApredict$class), GrpPredictProbs=LDApredict$posterior))


  }


  cores <- parallel::detectCores()
  clust <- parallel::makeCluster(cores[1]-1)
  doParallel::registerDoParallel(clust)

  a <- 1
  ParResults <- foreach::foreach(a = 1:EqualIter, .combine = ParOutput) %dopar% {
    ParEqualIterPredict(TrainingData, UnknownData, GrpMem=GroupMembership, ShapeGPA=ShapeGPA, Sliding=Sliding, SlidingLMindex=SlidingLMindex, PClim=PClim, SizeShape = SizeShape, SampleSize = SampleSize)
  }

  parallel::stopCluster(clust)

  return(ParResults)
}

