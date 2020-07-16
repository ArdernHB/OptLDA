

#' Plotting function for stepwise results
#'
#' This function plots the results of a stepwise analyses with resampling to equal sample size.
#' The function will plot the line connecting the means of each stepwise increment increase of PC
#' and will also plot the range around the mean. The default upper and lower limit is set to the 5th
#' and 95th percentiles.
#'
#' @param Percentiles a vector of two values denoting the upper and lower limit of the resampling to be plotted. Values must be between 0 and 1. Default is .05 and .95.
#' @param StepwiseResultsMat a matrix where the rows are the cross-validation percentage result for each iteration of resampling to equal sample size and the columns are the stepwise increase in PC.
#' @param PlotCol is the colour of the line and the polygon representing the range of cross-validation percentage at the user defined percentiles.
#' @param Add is a logical value to determine if the the plot should be added to a previous plot. For example different grouping approaches could be used and the results compared by adding the plots to the same graph.
#' @param Xlabel is the label to be used on the x axis. Default is 'PC'.
#' @return Plots a graph of mean CCV percentages with the percentiles plotted as a polygon range around the mean.
#'
#'
#' @keywords plotting
#' @author Ardern Hulme-Beaman
#'
#'
#' @export

PlotStepwise <- function(StepwiseResultsMat, Percentiles=c(.05, .95), PlotCol='darkblue', Add=FALSE, Xlabel='PC'){

  #StepwiseResultsMat = VoleStepLDA; PlotCol = 'darkblue'; Add = FALSE

  DataDim <- dim(StepwiseResultsMat)[2]

  if (Add==FALSE){
    graphics::plot(y = colMeans(StepwiseResultsMat), x = 1:DataDim, type = 'n', ylim = c(10,105), ylab = 'CCV %', xlab = Xlabel)
    graphics::abline(h = seq(from = 20, to = 100, by = 10), v = seq(from = 2, to = DataDim, by =2), lty = '1919')

  }


  ResRange <- apply(StepwiseResultsMat, MARGIN = 2, FUN = stats::quantile, probs = Percentiles, na.rm=TRUE)
  graphics::polygon(x = c(1:DataDim, DataDim:1),
                    y = c(ResRange[1,], ResRange[2,DataDim:1])*100,
                    col = transpar(PlotCol, alpha = 75),
                    border = NA)

  graphics::lines(y = colMeans(StepwiseResultsMat*100), x = 1:DataDim, col=PlotCol, lwd=3)

}




#' Internal function: Transparent named colour
#'
#' This function takes a named colour and returns the transparent equivalent
#' @param Colour A colour name from colours() function which is desired in transparent form.
#' @param alpha The level of transparency from 1 (completely transparent) to 100 (completely opaque) that the returned colour should be.
#' @return The transparent equivalent of a named colour
#' @keywords internal
#' @keywords colour
#' @keywords transparency
#' @author Ardern Hulme-Beaman


transpar<-function(Colour, alpha=100){
  newColour<-grDevices::col2rgb(Colour)
  apply(newColour, 2, function(curcoldata){grDevices::rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}





