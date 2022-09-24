#' Oneway Table Analysis
#'
#' @param DataTab A 1-way frequency table as produced by \code{xtabs()}.
#' @param alpha A real number between 0 and 1.
#' @importFrom grDevices col2rgb rgb
#' @importFrom graphics abline arrows barplot axis
#' @importFrom stats qnorm
#'
#' @description
#'    The function \code{oneWay()} takes a 1-way frequency table and
#'    some value \eqn{\alpha} to create a \eqn{(1-\alpha)100\%} confidence
#'    interval on each estimate \eqn{p} from the table. A plot (inspired by
#'    `s20x`) is made to show each \eqn{p} value along with their respective
#'    confidence intervals and the null \eqn{p_0} value. A similar process will
#'    be run on the differences of each proportion
#'
#'
#' @details
#'    The values for \eqn{p_i} are calculated as
#'    \eqn{\frac{\# \text{ in group } i}{n}}. Each confidence interval is
#'    made from a lower and upper bound calculated from Wald's method:
#'
#'    \eqn{p \pm z_{\alpha/2}\sqrt{\frac{p(1-p)}{n}}}
#'
#'
#' @return
#'    \code{oneWay()} will return a named list of following values:
#'    \itemize{
#'    \item \code{prob} the estimation for \eqn{p}.
#'    \item \code{Conf} confidence interval data for each \eqn{p}.
#'    \item \code{Diff} differences between each \eqn{p}.
#'    }
#'    A bar plot will also be generated that plots the values of \eqn{p} along with
#'    their respective confidence interval along with a horizontal line at
#'    \code{h = 1/length{DataTab}}. Each bar is shaded according to its \eqn{p} value;
#'    the larger \eqn{p} value corresponds to a darker bar.
#'
#'    A scatter plot and intervals are also returned.
#'
#' @export
#'
#' @examples
#' ## This is a contingency table in array form
#' DF <- as.data.frame(UCBAdmissions)
#'
#' ## Make Cross Table
#' tab <- xtabs(Freq ~ Dept, DF)
#'
#' ## Run the Function
#' oneWay(tab, alpha = 0.05)
oneWay <- function(DataTab, alpha = 0.05){
  ################# Calculations #################

  out <- chisq.test(DataTab, p = rep(1/length(DataTab),length(DataTab)))
  #This Functions creates a confidence interval based on Wald's Method
  #Two tailed confidence interval
  alph2 <- alpha/2
  z = qnorm(1-alph2)

  Data.f = data.frame(DataTab)
  TabSum = sum(DataTab)

  p = DataTab/TabSum #Vector of probabilities
  sd = sqrt(p *(1-p)/TabSum) #standard deviation

  #Lower and Upper bound of CI
  Lower = round(p-z*sd, 4)
  Upper = round(p+z*sd, 4)

  CIs <- data.frame(cbind(p = round(p, 4), Lower, Upper)) #bind by column

  ################# Calculations for diff #################

  #Two tailed Confidence Interval
  alph2 <- alpha/2
  z = qnorm(1-alph2)

  length = length(DataTab)
  n = sum(DataTab)

  #Vectors that will hold the upper and lower bounds of the confidence intervals
  Lower2 <- c()
  Upper2 <- c()
  Diff <- c()
  L <- c()

  #Double for Loop to calculate each difference
  #Not all of these are needed, but this was the easiest way I could think
  #to do the calculations.
  for (i in 1:length){
    for (j in 1:length){
      #Differences
      Diff[length*(i-1) + j] <- p[i] - p[j]

      #Calculate Lower Bound
      Lower2 <- c(Lower2, p[i] - p[j] - z*sqrt(
        (p[i]*(1-p[i]))/n + (p[j]*(1-p[j]))/n + (2*p[i]*p[j])/n))

      #Calculate Upper Bound
      Upper2 <- c(Upper2, p[i] - p[j] + z*sqrt(
        (p[i]*(1-p[i]))/n + (p[j]*(1-p[j]))/n + (2*p[i]*p[j])/n))

      #Labels
      L <- c(L, paste0("p[", i, "] - p[", j, "]"))
    }
  }

  ################# Refine Data ##################
  #This will refine the data to the top right triangle of a matrix that contians
  # the difference, lower bound, and upper bound
  RDif <- c()
  RLow <- c()
  RUp  <- c()
  Labels <- c()
  shift = 2
  for(i in 0:(length-2)){
    for(j in shift:length){
      RDif <- c(RDif, Diff[length*i + j])
      RLow <- c(RLow, Lower2[length*i + j])
      RUp <- c(RUp, Upper2[length*i + j])
      Labels <- c(Labels, L[length*i + j])
    }
    shift = shift + 1
  }

  ################# Plotting #####################

  par(mfrow = c(2,1), mar = c(4,2,3,2))
  # `ColPcent` helps make colors more transparent based on `p`.
  # The code is modified from:
  # https://www.dataanalytics.org.uk/make-transparent-colors-in-r

  ColPcent <- function(x = 0.5, cols = c('black') ){
    val <- col2rgb(cols)
    newcols <- rgb(val[1,],val[2,],val[3,], maxColorValue = 255, alpha = round((255*x)))
    invisible(newcols)
  }

  #Mids is the points where to place CIs
  mids <- barplot(CIs$p ~ Data.f[,1], col = ColPcent(CIs$p),
                  ylim = c(0, max(Upper)+.05),
                  main = "Proportions by Catagory \nwith Confidence Interval",
                  xlab = (colnames(Data.f))[1], ylab = "Porportion")

  #Add line of null p value
  abline(h = 1/length(DataTab), col = "blue4", lwd = 1.75, lty = 5)

  #plot the confidence intervals on bar plot
  arrows(mids, Lower, mids, Upper, length = 0.1, code = 3, angle = 90,
         col = "firebrick1", lwd = 2)


  plot(1:length(RDif), RDif, xaxt = 'n' ,col = "black", pch = 0, lwd = 3,
       ylim = c(min(RLow), max(RUp)), xlab = "Difference")
  axis(1,at = 1:length(RDif),labels = Labels)
  abline(h = 0, col = "blue4", lwd = 1.75, lty = 5)
  arrows(1:length(RDif), RLow, 1:length(RDif), RUp, length = 0.1, code = 3, angle = 90,
         col = "firebrick1", lwd = 2)

  my_list <- list(prob = as.vector(p), Conf = CIs, Diff = RDif, pval = out$p.value)

  print(my_list)

  return(invisible(my_list))
}
