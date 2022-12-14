#' Two Way Table Analysis
#'
#' @param DataTable A 1-way frequency table as produced by \code{xtabs()}.
#' @param alpha A real number between 0 and 1.
#' @importFrom graphics par points
#' @importFrom stats chisq.test
#'
#' @description
#'    The function \code{twoWay()} takes a 2-way frequency table and
#'    some value \eqn{\alpha} to study. A plot (inspired by
#'    `s20x::rowdistr`) is created to help compare the data and better
#'    understand if there is independece between the data.
#'
#' @return
#'    \code{twoWay()} will return a named list of following values:
#'    \itemize{
#'    \item \code{Props} the estimation for \eqn{p} on each row.
#'    \item \code{Counts} counts of each occurance.
#'    \item \code{Expected} expected count of each data, generated by \code{chisq.test}.
#'    \item \code{RowProps} \eqn{p} by row
#'    \item \code{pval} p-value to accept or reject null.
#'    }
#'    A bar plot will also be generated that plots the row proportions along with
#'    their respective confidence interval. Each bar is shaded according to
#'    its row.
#' @export
#'
#'
twoWay <- function(DataTable, alpha = 0.05){

  alph2 <- alpha/2
  z = qnorm(1-alph2)

  n <- sum(DataTable)
  length <- length(DataTable)
  length1 <- length(DataTable[1,]) #number of cols
  length2 <- length(DataTable[,1]) #number of rows

  Chi_out <- chisq.test(DataTable)
  p1 <- DataTable/n

  p <- c() #row probabilities
  nr <- c() #total in each row

  #Find total in each row
  for(i in 1:length1){
    nr <- c(nr, sum(DataTable[i,]))
  }

  #For loops to get label names and row props
  Label <- c() #Label names for plot
  for(i in 1:length1){
    for(j in 1:length2){
      p <- c(p, DataTable[j,i]/nr[j])
    }
    Label <- c(Label, names(p1[,1]))
  }

  #Lower and Upper bound of CI of rows
  sd = sqrt(p *(1-p)/n) #standard deviation
  Lower = (p-z*sd)
  Upper = (p+z*sd)

  #Row confidence interval
  CIs <- data.frame(cbind(p = round(p, 4), Lower, Upper))


  ################# Plotting #####################

  mid <- barplot(p, col = c("lightblue", "purple"),
                 ylim = c(0, max(p) + 0.05), main = "Row Proportions")

  axis(1,at = mid,labels = Label)
  abline(v = mean(mid))
  points(mid,Chi_out$expected/n, pch = 4, lwd = 3)
  arrows(mid, Lower, mid, Upper, length = 0.1, code = 3, angle = 90,
         col = "firebrick1", lwd = 2)

  #named list
  my_list <- list(Props = DataTable/n, Counts = Chi_out$observed,
                  data = addmargins(DataTable, margin = c(1,2)),
                  Expected = Chi_out$expected, RowProps = p, CI = CIs,
                  pval = Chi_out$p.value)

  print(addmargins(DataTable, margin = c(1,2)))
  cat("Row Props: \n", p, "\n")
  print(Chi_out)
  return(invisible(my_list))
}
