#' Fisher Exact Test
#'
#' @param DataTab A 1-way frequency table as produced by \code{xtabs()}.
#' @param alpha A real number between 0 and 1.
#' @importFrom stats addmargins fisher.test
#'
#' @description
#'    The function \code{fisherExactTest()} takes a 2-way, 2xC frequency table and
#'    some value \eqn{\alpha} to conduct a hypothesis test of independence on the data
#'    A plot (inspired by `s20x::rowdistr`) is made to better understand the data. The
#'    Main function used is `fisher.test()` which provides the p-value of our test.
#'
#'
#' @return
#'    \code{fisherExactTest()} will return a named list of following values:
#'    \itemize{
#'    \item \code{pval} p value indicating the reslult of hypothesis test.
#'    \item \code{data} The passed in data with added margin sums.
#'    }
#'    A bar plot will also be generated that plots the row proportions of the data
#'    so we can better understand the data.
#'
#' @export
#'
fisherExactTest <- function(DataTab, alpha = 0.05){
  fish_out <- fisher.test(DataTab, alpha)

  #Details of data
  n <- sum(DataTab)
  length <- length(DataTab)
  num_row <- length(DataTab[,1])
  num_col <- length(DataTab[1,])
  n_row <- c()
  n_col <- c()
  p1 <- DataTab/n

  #Count row totals
  for(i in 1:num_row){
    n_row <- c(n_row, sum(DataTab[i,]))
  }

  for(i in 1:num_col){
    n_col <- c(n_col, sum(DataTab[,i]))
  }

  p <- c()

  #For loops to get label names and row props
  Label <- c() #Label names for plot
  for(i in 1:num_row){
    for(j in 1:num_col){
      p <- c(p, DataTab[i,j]/n_row[i])
    }
  }
  Label <- c(Label, rep(names(p1[,1]),num_col))


  ################# Plotting #############
  barplot(p, col = c("lightblue", "purple"),
          ylim = c(0, max(p) + 0.05), main = "Row Proportions")
  axis(1,at = 1:length,labels = Label)

  p <- matrix(p, num_row, num_col, byrow = TRUE)

  #add sum margin to data
  withmargin <- addmargins(DataTab, margin = c(1,2))
  print(withmargin)
  print(fish_out)

  #return list
  my_listf <- list(pval = fish_out$p.value, data = withmargin)

    return(invisible(my_listf))
}
