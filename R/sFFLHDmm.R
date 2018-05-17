#' sFFLHD maximin
#'
#' sFFLHD R6 object that gives a batch of points at a time using maximin.
#' To do this it takes all batches for stage at beginning of stage
#' and then reorders them.
#' Not that great in practice.
#' Requires extra optimization and storage.
#'
#' @field D numeric. The number of dimensions for the design. Must be set.
#' @field L numeric. The number of points in each batch, also the number of
#'  levels of each dimension. Must be set.
#' @field b integer. The batch number.
#' @field s sFFLHD. The design it takes the points and then reorders them.
#' @field X matrix. The points given in the design.
#' @field Xchoices list. Batches taken from s and have been reordered,
#' but which have not been returned to the user yet.
#'
#'
#'
#' @return A sFFLHDmm object
#'
#' @export
#'
#' @importFrom stats runif
#' @importFrom methods new
#' @importFrom conf.design factorize
#' @importFrom DoE.base oa.design
#'
#' @examples
#' s <- sFFLHDmm$new(D=2,L=3)
#' s$get.batch()
#' s <- sFFLHDmm$new(D=2,L=4)
#' s$get.batch()
sFFLHDmm <- R6::R6Class(classname="sFFLHDmm",
    public = list(
      s = NULL, # keep a sFFLHD
      b = 0,
      # have_choices = FALSE,
      Xchoices = list(),
      X = NULL,
      D = NULL,
      L = NULL,
      initialize = function(D, L, ...) {
        self$s <- sFFLHD::sFFLHD(D=D, L=L, ...)
        self$X <- matrix(NA, nrow=0, ncol=D)
        self$D <- D
        self$L <- L
      },
      get.batch = function() {#browser()
        if (self$b == 0) {
          bat <- self$s$get.batch()
        } else {
          if (length(self$Xchoices) == 0) {
            Xnew <- self$s$get.batches.to.golden()
            Xnewchoices <- split_matrix(mat=Xnew, rowspergroup = self$s$L,
                                        shuffle = FALSE)
            self$Xchoices <- Xnewchoices #c(self$Xchoices, )
            #self$s$get.batches.to.golden()
            #self$have_choices <- TRUE
          }
          #choices <- (nrow(s$Xb)) / s$L - self$b
          if (length(self$Xchoices) == 1) {
            selection <- 1
          } else {
          # now we have choices, do maximin
            mindists <- sapply(self$Xchoices,
                               function(xx) {self$mindist.1(self$X, xx)}
                               )
            selection <- which.max(mindists)
          }
          bat <- self$Xchoices[[selection]]
          self$Xchoices[[selection]] <- NULL
        }
        self$X <- rbind(self$X, bat)
        self$b <- self$b + 1
        bat
      },
      maximin = function() {
        for (i in 1:choices) {
          mindists[i] <- self$mindist.1(self$s$Xb[1:(self$b * self$s$L)],
                                        self$s$Xb)
        }
      },
      #mindist.all = function(Xlist, X1) {
      #  min(sapply(Xlist, function(aa){self$mindist.1(aa, X1)}))
      #}
      mindist.1 = function(a,b) {#browser()
        # a and b are matrices, points along row
        min(outer(1:nrow(a), 1:nrow(b),
                  Vectorize(function(i,j) {sum((a[i,] - b[j,])^2)})))
      },

      get.batches = function(num) { # get multiple batches at once
        out <- matrix(NA, nrow=0,ncol=self$D)
        for (i in 1:num) {out <- rbind(out,self$get.batch())}
        return(out)
      } # end get.batches function
    )
)
