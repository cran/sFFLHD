#' sFFLHD with flexible L
#'
#' R6 object that gives uses a sFFLHD with L near the requested one,
#' but gives them back in the requested L
#'
#' @field D numeric. The number of dimensions for the design. Must be set.
#' @field L numeric. The number of points in each batch, also the number of
#'  levels of each dimension. Must be set.
#' @field b integer. The batch number.
#' @field s sFFLHD. The design it takes the points and then reorders them.
#' @field X matrix. The points given in the design.
#' @field X_choices matrix. Points taken from s and have been reordered,
#' but which have not been returned to the user yet.
#'
#'
#'
#' @return A sFFLHD_Lflex object
#'
#' @export
#'
#' @importFrom DoE.base oa.design
#'
#' @examples
#' s <- sFFLHD_Lflex$new(D=8,L=4)
#' s$get.batch()
#' # sFFLHD(D=7,L=10)$get.batch() doesn't work, needs L=7,8,9,11
#' s <- sFFLHD_Lflex$new(D=7,L=10) # Uses L=9
#' s$get.batch()
#' s <- sFFLHD_Lflex$new(D=7,L=10, prefer_L="up") # Should use 11
sFFLHD_Lflex <- R6::R6Class(
  classname="sFFLHD_Lflex",
  public = list(
    s = NULL, # sFFLHD object
    D = NULL, # Dimensions
    L = NULL, # Requested L
    L_used = NULL, # L actually used
    prefer_L = NULL, # "near", "down", "up"
    Xchoices = NULL,
    initialize = function(D, L, prefer_L = "near", ...) {
      self$D <- D
      self$L <- L
      self$Xchoices <- matrix(0, ncol=D, nrow=0)
      self$prefer_L <- prefer_L

      # Check for requested OA
      OA.avail <- DoE.base::show.oas(nruns=L^2, factors=list(nlevels=L, number=D+1), show=0)
      # If avail, then use that L
      if (!is.null(OA.avail)) {
        # OA <- DoE.base::oa.design(nruns=L^2, nfactors=D+1, nlevels=L, columns="min3")
        self$L_used = L
      }

      # If it wasn't available, then check other L to tell user what to try instead
      # if (inherits(OA, "try-error")) {
      else {
        #avail.oas <- DoE.base::show.oas(factors=list(nlevels=L, number=D+1))
        # Check L in 2 to 16
        avail.oas <- sapply(2:16,
                            function(i) {
                              capture.output(av <- DoE.base::show.oas(nruns=i^2, factors=list(nlevels=i, number=D+1), show=0))
                              !is.null(av)
                            }
        )

        if (all(!avail.oas)) {
          stop("No OA can be found for L in 2 to 16 for given D")
        } else {
          avail.Ls <- (2:16)[avail.oas]

          # stop(paste("Try L one of",paste(avail.Ls, collapse=', '),"instead"))

          # Select L according to prefer_L
          if (prefer_L == "down") {
            # Find available below, pick highest
            low.avail <- avail.Ls[avail.Ls<L]
            if (length(low.avail) > 0) {self$L_used <- low.avail[length(low.avail)]}
            # Otherwise pick smallest above
            else {self$L_used <- min(avail.Ls)}
          } else if (prefer_L == "up") {
            # Find available above, pick 1st
            high.avail <- avail.Ls[avail.Ls>L]
            if (length(high.avail) > 0) {self$L_used <- high.avail[1]}
            # Otherwise pick highest below
            else {self$L_used <- max(avail.Ls)}

          } else if (prefer_L == "near") {
            # Pick value nearest is abs val
            self$L_used <- avail.Ls[which.min(abs(avail.Ls-L))]
          } else {
            stop('prefer_L must be one of "near", "down", "up"')
          }

        }
        # print(paste("picking", self$L_used))
        message("sFFLHD_Lflex requested L=", self$L,", but using L=", self$L_used)
        # self$L_used <- L_used
        # OA <- DoE.base::oa.design(nruns=nruns, nfactors=D+1, nlevels=L, columns="min3")

      }
      # Create s using L_used
      self$s <- sFFLHD::sFFLHD(D=self$D, L=self$L_used, ...)
    },
    get.batch = function() {#browser()
      # If not enough points in Xchoices
      while (nrow(self$Xchoices) < self$L) {
        # Add points to Xchoices
        newbatch <- self$s$get.batch()
        self$Xchoices <- rbind(self$Xchoices, newbatch)
      }
      returnbatch <- self$Xchoices[1:self$L, , drop=FALSE]
      self$Xchoices <- self$Xchoices[-(1:self$L), , drop=FALSE]
      returnbatch
    },

    get.batches = function(num) { # get multiple batches at once
      out <- matrix(NA, nrow=0,ncol=self$D)
      for (i in 1:num) {out <- rbind(out,self$get.batch())}
      return(out)
    } # end get.batches function
  )
)
