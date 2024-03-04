#' Abstract class containing time-dependent effects for MR analysis
#'
TimeDependentModel <- setClass(
  "TimeDependentModel",
  contains = "VIRTUAL",
  slots = c(
    model = "ANY"
  )
)

#' Return the underlying model of this instance
#'
#' @param this The object from which to extract the model from.
#'
#' @export
setGeneric("model", "TimeDependentModel", function(this) standardGeneric("model"))
setMethod("model", "TimeDependentModel", function(this){
  this@model
})

#' Calculate the total genetic effect at a certain timepoint
#'
#' @param this The model to use for calculation
#' @param age The ages (or time points) for which to calculate the genetic effect
#'
#' @export
setGeneric("totalEffect", "TimeDependentModel", function(this, age) standardGeneric("totalEffect"))
setMethod("totalEffect", "TimeDependentModel", function(this, age){
  rlang::abort("Not implemented for base class TimeDependentModel!")
})
