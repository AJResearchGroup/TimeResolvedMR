#' Abstract class containing time-dependent effects for MR analysis
#'
#' @docType class
#' @slot model The underlying model
TimeDependentModel <- setClass(
  "TimeDependentModel",
  contains = "VIRTUAL",
  slots = c(
    model = "ANY"
  )
)

#' Return the underlying model of this instance
#'
#' @param this The object to get the model from.
#' @docType methods
#' @rdname model-methods
#' @export
setGeneric("model", function(this) {standardGeneric("model")})
#' @rdname model-methods
#' @aliases model,TimeDependentModel-method
setMethod("model", "TimeDependentModel", function(this){
  this@model
})

#' Calculate the total genetic effect at a certain timepoint
#'
#' @param this The model to use for calculation
#' @param age The ages (or time points) for which to calculate the genetic effect
#'
#' @export
#' @docType methods
#' @rdname totalEffect-methods
setGeneric("totalEffect", function(this, age) standardGeneric("totalEffect"))
#' @rdname totalEffect-methods
#' @aliases totalEffect,TimeDependentModel-method
setMethod("totalEffect", "TimeDependentModel", function(this, age){
  rlang::abort("Not implemented for base class TimeDependentModel!")
})
