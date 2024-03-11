#' Models the relationship between age and time-dependent genetic effects using
#' a glm with interaction terms
#'
#' @slot model The underlting model ([`glm`][stats::glm()])
#' @slot fixedEffect The time-independent fixed genetic effect
#'   (`numeric` length 1)
#' @slot interactionEffects The coefficients of the interaction terms (`numeric`)
#' @slot exponents The exponent of age in the interaction terms (`numeric`)
#' @include TimeDependentModel-class.R
#' @docType class
#' @export
GlmTimeDependentModel <- setClass(
  "GlmTimeDependentModel",
  contains = "TimeDependentModel",
  slots = c(
    model = "glm",
    fixedEffect = "numeric",
    interactionEffects = "numeric",
    exponents = "numeric"
  )
)

setValidity("GlmTimeDependentModel", function(this){
  with(this, {
    if (length(fixedEffect) != 1)
      return(paste0("fixedEffect must have length 1. Is ", length(fixedEffect)))
    if (length(interactionEffects) != length(exponents))
      return(gettextf(
        "interactionEffects and exponents must be the same length. Are %d and %d",
        length(interactionEffects), length(exponents)
      ))
    TRUE
  })
})

#' Strength of the fixed effect
#'
#' @param this The model to get the fixed effect strength from
#' @rdname fixedEffect-methods
#' @docType methods
#' @export
setGeneric("fixedEffect", function(this) standardGeneric("fixedEffect"))
#' @rdname fixedEffect-methods
#' @aliases fixedEffect,GlmTimeDependentModel-method
setMethod("fixedEffect", "GlmTimeDependentModel", function(this) this@fixedEffect )

#' Effect strengths of the interaction terms. Same numer and order as `exponents`
#'
#' @param this The object to get the interaction effect strengths from
#' @docType methods
#' @rdname interactionEffects-methods
#' @export
setGeneric("interactionEffects", function(this) standardGeneric("interactionEffects"))
#' @rdname interactionEffects-methods
#' @aliases interactionEffects,GlmTimeDependentModel-method
setMethod("interactionEffects", "GlmTimeDependentModel", function(this) this@interactionEffects)

#' The exponents of age in the interaction terms. Same number and order as
#' `interactionEffects`
#'
#' @param this The object to get the exponents from
#' @docType methods
#' @rdname exponents-methods
#' @export
setGeneric("exponents", function(this) standardGeneric("exponents"))
#' @rdname exponents-methods
#' @aliases exponents,GlmTimeDependentModel-method
setMethod("exponents", "GlmTimeDependentModel", function(this) this@exponents)

#' @rdname totalEffect-methods
#' @aliases totalEffect,GlmTimeDependentModel-method
setMethod("totalEffect", "GlmTimeDependentModel", function(this, age){
  vapply(age, \(a) {
    fixedEffect(this) + sum(interactionEffects(this) * (a ^ exponents(this)))
  }, numeric(1))
})
