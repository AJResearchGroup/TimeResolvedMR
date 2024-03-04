#' Models the relationship between age and time-dependent genetic effects using
#' a glm with interaction terms
#'
#' @slot model The underlting model ([`glm`][stats::glm()])
#' @slot fixedEffect The time-independent fixed genetic effect
#'   (`numeric` length 1)
#' @slot interactionEffects The coefficients of the interaction terms (`numeric`)
#' @slot exponents The exponent of age in the interaction terms (`numeric`)
#' @include TimeDependentModel-class.R
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
#' @keywords internal
#' @rdname GlmTimeDependentModel-methods
#' @export
setGeneric("fixedEffect", function(this) standardGeneric("fixedEffect"))
setMethod("fixedEffect", "GlmTimeDependentModel", function(this) this@fixedEffect )

#' Effect strength of the interaction terms. Same numer and order as `exponents`
#' @keywords internal
#' @rdname GlmTimeDependentModel-methods
#' @export
setGeneric("interactionEffects", function(this) standardGeneric("interactionEffects"))
setMethod("interactionEffects", "GlmTimeDependentModel", function(this) this@interactionEffects)

#' The exponents of age in the interaction terms. Same number and order as
#' `interactionEffects`
#' @keywords internal
#' @rdname GlmTimeDependentModel-methods
#' @export
setGeneric("exponents", function(this) standardGeneric("exponents"))
setMethod("exponents", "GlmTimeDependentModel", function(this) this@exponents)

setMethod("totalEffect", "GlmTimeDependentModel", function(this, age){
  vapply(age, \(a) {
    fixedEffect(this) + interactionEffects(this) * (a ^ exponents(this))
  }, numeric(1))
})
