#' @include internal-asserts.R
#' @importFrom Rdpack reprompt
NULL

# -----------------------------------------------------------------------------
# Conjugate NIW / NIX Parameter Conversions
# -----------------------------------------------------------------------------

#' Convert Between NIW/NIX Hyperparameters and Expected/Marginal Moments
#'
#' Functions to convert between conjugate Normal-Inverse-Wishart (NIW) or
#' Normal-Inverse-Chi-Squared (NIX) hyperparameters (\eqn{m, S}) and the
#' expected or marginal distribution moments (\eqn{\mu, \Sigma})
#' \insertCite{@see @murphy2012 p. 134}{MVBeliefUpdatr}.
#'
#' Under the NIW model with \eqn{D} dimensions, scatter matrix \eqn{S}, and
#' degrees of freedom / pseudocounts \eqn{\nu > D + 1}, the expected
#' covariance matrix is:
#' \deqn{\Sigma = \frac{S}{\nu - D - 1}}
#' Conversely, the scatter matrix is obtained from \eqn{\Sigma} by:
#' \deqn{S = \Sigma \cdot (\nu - D - 1)}.
#' For conjugate NIW and NIX models, the expected mean \eqn{\mu} is identically
#' equal to the location parameter \eqn{m}:
#' \deqn{\mu = m}
#'
#' @param S Scatter matrix \eqn{S} (numeric matrix, scalar, or list of
#'   matrices).
#' @param Sigma Expected category covariance matrix \eqn{\Sigma} (numeric
#'   matrix, scalar, or list of matrices).
#' @param nu Strength of belief (pseudocount) about \eqn{\Sigma}
#'   (\eqn{\nu > D + 1}).
#' @param kappa Strength of belief (pseudocount) about \eqn{\mu}.
#' @param m Location hyperparameter \eqn{m} (numeric vector, matrix, or list).
#' @param mu Expected category mean vector \eqn{\mu} (numeric vector, matrix,
#'   or list).
#' @return For \code{get_expected_Sigma_from_S}, the expected covariance
#'   matrix \eqn{\Sigma = S / (\nu - D - 1)}. For
#'   \code{get_S_from_expected_Sigma}, the scatter matrix
#'   \eqn{S = \Sigma \cdot (\nu - D - 1)}. For
#'   \code{get_marginal_Sigma_from_S}, the marginal category covariance matrix
#'   \eqn{\Sigma = (\kappa + 1) / \kappa \cdot S / (\nu - D - 1)}. For
#'   \code{get_expected_mu_from_m} and \code{get_m_from_expected_mu}, the
#'   mean vector or location hyperparameter.
#'
#' @references \insertRef{murphy2012}{MVBeliefUpdatr}
#' @seealso \code{\link{get_expected_sigma}}, \code{\link{get_expected_mu}}
#' @rdname get_expected_Sigma_from_S
#' @export
#' @examples
#' # 2D example
#' S <- matrix(c(4, 1, 1, 3), nrow = 2)
#' nu <- 10
#' kappa <- 5
#'
#' # Convert scatter matrix to expected covariance matrix (D = 2)
#' Sigma <- get_expected_Sigma_from_S(S, nu = nu)
#' Sigma
#'
#' # Convert expected covariance matrix back to scatter matrix
#' get_S_from_expected_Sigma(Sigma, nu = nu)
#'
#' # Marginal covariance matrix
#' get_marginal_Sigma_from_S(S, nu = nu, kappa = kappa)
#'
#' # Mean and location conversions
#' m <- c(0, 1)
#' mu <- get_expected_mu_from_m(m)
#' get_m_from_expected_mu(mu)
get_expected_Sigma_from_S <- function(S, nu) {
  is_single <- !is.list(S)
  if (is_single) {
    S <- list(S)
    nu <- list(nu)
  }

  Sigma <- mapply(
    function(S_i, nu_i) {
      if (is.null(S_i) || is.null(nu_i)) {
        return(NULL)
      }
      D <- if (is.matrix(S_i)) {
        ncol(S_i)
      } else if (is.list(S_i)) {
        ncol(as.matrix(S_i[[1L]]))
      } else {
        length(S_i)
      }
      if (is.na(nu_i) || nu_i <= D + 1) {
        if (is.matrix(S_i)) {
          matrix(NA_real_, nrow = nrow(S_i), ncol = ncol(S_i))
        } else {
          NA_real_
        }
      } else {
        S_i / (nu_i - D - 1)
      }
    },
    S, nu,
    SIMPLIFY = FALSE
  )

  if (all(vapply(Sigma, function(x) length(x) == 1L, logical(1L)))) {
    Sigma <- unlist(Sigma)
  }
  if (is_single) Sigma[[1L]] else Sigma
}

#' @rdname get_expected_Sigma_from_S
#' @export
get_S_from_expected_Sigma <- function(Sigma, nu) {
  is_single <- !is.list(Sigma)
  if (is_single) {
    Sigma <- list(Sigma)
    nu <- list(nu)
  }

  S <- mapply(
    function(Sigma_i, nu_i) {
      if (is.null(Sigma_i) || is.null(nu_i)) {
        return(NULL)
      }
      D <- if (is.matrix(Sigma_i)) {
        ncol(Sigma_i)
      } else if (is.list(Sigma_i)) {
        ncol(as.matrix(Sigma_i[[1L]]))
      } else {
        length(Sigma_i)
      }
      Sigma_i * (nu_i - D - 1)
    },
    Sigma, nu,
    SIMPLIFY = FALSE
  )

  if (all(vapply(S, function(x) length(x) == 1L, logical(1L)))) {
    S <- unlist(S)
  }
  if (is_single) S[[1L]] else S
}

#' @rdname get_expected_Sigma_from_S
#' @export
get_marginal_Sigma_from_S <- function(S, nu, kappa) {
  is_single <- (!is.list(S) && length(nu) == 1L && length(kappa) == 1L) ||
    (is.list(S) && length(S) == 1L && length(nu) == 1L && length(kappa) == 1L)

  Sigma_exp <- get_expected_Sigma_from_S(S, nu)
  if (!is.list(Sigma_exp)) {
    Sigma_exp <- as.list(Sigma_exp)
  }
  if (!is.list(kappa)) {
    kappa <- as.list(kappa)
  }

  Sigma_marg <- mapply(
    function(sig_i, kap_i) {
      if (is.null(sig_i) || is.null(kap_i) ||
          any(is.na(kap_i)) || any(kap_i <= 0)) {
        return(sig_i)
      }
      ((as.numeric(kap_i)[1L] + 1) / as.numeric(kap_i)[1L]) * sig_i
    },
    Sigma_exp, kappa,
    SIMPLIFY = FALSE
  )
  if (is_single) {
    Sigma_marg[[1L]]
  } else {
    Sigma_marg
  }
}

#' @rdname get_expected_Sigma_from_S
#' @export
get_expected_mu_from_m <- function(m) {
  if (!is.list(m)) {
    return(m)
  }
  if (all(vapply(m, function(x) length(x) == 1L, logical(1L)))) {
    mu <- unlist(m)
  } else {
    mu <- m
  }
  if (length(mu) == 1L) mu <- mu[[1L]]
  mu
}

#' @rdname get_expected_Sigma_from_S
#' @export
get_m_from_expected_mu <- function(mu) {
  if (!is.list(mu)) {
    return(mu)
  }
  if (all(vapply(mu, function(x) length(x) == 1L, logical(1L)))) {
    m <- unlist(mu)
  } else {
    m <- mu
  }
  if (length(m) == 1L) m <- m[[1L]]
  m
}

