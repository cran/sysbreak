## Tests for robust = "hc": the heteroskedasticity-consistent estimate is checked
## against direct computations; "hac" is checked against values of the 1.1.0 build.

stack_sys <- function(y, z) {
  n <- ncol(y)
  list(maty = matrix(as.vector(t(y)), ncol = 1),
       matz = do.call(rbind, lapply(seq_len(nrow(y)), function(t)
         kronecker(diag(n), z[t, , drop = FALSE]))),
       n = n, bigt = nrow(y))
}

sim_het1 <- function(bigt = 80, br = 40) {
  set.seed(21)
  x <- rnorm(bigt)
  u <- rnorm(bigt) * sqrt(0.2 + 0.8 * x^2)
  y <- 1 + 0.5 * x + ifelse(seq_len(bigt) > br, 0.8 * x, 0) + u
  list(y = matrix(y), z = cbind(1, x))
}

sim_het2 <- function(bigt = 120, br = 60) {
  set.seed(22)
  x <- rnorm(bigt)
  e <- matrix(rnorm(2 * bigt), bigt) * sqrt(0.2 + 0.8 * x^2)
  y1 <- 1 + 0.5 * x + ifelse(seq_len(bigt) > br, 1.2, 0) + e[, 1]
  y2 <- -0.5 + 0.3 * x + ifelse(seq_len(bigt) > br, -0.8, 0) + 0.5 * e[, 1] + e[, 2]
  data.frame(y1 = y1, y2 = y2, x = x)
}

hc0_wald <- function(y, z, br) {
  bigt <- nrow(z); q <- ncol(z)
  d1 <- seq_len(bigt) <= br
  X <- cbind(z * d1, z * (!d1))
  b <- solve(crossprod(X), crossprod(X, y))
  u <- as.vector(y - X %*% b)
  xi <- solve(crossprod(X))
  V <- xi %*% crossprod(X * u) %*% xi
  Rm <- cbind(-diag(q), diag(q))
  list(b = b, u = u, X = X, V = V,
       wald = as.numeric(t(Rm %*% b) %*% solve(Rm %*% V %*% t(Rm)) %*% (Rm %*% b)))
}

collect_warnings <- function(expr) {
  w <- character()
  val <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}

test_that("correct() with robust = \"hc\" is the outer-product estimate", {
  set.seed(1)
  v <- matrix(rnorm(300), 100, 3)
  hc <- sysbreak:::correct(v, prewhit = FALSE, robust = "hc")
  expect_equal(hc, crossprod(v) / (nrow(v) - ncol(v)))
  expect_identical(sysbreak:::correct(v, prewhit = TRUE, robust = "hc"), hc)
})

test_that("coefse() and supft() with robust = \"hc\" equal the White formulas (n = 1)", {
  d <- sim_het1()
  s <- stack_sys(d$y, d$z)
  ref <- hc0_wald(d$y, d$z, 40)
  bigt <- s$bigt; k <- 2 * ncol(d$z)
  vv <- rbind(mean(ref$u^2), mean(ref$u^2))
  se <- sysbreak:::coefse(s$maty, s$matz, n = 1, m = 1, bigt = bigt, br = 40,
                          beta = ref$b, vv = vv, R = NULL, robust = "hc",
                          brv = FALSE, prewhit = FALSE)
  expect_equal(unname(se$vcov), unname(bigt / (bigt - k) * ref$V), tolerance = 1e-10)
  f <- sysbreak:::supft(s$maty, s$matz, 1, 1, 40, ref$b, bigt,
                        prewhit = FALSE, robust = "hc")
  # the T - d divisor of correct() and the (T - (m + 1) q) / T factor of supft()
  expect_equal(f, ((bigt - k) / bigt)^2 * ref$wald, tolerance = 1e-10)
})

test_that("seqf() with robust = \"hc\" is the maximum of the White statistics", {
  d <- sim_het1()
  s <- stack_sys(d$y, d$z)
  h <- 16
  got <- sysbreak:::seqf(s$maty, s$matz, 1, integer(0), s$bigt, h,
                         prewhit = FALSE, robust = "hc", trm = h / s$bigt)
  k <- 2 * ncol(d$z)
  want <- ((s$bigt - k) / s$bigt)^2 *
    max(vapply(h:(s$bigt - h), function(tau) hc0_wald(d$y, d$z, tau)$wald, numeric(1)))
  expect_equal(got, want, tolerance = 1e-6)
})

test_that("intpc() with robust = \"hc\" equals the Theorem 3 moments", {
  dat <- sim_het2(bigt = 200, br = 100)
  y <- as.matrix(dat[, c("y1", "y2")]); z <- cbind(1, dat$x)
  s <- stack_sys(y, z); n <- 2; q <- ncol(s$matz); bigt <- s$bigt; br <- 100
  rsq <- function(S) { e <- eigen(S, symmetric = TRUE); e$vectors %*% diag(sqrt(e$values)) %*% t(e$vectors) }
  for (cfg in list(c(FALSE, FALSE), c(TRUE, FALSE), c(FALSE, TRUE), c(TRUE, TRUE))) {
    hetq <- cfg[1]; brv <- cfg[2]
    fit <- sysbreak:::r_estim(s$maty, s$matz, bigt, n, 1, br, NULL,
                              brbeta = TRUE, brv = brv)
    beta <- fit$nbeta; vv <- fit$nvv
    S1 <- vv[1:n, ]; S2 <- vv[(n + 1):(2 * n), ]
    u <- matrix(as.vector(s$maty - sysbreak:::pzbar(s$matz, 1, br, bigt) %*% beta),
                ncol = n, byrow = TRUE)
    res <- u
    res[1:br, ] <- u[1:br, ] %*% solve(rsq(S1))
    res[(br + 1):bigt, ] <- u[(br + 1):bigt, ] %*% solve(rsq(S2))
    A1 <- rsq(S1) %*% solve(S2) %*% (S2 - S1) %*% solve(rsq(S1))
    A2 <- rsq(S2) %*% solve(S1) %*% (S2 - S1) %*% solve(rsq(S2))
    db <- beta[(q + 1):(2 * q), , drop = FALSE] - beta[1:q, , drop = FALSE]
    side <- function(rows, Sown, Soth) {
      Pi <- matrix(0, q, q); Q <- matrix(0, q, q); Om <- matrix(0, n^2, n^2)
      for (t in rows) {
        X <- s$matz[((t - 1) * n + 1):(t * n), , drop = FALSE]
        sc <- t(X) %*% solve(Soth) %*% u[t, ]
        Pi <- Pi + sc %*% t(sc)
        Q <- Q + t(X) %*% solve(Soth) %*% X
        eta <- solve(rsq(Sown)) %*% u[t, ]
        w <- as.vector(eta %*% t(eta) - diag(n))
        Om <- Om + w %*% t(w)
      }
      Ts <- length(rows)
      list(Pi = Pi / (Ts - q), Q = Q / Ts, Om = Om / (Ts - n^2))
    }
    if (!hetq && !brv) {
      L <- side(1:bigt, S1, S2)
      want <- c(sum(db * (L$Q %*% db)), sum(db * (L$Q %*% db)),
                sqrt(sum(db * (L$Pi %*% db))), sqrt(sum(db * (L$Pi %*% db))))
    } else {
      L <- side(1:br, S1, S2); Rr <- side((br + 1):bigt, S2, S1)
      want <- c(0.5 * sum(diag(A1 %*% A1)) + sum(db * (L$Q %*% db)),
                0.5 * sum(diag(A2 %*% A2)) + sum(db * (Rr$Q %*% db)),
                sqrt(0.25 * sum(as.vector(A1) * (L$Om %*% as.vector(A1))) + sum(db * (L$Pi %*% db))),
                sqrt(0.25 * sum(as.vector(A2) * (Rr$Om %*% as.vector(A2))) + sum(db * (Rr$Pi %*% db))))
    }
    got <- sysbreak:::intpc(s$matz, res, n, q, S1, S2, A1, A2, 1, br, br + 1, bigt,
                            bigt, hetq = hetq, robust = "hc", brv = brv,
                            brbeta = TRUE, prewhit = FALSE, diffb = db)
    expect_equal(unname(unlist(got[c("Delt1", "Delt2", "Gam1", "Gam2")])), want,
                 tolerance = 1e-10)
  }
})

test_that("robust = \"hac\" reproduces the 1.1.0 hac = TRUE statistics (Application 1)", {
  skip_on_cran()
  hd <- sysbreak::housing; nn <- nrow(hd)
  dat <- data.frame(hd[-1, ], L1.dGDP = hd$dGDP[-nn], L1.dHP = hd$dHP[-nn],
                    L1.FFR = hd$FFR[-nn])
  fit <- suppressWarnings(sysbreak(cbind(dGDP, dHP, dI) ~ L1.dGDP + L1.dHP + L1.FFR,
                                   data = dat, m = 3, trim = 0.2, robust = "hac"))
  expect_equal(fit$testing$statistic, c(63.11702627, 132.73866800, 136.93722106),
               tolerance = 1e-7, ignore_attr = TRUE)
  expect_identical(fit$testing$statistic_type, "SupF")
  expect_equal(fit$br_selected$breaks, 1)
})

test_that("robust tests are for coefficient breaks only", {
  dat <- sim_het2()
  r <- collect_warnings(sysbreak(cbind(y1, y2) ~ x, data = dat, m = 1, trim = 0.2,
                                 brv = TRUE, robust = "hc"))
  expect_true(any(grepl("coefficient breaks only", r$warnings)))
  r <- collect_warnings(sysbreak(cbind(y1, y2) ~ x, data = dat, m = 1, trim = 0.2,
                                 coef_breaks = "none", brv = TRUE, robust = "hc"))
  expect_true(any(grepl("not available", r$warnings)))
  expect_true(all(is.na(r$value$testing$statistic)))
})

test_that("invalid values of robust are rejected", {
  dat <- sim_het2()
  expect_error(sysbreak(cbind(y1, y2) ~ x, data = dat, m = 1, trim = 0.2,
                        robust = "white"))
  expect_error(sysbreak(cbind(y1, y2) ~ x, data = dat, m = 1, trim = 0.2,
                        robust = TRUE), "replaces the logical")
})

test_that("the choice is stored and printed, and prewhit is used only with \"hac\"", {
  dat <- sim_het2()
  fit <- sysbreak(cbind(y1, y2) ~ x, data = dat, m = 1, trim = 0.2,
                  robust = "hc", prewhit = TRUE)
  expect_identical(fit$spec$robust, "hc")
  expect_false(fit$spec$prewhit)
  expect_output(print(fit), "robust = hc")
  fit0 <- sysbreak(cbind(y1, y2) ~ x, data = dat, m = 1, trim = 0.2, robust = "hc")
  expect_identical(fit$testing, fit0$testing)
  expect_identical(fit$estimation, fit0$estimation)
})

test_that("methods agree with the fit under robust = \"hc\"", {
  dat <- sim_het2()
  fit <- sysbreak(cbind(y1, y2) ~ x, data = dat, m = 1, trim = 0.2, robust = "hc")
  check_methods(fit)
  y <- as.matrix(dat[, c("y1", "y2")])
  fitv <- varbreak(y, lags = 1, m = 1, trim = 0.2, robust = "hc")
  check_methods(fitv)
})
