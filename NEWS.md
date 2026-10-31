# sysbreak (version 1.1.0)

* `coef()` now returns a named vector in the same order and with the same names as `vcov()` (e.g. `regime1:y1:x`); the per-regime tables remain in `summary()`.

* New argument `robust = c("none", "hc", "hac")` in `sysbreak()`, `varbreak()` and `netbreak()` replaces the logical `hac`; `"hac"` is the former `hac = TRUE`. `"hc"` computes the robust SupF tests, coefficient standard errors and break-date confidence intervals with a heteroskedasticity-consistent covariance estimate without autocovariance terms.

* The automatic HAC bandwidth is capped at T / (2d), where d is the number of score columns, so that the HAC estimate does not become near singular in short segments.

* The internal pseudo-inverse discards only numerically zero singular values (relative tolerance 1e-14 instead of the square root of machine epsilon). Results no longer depend on the units of the variables; the Application 1 estimates and statistics change (the dI-equation coefficients, and the GDP break from 1983:Q1 to 1982:Q4).

* `print(summary(fit))` has a new layout. The sequential tests are labelled `SEQ(l+1|l)`, as in Qu and Perron (2007), also in the row names of `summary(fit)$tests`.
