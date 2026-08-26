run_slow_tests <- function() {
  identical(
    tolower(Sys.getenv("BAYESICS_RUN_SLOW_TESTS")),
    "true"
  )
}

run_parallel_tests <- function() {
  identical(
    tolower(Sys.getenv("BAYESICS_RUN_PARALLEL_TESTS")),
    "true"
  )
}

run_is_tests <- function() {
  identical(
    tolower(Sys.getenv("BAYESICS_RUN_IS_TESTS")),
    "true"
  )
}
  