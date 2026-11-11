root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE); lib <- file.path(root,"library")
.libPaths(c(lib,.libPaths()));library(SGCB,lib.loc=lib)
stopifnot(as.character(packageVersion("SGCB",lib.loc=lib))=="0.6.0")
tests <- testthat::test_dir(file.path(root,"software/SGCB/tests/testthat"),reporter="summary",package="SGCB",
    load_package="installed",stop_on_failure=TRUE,stop_on_warning=FALSE)
stats <- as.data.frame(tests)
data.table::fwrite(stats[,intersect(names(stats),c("file","context","test","nb","failed","skipped","error","warning","passed")),drop=FALSE],
    file.path(root,"results","package_tests.tsv"),sep="\t")