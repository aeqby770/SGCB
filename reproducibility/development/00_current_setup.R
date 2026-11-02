# Load the installed native build and the current public R entry point.
LIBRARY <- normalizePath(Sys.getenv("SGCB_LIBRARY",file.path(ROOT,"library_factor")),winslash="/",mustWork=TRUE)
.libPaths(c(LIBRARY,.libPaths()))
library(SGCB)
stopifnot(normalizePath(find.package("SGCB"),winslash="/")==file.path(LIBRARY,"SGCB"))
NATIVE_ID <- attr(SGCB:::hurdle_gg_cpp(matrix(0,0,4),c(0L,0L,1L,1L),1,1L),"native_source_id")
SOURCE_ID <- NATIVE_ID
REFERENCE <- new.env(parent=asNamespace("SGCB"))
sys.source(file.path(ROOT,"software","SGCB","R","differential.R"),envir=REFERENCE)
stopifnot(identical(body(REFERENCE$sgcbDE),body(SGCB::sgcbDE)))
