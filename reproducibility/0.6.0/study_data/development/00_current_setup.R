# Load the installed native build and the current public R entry point.
LIBRARY <- normalizePath(Sys.getenv("SGCB_LIBRARY",file.path(dirname(ROOT),"library")),winslash="/",mustWork=TRUE)
.libPaths(c(LIBRARY,.libPaths()))
library(SGCB)
stopifnot(as.character(packageVersion("SGCB",lib.loc=LIBRARY))=="0.6.0")
stopifnot(normalizePath(find.package("SGCB"),winslash="/")==file.path(LIBRARY,"SGCB"))
NATIVE_ID <- attr(SGCB:::hurdle_gg_cpp(matrix(0,0,4),c(0L,0L,1L,1L),1,1L),"native_source_id")
SOURCE_ID <- NATIVE_ID
REFERENCE <- new.env(parent=asNamespace("SGCB"))
sys.source(file.path(dirname(ROOT),"software","SGCB","R","differential.R"),envir=REFERENCE)
stopifnot(identical(body(REFERENCE$sgcbDE),body(SGCB::sgcbDE)))

# Native IDs are retained author build labels, not content fingerprints.
# Stored fits are identified by the released record checksum; new fits use the current library.
verify_study_record <- function(path, source_id, record_id) {
    switch(as.character(record_id==as.character(packageVersion("SGCB",lib.loc=LIBRARY))),
        "TRUE"=stopifnot(identical(source_id,NATIVE_ID)),
        "FALSE"={
            index <- data.table::fread(file.path(dirname(ROOT),"results/r1_file_index.tsv"))
            expected <- index[index$file==paste0("study_data/results/",basename(path)),]
            stopifnot(nrow(expected)==1L,identical(expected$author_build_id,record_id),
                identical(expected$content_md5,unname(tools::md5sum(path))))
        })
    invisible(source_id)
}
