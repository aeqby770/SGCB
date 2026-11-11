root <- normalizePath(Sys.getenv("SGCB_R2_ROOT"),mustWork=TRUE,winslash="/")
lib <- file.path(root,"library"); stopifnot(as.character(packageVersion("SGCB",lib.loc=lib))=="0.6.0")
target <- file.path(root,"release","SGCB_0.6.0.zip")
zip::zipr(target,files="SGCB",root=lib,include_directories=TRUE,compression_level=7)
test <- file.path(root,"tmp","windows_install_check");dir.create(test,showWarnings=FALSE)
install.packages(target,repos=NULL,type="win.binary",lib=test)
library(SGCB,lib.loc=test);stopifnot(as.character(packageVersion("SGCB",lib.loc=test))=="0.6.0",
    "endpoint"%in%names(formals(sgcbDE)))
cat("Windows binary installed successfully:",find.package("SGCB"),"\n")