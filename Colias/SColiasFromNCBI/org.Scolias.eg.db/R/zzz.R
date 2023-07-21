datacache <- new.env(hash=TRUE, parent=emptyenv())

org.Scolias.eg <- function() showQCData("org.Scolias.eg", datacache)
org.Scolias.eg_dbconn <- function() dbconn(datacache)
org.Scolias.eg_dbfile <- function() dbfile(datacache)
org.Scolias.eg_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
org.Scolias.eg_dbInfo <- function() dbInfo(datacache)

org.Scolias.egORGANISM <- "Scomber colias"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "org.Scolias.eg.sqlite", package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)

    ## Create the OrgDb object
    sPkgname <- sub(".db$","",pkgname)
    db <- loadDb(system.file("extdata", paste(sPkgname,
      ".sqlite",sep=""), package=pkgname, lib.loc=libname),
                   packageName=pkgname)    
    dbNewname <- AnnotationDbi:::dbObjectName(pkgname,"OrgDb")
    ns <- asNamespace(pkgname)
    assign(dbNewname, db, envir=ns)
    namespaceExport(ns, dbNewname)
        
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("org.Scolias.eg.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(org.Scolias.eg_dbconn())
}

