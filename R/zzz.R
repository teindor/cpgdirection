.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        "cpgdirection ", utils::packageVersion("cpgdirection"), "\n",
        "  cpg_expression_direction(my_cpgs)   all methods; see best_direction\n",
        "  cpgd_selftest()                     verify the installation\n",
        "  Predictions are conditional on the CpG being an eQTM; see ?cpg_expression_direction"
    )
}
