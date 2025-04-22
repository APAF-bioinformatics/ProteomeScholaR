# R/zzz.R

.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        "-------------------------------------------------------------------\n",
        "Welcome to ProteomeScholaR!\n\n",
        "IMPORTANT: Please run ProteomeScholaR::loadDependencies() \n",
        "           to ensure all required packages are installed and loaded.\n",
        "-------------------------------------------------------------------"
    )
}
