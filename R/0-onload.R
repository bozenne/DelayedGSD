### 0-onload.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 17 2022 (16:45) 
## Version: 
## Last-Updated: okt 11 2023 (13:54) 
##           By: Brice Ozenne
##     Update #: 3
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:
DelayedGSD.env <- new.env() # create a specific environment for the package

.onAttach <- function(lib, pkg="DelayedGSD") {
    desc <- utils::packageDescription(pkg)
    packageStartupMessage(desc$Package, " version ",desc$Version)
    DelayedGSD.options(reinitialise = TRUE) # generate .LMMstar-options when loading the package   
}

##----------------------------------------------------------------------
### 0-onload.R ends here
