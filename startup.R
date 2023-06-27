rm(list = ls()); graphics.off(); cat("\014")

os = Sys.info()['sysname'] 

# load libraries
if(os == 'Linux') 
{    
    libs = c('colorout', 'parallel', 'R.matlab', 'statmod', 'sfsmisc', 'fields', 'JOPS', 'matrixStats', 'scales', 'ggplot2') 
    install.packages(setdiff(libs, rownames(installed.packages())))  
    pkgLoad = lapply(libs, function( x ) require(x, character.only = T))    
} 

if(os == 'Windows')
{
    libs = c('R.matlab', 'statmod', 'sfsmisc', 'fields', 'JOPS', 'devtools', 'matrixStats', 'scales', 'parallelsugar', 'ggplot2') 
    install.packages(setdiff(libs, rownames(installed.packages())))  
    winDiff = setdiff(libs, rownames(installed.packages()))
    if((winDiff == 'parallelsugar') & (('parallelsugar' %in% setdiff(libs, rownames(installed.packages()))))) 
    {  
        devtools::install_github('nathanvan/parallelsugar')
    }
    pkgLoad = lapply(libs, function(x) require(x, character.only = T))
}


