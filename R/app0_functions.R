# Write code to install all the packages 
### Function to check if packages are installed and if not, install them
install_and_load <- function(packages) {
  # Modified from https://www.listendata.com/2018/12/install-load-multiple-r-packages.html
  # The function below performs the following operations -
  #  - First it finds all the already installed R packages
  #  - Check packages which we want to install are already installed or not.
  #  - If package is already installed, it does not install it again.
  #  - If package is missing (not installed), it installs the package.
  #  - Loop through steps 2, 3 and 4 for multiple packages we want to install
  #  - Load all the packages (both already available and new ones).
  k <- packages[!(packages %in% installed.packages()[, "Package"])]
  if(length(k)){
    install.packages(k, 
                     repos = "https://cran.rstudio.com/", 
                     dependencies = TRUE)
  }
  
  for(package_name in packages){
    library(package_name,
            character.only = TRUE, 
            quietly = TRUE)
  }
}






