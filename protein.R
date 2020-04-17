
library(readxl)

getProtein <- function(){

protein <- read_excel(path="NIRdata.xlsx",
                      sheet="Sorted Barcodes",
                      range="G1:M100")

nirkeys1 <- read_excel(path="NIR.xlsx",
                      range="A1:B51")

nirkeys2 <- read_excel(path="NIR.xlsx",
                       range="D1:E45")


nirkeys <- rbind(nirkeys1, nirkeys2)
names(nirkeys) <- c( "ID" , "Sample ID" )

protein <- merge( protein, nirkeys, all = TRUE, by="Sample ID")
protein <- tidyr::separate( protein , "ID", c("Cross", "ID", "Genotype"), " +")
protein
}

protein <- getProtein()
