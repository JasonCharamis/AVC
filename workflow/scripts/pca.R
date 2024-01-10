## Principal components analysis based on SNPs identified through DNA-seq ##
## Adopted pipeline from here: https://rpubs.com/madisondougherty/980777 ##


dependencies <- c("optparse", "vcfR", "vegan", "tidyverse", "edgeR", "ggfortify", "ggrepel", "ggthemes", "ggplot2")

load_and_suppress <- function ( pkg ) {
		      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
		    }

for ( pkg in dependencies ) {
  load_and_suppress( pkg )
}

# Define the option parser
option_list <- list(
  make_option(
    c("-i", "--input-vcf"),
    type = "character",
    help = "Path to input VCF file."
  )
)


# Parse the command-line arguments
opt_parser <- OptionParser(usage = "Usage: %prog -i INPUT_FILE", option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the input file option is provided
if (is.null(opt$input)) {
  stop("Error: Input file option is required. Use -i or --input-vcf to specify the input file.")
}

# Check if the input file exists
if (!file.exists(opt$input)) {
  stop(paste("Error: The file", opt$input, "does not exist."))
}



## OPEN VCF FILE AND EXTRACT GENOTYPES USING THE vcfR PACKAGE ##

snps <- vcfR::read.vcfR(opt$input, convertNA  = TRUE)

snps_num <- vcfR::extract.gt(snps, 
	    		     element = "GT",
			     IDtoRowNames  = F,
			     as.numeric = T,
			     convertNA = T,
			     return.alleles = F)

snps_num_t <- t(snps_num) 
snps_num_df <- data.frame(snps_num_t) 

find_NAs <- function(x){
	      NAs_TF <- is.na(x)
	      i_NA <- which(NAs_TF == TRUE)
	      N_NA <- length(i_NA)
	      
	      cat("Results:",N_NA, "NAs present\n.")
	    return(i_NA)
}

# N_rows
# number of rows (individuals)
N_rows <- nrow(snps_num_t)

# N_NA
# vector to hold output (number of NAs)
N_NA   <- rep(x = 0, times = N_rows)

# N_SNPs
# total number of columns (SNPs)
N_SNPs <- ncol(snps_num_t)

for(i in 1:N_rows) {
  
  # for each row, find the location of
  ## NAs with snps_num_t()
  i_NA <- find_NAs(snps_num_t[i,]) 
  
  # then determine how many NAs
  ## with length()
  N_NA_i <- length(i_NA)
  
  # then save the output to 
  ## our storage vector
  N_NA[i] <- N_NA_i
}

percent_NA <- N_NA/N_SNPs*100

# Call which() on percent_NA
i_NA_50percent <- which(percent_NA > 50) 

row_names <- row.names(snps_num_t)
row_names02 <- gsub(".bam","",row_names)

invar_omit <- function(x){
  cat("Dataframe of dim",dim(x), "processed...\n")
  sds <- apply(x, 2, sd, na.rm = TRUE)
  i_var0 <- which(sds == 0)
 
  cat(length(i_var0),"columns removed\n")
  
  if(length(i_var0) > 0){
     x <- x[, -i_var0]
  }
  
  ## add return()  with x in it
  return(x)                      
}


snps_no_invar <- invar_omit(snps_num_t) 
snps_noNAs <- snps_no_invar
N_col <- ncol(snps_no_invar)

for(i in 1:N_col){
  
  # get the current column
  column_i <- snps_noNAs[, i]
  
  # get the mean of the current column
  mean_i <- mean(column_i, na.rm = TRUE)
  
  # get the NAs in the current column
  NAs_i <- which(is.na(column_i))
  
  # record the number of NAs
  N_NAs <- length(NAs_i)

  # replace the NAs in the current column
  column_i[NAs_i] <- mean_i
  
  # replace the original column with the updated columns
  snps_noNAs[, i] <- column_i
  
}
 
if (is.null(snps_noNAs)) {
  stop("Error: Matrix has zero dimensions after invariant removal.")
  } else {
    print("Invariance removed successfully.")
}

write.csv(snps_noNAs, file = "SNPs_cleaned.csv", row.names = F)

SNPs_cleaned <- read.csv("SNPs_cleaned.csv")
pca_scaled <- prcomp(SNPs_cleaned)

## generate sample names for grouping replicates ##
xt <- as.data.frame(vegan::scores(pca_scaled))

groups <- rownames(xt) <- row_names02

## add column with sample name to group replicates ##
xtl <- xt %>% add_column(Sample = groups)

pdf("PCA.pdf")

## draw auto-PCA plot with color mappings for groups ##
print ( autoplot(pca_scaled, data=xtl, colour='Sample', legend.size=5) +  labs(title = "Principal Component Analysis" ) + theme_bw() + geom_text_repel(aes(label=rownames(xt),color=Sample), show.legend = FALSE  ) + theme(plot.title=element_text(face="bold",hjust=0.5), legend.title = element_text(size=12), legend.text = element_text(size=12), axis.title=element_text(size=12))   )

dev.off()

