#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title Convert Markview Metabolomics file (alignment table) into 
#' a artMS compatible format
#' 
#' @description `artMS` enables the relative quantification of untargeted 
#' polar metabolites using the alignment table generated by Markview. 

#' MarkerView is an ABSciex software that supports the files 
#' generated by Analyst software (`.wiff`) used to run our specific mass 
#' spectrometer (ABSciex Triple TOF 5600+). 
#' It also supports `.t2d` files generated by the 
#' Applied Biosystems 4700/4800 MALDI-TOF. 

#' MarkerView software is used to align mass spectrometry data from several 
#' samples for comparison. Using the import feature in the software, `.wiff` 
#' files (also `.t2d` MALDI-TOF files and tab-delimited `.txt` mass spectra data 
#' in mass-intensity format) are loaded for retention time alignment. 
#' Once the data files are selected, a series of windows will appear wherein 
#' peak finding, alignment, and filtering options can be entered and selected. 
#' These options include minimum spectral peak width, minimum retention time 
#' peak width, retention time and mass tolerance, and the ability to filter 
#' out peaks that do not appear in more than a user selected number of samples.
#' 
#' `artmsConvertMetabolomics`` processes the markview file to enable 
#' QC analysis and relative quantification using the artMS functions
#' 
#' @param input_file (char) Markview input file
#' @param out_file (char) Output file name
#' @param id_file (char) KEGG database
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @return (text file) Outputs the converted output name
#' @keywords metabolomics, convert
#' @examples
#' # Testing that the arguments cannot be null
#' artmsConvertMetabolomics(input_file = NULL, 
#'                          out_file = NULL)
#' @export
artmsConvertMetabolomics <- function(input_file, 
                                out_file, 
                                id_file = NULL,
                                verbose = TRUE){
  
  if(missing(input_file))
    stop("Input file name <input_file> is missed")
  if(missing(out_file))
    stop("<out_file> is missed")
  
  if(is.null(input_file) & is.null(out_file)){
    return("Both <input_file> and <out_file> are required")
  }
  
  if(verbose) message(">> Reading in data from: ", input_file," ")
  
  x <- fread(input_file)
  tmp <- data.table::melt(data=x, id=c(seq_len(7)), 
                          variable.name="RawFile", 
                          value.name="Intensity")
  
  tmp[,'Row'] = NULL
  tmp[,'Index'] = NULL
  setnames(tmp, c(1,3), c('Modified.sequence', 'Retention time'))
  tmp$Modified.sequence <- gsub(' \\([0-9]+\\)', '', tmp$Modified.sequence)
  tmp$Proteins <- tmp$Modified.sequence
  tmp$Charge <- 1
  
  tmp <- as.data.frame(tmp, stringsAsFactors = FALSE)
  
  # annotate with  known id's
  if(!is.null(id_file)){
    ids <- read.delim(id_file, stringsAsFactors = FALSE, sep='\t')
    ids$KEGG <- gsub('\\\xca','',ids$KEGG)
    #### This could be optimized!!!!!
    for( i in seq_len(length(ids$KEGG)) ){  
      idx <- tmp$'m/z' %in% ids$m.z[i]
      tmp$Proteins[idx] = ids$KEGG[i]
    }
  }
  
  if(verbose) message(">> Writing out data to:", out_file, "... ")
  write.table(as.data.frame(tmp, stringsAsFactors = FALSE), 
              out_file, 
              row.names = FALSE, 
              col.names = TRUE, 
              sep ='\t', 
              quote = FALSE)
  if(verbose) message('>> CONVERSION COMPLETE! ')
}


