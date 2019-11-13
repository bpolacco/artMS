# I propose this file as the collection of all data loading.  This will
# provide better code reusability and easier debugging of data-loading
# and formatting issues.  -Ben Polacco Nov/13/2019



# ------------------------------------------------------------------------------
# @title Load a keys file or return a validated keys data object if a table is passed in
# @description This function is used in order to make it so a user can submit
#  either a path to a data file or a data object in data.frame or data.table
#  form.
# @param input_file (object or data.frame) The filepath/keys object to be checked.
# @param fractions, usually from config$fractions$enabled (I think...)
# @return A valid keys table
# @keywords internal, file, keys, input
.artms_loadAndValidateKeys <- function(input_file, fractions = 0){
  # OPEN KEYS
  keys <- .artms_checkIfFile(input_file)
  keys <- .artms_checkRawFileColumnName(keys)

  if (fractions) {
     if (any(!'FractionKey' %in% colnames(keys))) {
      stop(' <fractions> was activated but <fractionkey> column not found in the keys data ')
    }
  }
  
  keys <- .artms_validateBioReplicateColumn(keys)
  
  return(keys)
}


.artms_validateBioReplicateColumn = function(keys){
  #Bioreplicate is expected to be a character vector. It may load as integer or numeric if in format 1 or 1.2
  if (class(keys$BioReplicate) != "character"){
    keys$BioReplicate = as.character(keys$BioReplicate)
  }
  
  #BioReplicate is expected to be distinct across conditions.
  # i.e. with conditions:  Endo, Lys and 2 bioreplicates, the BioReplicate column should be
  # c('Endo-1', 'Endo-2', 'Lys-1', 'Lys-2') and not c('1','2','1','2')
  numConditions = length(unique(keys$Condition))
  numCondPerBR =  unlist(lapply (split (keys$Condition, keys$BioReplicate), function(x){length(unique(x))}))
  if (any(numCondPerBR > 1)){
    # we need to apply a fix
    if (any(numCondPerBR != numConditions)){
      #slightly unexpected, different numbers of BR per condition.  issue a warning
      message ("While building BioReplicate labels, different numbers of replicates per condition were detected")
    }
    keys$BioReplicate = paste (keys$Condition, keys$BioReplicate, sep="-")
  }
  return (keys)
}
