
# ------------------------------------------------------------------------------
#' @title Converts the Protein ID column of the evidence 
#' file selected by the user to mod-site-specific notation: 
#' `ProteinID` to `ProteinID_AAnumber` notation
#'
#' @description It enables the modified-peptide specific quantification by
#' converting the Protein column of the evidence file selected by the user  
#' to an `ProteinID_AAnumber`notation. 
#' In this way, each of the modified peptides can be quantified
#' independently across conditions. 
#' 
#' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#' 
#' WARNING: we have detected a version of MaxQuant (>1.6.3.0) outputs a`
#' "Modified sequence" column of the evidence file that has two important 
#' changes for the annotation of phosphorylation:
#' - Uses `p` instead of `(ph)`
#' - The modified residue (i.e. `STY`) is the residue on the right of the `p`, 
#' instead of the residue to the left of `(ph)`, as usual.
#' We have introduced a modification to detect and address this issue, but
#' we advice the user to double check both the new evidence file with the
#' introduce new notation and the `-mapping.txt` file and check that there
#' are no NA values for the notation of phophopeptides.
#' 
#' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#' 
#' @param evidence_file (char) The evidence file name and location
#' @param column_name (char) The Protein Column Name to map. Options:
#' - `Leadind razor protein` (default)
#' - `Leading protein`
#' - `Proteins`
#' It only supports Uniprot Entry IDs and RefSeq, but it might work for
#' other database IDs
#' @param ref_proteome_file (char) The reference proteome used as database
#' to search the `evidence.txt` file with MaxQuant. It will be used to map the
#' modified peptide to the protein sequence and find the site location.
#' Therefore, it does not use the MaxQuant's `Phospho (STY)Sites.txt`
#' @param output_file (char) Output file name 
#' (`ptmsites-evidence.txt` recommended)
#' @param overwrite_evidence (logical) if <output_file> is the same 
#' as <evidence_file>, `overwrite_evidence = FALSE` (default) doesn't allow to
#' overwrite the evidence file. Otherwise, `overwrite_evidence = TRUE` allows
#' to overwrite the evidence_file (this option might be activated if the user
#' allows to use the same `ptm-sites-evidence.txt` file to re-annotate all
#' the Protein IDs columns)
#' @param mod_type (char) The posttranslational modification. Options:
#' - `UB`: Protein Ubiquitination
#' - `PH`: Protein Phosphorylation
#' - `AC`: Protein Acetylation
#' @param verbose (logical) `TRUE` (default) shows function messages
#' @param keep_all_rows (logical) should unmodified and un-mapped rows be kept in the output
#' @param label_unmod_sites (logical) groups basd on modified as well as unmodified
#'  (if ever found modified) sites.  Example output (see return) now includes `A34890_S3;A34890_S5` and
#'   `A34890_S3;A34890_nS5`. The `n` denotes no mod on `S5` and will only be labeled
#'   as such if site 5 is ever found phosphorylated in the whole dataset.
#' @return (file) Return a new evidence file with the specified Protein id 
#' column modified by adding the sequence site location(s) + postranslational
#' modification(s) to the uniprot entry / refseq id.
#' 
#' Output ID examples: `A34890_ph3`; `Q64890_ph24_ph456`; 
#' `Q64890_ub34_ub129_ub234`; `Q64890_ac35`.
#' @keywords evidence, convert, ptm, ph, ub, ac
#' @examples
#' # Testing warning if files are not submitted. 
#' artmsProtein2SiteConversion(evidence_file = NULL, ref_proteome_file = NULL, 
#' output_file = NULL)
#' @export
artmsProtein2SiteConversion <- function (evidence_file,
                                         ref_proteome_file,
                                         column_name = c('Leading razor protein', 
                                                         'Leading proteins', 
                                                         'Proteins'),
                                         output_file,
                                         mod_type,
                                         overwrite_evidence = FALSE,
                                         verbose = TRUE,
                                         keep_all_rows = FALSE,
                                         label_unmod_sites=TRUE) {
  
  if(is.null(evidence_file) & 
     is.null(ref_proteome_file) & 
     is.null(output_file)){
    return("Files must not be NULL")
  }
  
  if(any(missing(evidence_file) | 
         missing(ref_proteome_file) |
         missing(output_file) | 
         missing(mod_type)))
    stop("Missed (one or many) required argument(s)
         Please, check the help of this function to find out more")
  
  mod_type <- toupper(mod_type)
  if(!mod_type %in% c("PH", "UB", "AC"))
    stop("the mod_type ", mod_type, " is not supported")
  
  # CHECK PROTEIN COLUMN
  column_name <- match.arg(column_name)
  
  # When the evidence file is not a data.frame, check that it does not have
  # the same as the output file
  if(is.vector(evidence_file)){
    if(!overwrite_evidence){
      if(evidence_file == output_file) 
        stop("<output_file> cannot be the same as <evidence_file>. 
If you are confident about overwritting the evidence file, 
then make the argument 'overwrite_evidence = TRUE'")
    }
  }
  
  if(verbose){
    message("-------------------------------------------------------------")
    message("artMS: Annotate Protein ID column to ProteinID_Site notation ")
    message("-------------------------------------------------------------")
  }
  
  if(verbose)
    message(">> FILTERING MODIFIED PEPTIDES")
  
  if (mod_type == 'UB') {
    if(verbose) message('--- SELECTING << UB >> MODIFIED PEPTIDES ')
    maxq_mod_residue <- 'K\\(gl\\)'
    mod_residue <- 'K'
  } else if (mod_type == 'PH') {
    if(verbose) message('--- SELECTING << PH >> MODIFIED PEPTIDES ')
    maxq_mod_residue <- '(S|T|Y)\\(ph\\)'
    mod_residue <- 'S|T|Y'
  } else if (mod_type == 'AC') {
    if(verbose) message('--- SELECTING << AC >> MODIFIED PEPTIDES ')
    maxq_mod_residue <- 'K\\(ac\\)'
    mod_residue <- 'K'
  } else{
    stop(
      mod_type, " is not supported. 
      Check help to get the list of supported PTMs"
    )
  }
  
  ## map mod sites in data to index
  if(verbose) message("--- OPENING EVIDENCE FILE ")
  ## read in maxq. data
  maxq_data <- .artms_checkIfFile(evidence_file, dont_check_names = TRUE)
  maxq_data <- as.data.table(maxq_data)
  
  # Check maxquant version:
  if("Leading Razor Protein" %in% colnames(maxq_data))
    setnames(maxq_data, 'Leading Razor Protein', 'Leading razor protein')
  
  if("Leading Proteins" %in% colnames(maxq_data))
    setnames(maxq_data, 'Leading Proteins', 'Leading proteins')
  
  # remove contaminants, keep unique sequences, fix names
  maxq_data <- maxq_data[grep("CON__|REV__", maxq_data$Proteins, invert = TRUE), ]
  maxq_data <- maxq_data[grep("CON__|REV__", maxq_data$`Leading proteins`, invert = TRUE), ]
  
  # In case the evidence file comes as data.frame
  if("Modified sequence" %in% colnames(maxq_data)){
    if(verbose) message("---(+) 'Modified sequence' column available") 
  }else if("Modified.sequence" %in% colnames(maxq_data)){
    setnames(maxq_data, 'Modified.sequence', 'Modified sequence')  
  }else{
    stop("Problem with evidence file: <Modified sequence> column not found")
  }
  
  # Fix PROTEIN id when coming with the full uniprot id
  
  if(any(grepl("(^sp\\|)(.*)(\\|.*)", maxq_data$Proteins)))
    maxq_data$Proteins <- gsub("(sp\\|)(.*)(\\|.*)", "\\2", maxq_data$Proteins)
  
  if(any(grepl("(^sp\\|)(.*)(\\|.*)", maxq_data$`Leading proteins`)))
    maxq_data$`Leading proteins` <- gsub("(sp\\|)(.*)(\\|.*)", "\\2", maxq_data$`Leading proteins`)
  
  if(any(grepl("(^sp\\|)(.*)(\\|.*)", maxq_data$`Leading razor protein`)))
    maxq_data$`Leading razor protein` <- gsub("(sp\\|)(.*)(\\|.*)", "\\2", maxq_data$`Leading razor protein`)
  
  if(any(grepl("(^tr\\|)(.*)(\\|.*)", maxq_data$Proteins)))
    maxq_data$Proteins <- gsub("(tr\\|)(.*)(\\|.*)", "\\2", maxq_data$Proteins)
  
  if(any(grepl("(^tr\\|)(.*)(\\|.*)", maxq_data$`Leading proteins`)))
    maxq_data$`Leading proteins` <- gsub("(tr\\|)(.*)(\\|.*)", "\\2", maxq_data$`Leading proteins`)
  
  if(any(grepl("(^tr\\|)(.*)(\\|.*)", maxq_data$`Leading razor protein`)))
    maxq_data$`Leading razor protein` <- gsub("(tr\\|)(.*)(\\|.*)", "\\2", maxq_data$`Leading razor protein`)
  
  # PH: Check if the protein has already been annotated;
  if( length( grep("_S\\d", maxq_data[[column_name]]) ) > 10 ) 
    stop("The protein column <", column_name, "> seems to be already converted 
to ptm-site/peptide specific notation.
Otherwise, notice that this function only support Uniprot Entry Id or Refseq")
  
  # UB: check annotated already
  if( length( grep("_K\\d", maxq_data[[column_name]]) ) > 10 ) 
    stop("The protein column <", column_name, "> seems to be already converted 
to ptm-site/peptide specific notation.
Otherwise, notice that this function only support Uniprot Entry Id or Refseq.
If the proteins are still Uniprot Entry IDs and the file has not been converted already, please, contact artMS developers")
  
  
  if (mod_type == 'PH') {
    if(length(grep("(ph)", maxq_data$`Modified sequence`)) == 1){
      message("\n____________________________________________________________")
      message("| WARNING: ")
      message("| The version of MaxQuant that generated this evidence file")
      message("| introduced changes in the <Modified Sequence> column that requires some changes.")
      message("| Basically, the phosphorylation site is indicated as 'XXXpSXXX' instead of 'XXXS(ph)XXX> as before.")
      message("| Next: Applyng changes to return to the previous version")
      message("____________________________________________________________\n")
      maxq_data$`Modified sequence` <- gsub("pS", "S(ph)", maxq_data$`Modified sequence`)
      maxq_data$`Modified sequence` <- gsub("pT", "T(ph)", maxq_data$`Modified sequence`)
      maxq_data$`Modified sequence` <- gsub("pY", "Y(ph)", maxq_data$`Modified sequence`)
    }
  }
  
  if(verbose) message("--- READING REFERENCE PROTEOME ")
  ## read in reference proteome
  ref_proteome <- read.fasta(file = ref_proteome_file,
                             seqtype = "AA",
                             as.string = TRUE,
                             set.attributes = TRUE,
                             legacy.mode = TRUE,
                             seqonly = FALSE,
                             strip.desc = FALSE)
    
  
  ## make mod-site index
  p_seqs <- c()
  p_names <- c()
  p_annots <- c()
  
  for (e in ref_proteome) {
    p_seqs <- c(p_seqs, e[1])
    p_names <- c(p_names, attr(e, 'name'))
    p_annots <- c(p_annots, attr(e, 'Annot'))
  }
  
  ref_table <- data.table(names = p_names,
                          annots = p_annots,
                          seqs = p_seqs)
    
  
  # REFSEQ database
  ref_table[, uniprot_ac := gsub('([a-z,0-9,A-Z,\\.,\\_]+)', '\\1', names)]
  
  # UNIPROT
  ref_table[, uniprot_ac := gsub('([a-z,0-9,A-Z]+\\|{1})([A-Z,0-9,\\_]+)(\\|[A-Z,a-z,0-9,_]+)',
                                 '\\2', names)]
            
  
  # get all indicies/locations of the mod_residue S|T|Y
  indices <- lapply(ref_table$seqs, function(x)
    as.vector(str_locate_all(x, pattern = mod_residue)[[1]][, 1]))
    
  # list which residue it actually finds at each location
  ptm_sites <- lapply(ref_table$seqs, function(x)
    as.vector(str_match_all(x, pattern = mod_residue)))
    
  # get the number of residue locations per string
  lengths <- unlist(lapply(indices, FUN = length))
  # repeate the uniprot_ac as many times as a residue (S|T|Y) is found
  keys <- rep(ref_table$uniprot_ac, lengths)
  # combine the list of the exploded proteins (above) and add in which sites
  # are found as well as the location the residue was found in
  protein_indices <-   data.table(uniprot_ac = keys,
                                  ptm_site = unlist(ptm_sites),
                                  res_index = unlist(indices))
  setkey(protein_indices, uniprot_ac, res_index) # for fast subsetting by keys below
    
  
  
  # ---------------------------------------------------------------------------
  # EXTRACT PTM POSITIONS FROM MODIFIED PEPTIDES IN EVIDENCE FILE
  # ---------------------------------------------------------------------------
  
  
  unique_peptides_in_data <- unique(maxq_data[, c(column_name, 'Modified sequence'), with = FALSE])
    
  setnames(unique_peptides_in_data, 'Modified sequence', 'sequence') 
  
  #these track modifications
  mod_sites <- c()
  mod_seqs <- c()
  
  #this tracks all, modified or not, only used by label_unmod_sites == TRUE
  all_sites <- list()
  
  if(verbose) message("--- EXTRACTING PTM POSITIONS FROM THE MODIFIED PEPTIDES 
   (This might take a long time depending on the size of the fasta/evidence file) ")
  pb <- txtProgressBar(min = 0, max = nrow(unique_peptides_in_data), style=3)
  for (i in seq_len(nrow(unique_peptides_in_data))) {
    setTxtProgressBar(pb,i)
    entry <- unique_peptides_in_data[i, ]
    peptide_seq <- entry$sequence
    ## cleanup the sequence (removing all modifications) for matching the 
    ## protein sequence
    peptide_seq_clean <- gsub('[a-z,0-9,\\(,\\),_]', '', peptide_seq)
    mod_sites_in_peptide <-
      str_locate_all(string = peptide_seq, pattern = maxq_mod_residue)[[1]][, 1]
    
    if (length(mod_sites_in_peptide) > 0) {
      uniprot_acs <- entry[[column_name]]
      # separates the ambiguous cases (;) and appends the site info to all 
      # the proteins
      uniprot_acs <- str_split(string = uniprot_acs, pattern = ';')[[1]]
      
      for (uac in uniprot_acs) {
        # find the protein's full sequence based on the uniprot_ac
        protein_seq <- ref_table[uniprot_ac == uac, ]$seqs
        if (length(protein_seq) > 0) {
          ## get the position of the peptide in the protein sequence
          peptide_index_in_protein <- str_locate(protein_seq, peptide_seq_clean)[[1]][1]
          
          if (label_unmod_sites){
            #for purposes of counting unmodified in same sequence, build a map of peptide
            # to all possible sites, and prepare to tag those modified in this peptide
            all_sites[[peptide_seq]][[uac]]<- protein_indices[.(uac, 
                                                                peptide_index_in_protein:(peptide_index_in_protein+str_length(peptide_seq_clean)))  #subset using keys
                                                              ][!is.na(uac)   #remove the NA rows we get by asking for res_index  keys that don't exist
                                                                ][,is_modified:=FALSE] # initialize to FALSE; will set to TRUE below after locating mods
          }
          
          
          for (m in seq_len(length(mod_sites_in_peptide))) {
            mod_site <- mod_sites_in_peptide[m]
            peptide_seq_before_site <-
              str_sub(peptide_seq, 1, mod_site - 1)
            ## count all AA (not counting all modifications) before the
            ## modification to get the relative position of the modification
            ## in the peptide sequence
            residues_before_site <-
              str_count(string = peptide_seq_before_site, 
                        pattern = 'A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y')
            mod_site_index_in_protein <-
              peptide_index_in_protein + residues_before_site
            protein_mod_sites <- protein_indices[uniprot_ac == uac, ]
            if (!is.na(mod_site_index_in_protein)) {
              mod_res <-
                protein_mod_sites[res_index == mod_site_index_in_protein, 
                                  ptm_site]
              mod_site_id <-
                sprintf(
                  '%s_%s%s',
                  uac,
                  str_sub(
                    protein_seq,
                    mod_site_index_in_protein,
                    mod_site_index_in_protein
                  ),
                  mod_site_index_in_protein
                )
              mod_sites <- c(mod_sites, mod_site_id)
              mod_seqs <- c(mod_seqs, peptide_seq)
              if(label_unmod_sites)all_sites[[peptide_seq]][[uac]][mod_site_index_in_protein == res_index,
                                                                   is_modified :=TRUE]
              # message("mod_site_id ", mod_site_id, " peptide_seq ",peptide_seq)
              # stopifnot(length(mod_sites) == length(mod_seqs))
            } else{
              message(
                sprintf(
                  'MISMATCH\t%s\n\tPEPTIDE_SEQ\t%s\n\tMOD_SITE\t%s\n\t
                  PEPTIDE_IDX_IN_PROTEIN\t%s\n\tRESIDUES_BEFORE_SITE
                  \t%s\n\tPROTEIN_SEQ\t%s\n',
                  mod_site_id,
                  peptide_seq,
                  mod_site,
                  peptide_index_in_protein,
                  residues_before_site,
                  protein_seq
                )
              )
            }
          }
        }else{
          message("- Protein <", uac, "> not in the sequence database")
        }
      }
    }
  }
  close(pb)
  
  # CHECK that the mapping went well
  if(is.null(mod_seqs))
    stop("Protein IDs from evidence file and sequence database do not match. 
         Are you sure that you are using the right sequence database?")
  
  if(is.null(mod_sites))
    stop("Protein IDs from evidence file and sequence database do not match. 
         Are you sure that you are using the right sequence database?")
  
  if (!label_unmod_sites){
    mod_site_mapping <- data.table(mod_sites, mod_seqs)
  }else{
    # We need to keep unmodified sites in our labeles, but only for those
    # sites that we ever detect modified.
    # Example AAAS(ph)LLLSR should match both mod site protID_S4 and protID_nS8
    # if we ever have protID_S8 in this table
    
    #assemble table of all sites from the nested list of lists of data.tables
    all_sites_table <- rbindlist(lapply(all_sites, function(x)rbindlist(x,idcol="protein")),
                                 idcol="mod_seqs")
    #add convenient key columns for mapping to our list of detected modifications
    all_sites_table[,mod_site_key := sprintf('%s_%s%s', protein, ptm_site, res_index )]
    #reduce to just those sites we ever find modiifed
    all_sites_table <- all_sites_table[mod_site_key %in% mod_sites]
    all_sites_table[,unmod_site_key := sprintf('%s_n%s%s', protein, ptm_site, res_index )]
    # get the correct label depending on mod state
    all_sites_table[, mod_site_label := ifelse(is_modified, mod_site_key, unmod_site_key)]
    
    mod_site_mapping <- all_sites_table[,.(mod_sites = mod_site_label, mod_seqs = mod_seqs)]
  }
  
  mod_site_mapping_agg <- aggregate(mod_sites ~ mod_seqs,
                                    mod_site_mapping,
                                    FUN = function(x) paste(x, collapse = ';'))
    
  setnames(maxq_data, 'Modified sequence', 'mod_seqs')
  unmapped_mod_seqs <- maxq_data[!(mod_seqs %in% mod_site_mapping_agg$mod_seqs) &
                                   grepl('(gl)', mod_seqs) & !grepl('REV__|CON__', eval(column_name)), ]
    
  unmapped_mod_seqs <- unique(unmapped_mod_seqs[, c('mod_seqs', column_name), with = FALSE])
  
  if (dim(unmapped_mod_seqs)[1] > 0) {
    if(verbose) message('>> UNABLE TO MAP \t')
    if(verbose) print(unmapped_mod_seqs)
  } else{
    if(verbose) message("--- ALL SEQUENCES MAPPED ")
  }
  
  final_data <- merge(maxq_data, mod_site_mapping_agg, by = 'mod_seqs', all.x=keep_all_rows)
  ref_protein_column <- paste0(column_name, " Ref")
  setnames(
    final_data,
    c(column_name, 'mod_sites', 'mod_seqs'),
    c(ref_protein_column, column_name, 'Modified sequence')
  )
  
  final_data <- as.data.frame(final_data)
  final_data <- final_data %>% dplyr::select(column_name, everything())
  
  write.table(
    final_data,
    file = output_file,
    eol = '\n',
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  
  ## write a mapping table
  protein_seq_mapping <-
    unique(maxq_data[, c(column_name, 'mod_seqs'), with = FALSE])
  # setnames(protein_seq_mapping, eval(column_name), 'Protein')
  mapping_table <- merge(protein_seq_mapping,
                         mod_site_mapping_agg,
                         by = 'mod_seqs',
                         all = TRUE)
  write.table(
    mapping_table,
    file = gsub('.txt', '-mapping.txt', output_file),
    eol = '\n',
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  
  if(verbose) 
    message(
      ">> FILES OUT:\n",
      "--- New evidence-site file: ",
      output_file,
      "\n",
      "--- Details of the Mappings: ",
      gsub('.txt', '-mapping.txt', output_file)
    )
  if(verbose) message(">> CONVERSION COMPLETED")
}
