# Complete R Script: Plasmid Alignment and Variant Analysis
# includes optional Macrogen correction (that must be verified by Sanger sequencing)
# Developed by: Michael J Courtney, Turku Bioscience Centre, University of Turku
# Tested by: as above, add your name here 
# Funding : SynGAP Research Fund, Research Council of Finland, EU-IMPULSE/EU Horizon Europe

input_path <- "C:/Seafile/LabLibrary2/Plasmids/Sequences&plasmids/_FPSanalysis/Nanoluc_SynGAP"
SampleListName <- "SeqToAnalyse_Nluc_v2.txt"
# expects to be able to find the following columns, any order
# SequenceID, MacrogenCorrection, Filename, Path


#some constants
orf_start <- 613 # should be in the ref seq name? - this is for nanoluc set
tag_offset <- 186 -14 #  number of codons of the tag compensated by our A1 start site so that SynGAP p. are correct
Macrogen_insert_pos <-288+14*3+1 # with reference to start of SynGAP 
Macrogen_base_to_insert <- "G" # if Macrogen sequenced always add this

#*** USAGE ***
#Automated generation of systematic reports that are easy to read and check. 

# intended for when we have multiple full sequence results from >=1 sequencing runs and the differences are supposed to be a single base (or a few) or short deletions/inserts. 
# needs a file listing paths and filenames of sequence data, an ID and indicator if we need to correct for the Macrogen misread (always needs confirming by Sanger seq, but that’s easy)
#*** that's the SampleListName and input_path above ***###
# The first entry should be a reference, the wild type plasmid sequence as it should be
# The script needs to know where the ORF of interest starts in the reference and size of the tag so that the reported mutation positions correspond to the insert
# It will list all differences as protein (amino_acid_variants.tsv) and, including outside the coding sequence, as nucleotide (nucleotide_variants.tsv). The numbering will be for the protein of interest as well as the location in the plasmid and full alignments across the dataset (protein and nucleotide) are added to the headers of all pdw generated. 
# The pdw are saved in the same folder where the list of plasmids is specified as well as in the original sequencing folders.
# If one sequence is completely different or has major problems, this messes up the alignment, but it is easy to remove that from the list and run again, and maybe run this one separately if needed.

#** EXAMPLE SampleListName FILE **#
#SequenceID,MacrogenCorrection,Filename,Path
#Nanoluc_Reference,template.fasta,C:\Seafile\LabLibrary2\Plasmids\Sequences&plasmids\_FPSanalysis	
#36,yes,consensus.fasta,C:\Seafile\LabLibrary2\Plasmids\Sequences&plasmids\_Sequencing Results\wholeplasmidsequencing\MacrogenResults\EN00004042\Results\6M


options(repos = BiocManager::repositories())

install_if_missing <- function(pkgs) {
  to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(to_install)) install.packages(to_install)
}
install_bioc_if_missing <- function(pkgs) {
  to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(to_install)) BiocManager::install(to_install)
}

cran_packages <- c("BiocManager", "remotes",
                   "stringr",
                   "seqinr", "ape")
                   # "moments","rlang",
                   # "ggplot2",
                   # "viridis", "RColorBrewer",
                   # "scales",  "pheatmap",
                   # "ComplexHeatmap","circlize",
                   # "data.table", "purrr",
                   # "MASS","caret", "glmnet",#"pROC",
                   # "FactoMineR", "factoextra",
                   # #"Metrics", "randomForest", "class", "infotheo",
                   # #"ijtiff","stringr", 
                   # "tidyr", "dplyr")

bioc_packages <- c("Biostrings","BiocGenerics", "msa") #,"ropls") 

install_if_missing(cran_packages)
install_bioc_if_missing(bioc_packages)

#Load each package
invisible(lapply(unique(c(cran_packages, bioc_packages)), function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Package '%s' is not installed.", pkg))
  } else {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}))



#FUNCTIONS

# Function to rotate sequence to best match reference - note align to orf start or indels will get us lost 
rotate_to_match <- function(seq, ref, window_start = 1, window_size = 50) {
  seq_str_fwd <- as.character(seq)
  seq_str_rev <- as.character(Biostrings::reverseComplement(seq))
  ref_str <- as.character(ref)
  ref_window <- substr(ref_str, window_start, window_start + window_size - 1)
  
  best_score <- -Inf
  best_offset <- 0
  best_orientation <- "forward"
  best_rotated <- seq_str_fwd
  
  for (orientation in c("forward", "reverse")) {
    seq_str <- if (orientation == "forward") seq_str_fwd else seq_str_rev
    
    for (i in 0:(nchar(seq_str) - 1)) {
      rotated <- paste0(substr(seq_str, i + 1, nchar(seq_str)), substr(seq_str, 1, i))
      rotated_window <- substr(rotated, 1, window_size)
      score <- sum(strsplit(rotated_window, "")[[1]] == strsplit(ref_window, "")[[1]])
      
      if (score > best_score) {
        best_score <- score
        best_offset <- i
        best_orientation <- orientation
        best_rotated <- rotated
        best_rotated_aligned_to_start <- paste0(substr(rotated, nchar(rotated) - window_start + 2, nchar(rotated)), substr(rotated, 1, nchar(rotated) - window_start + 1))
        
      }
    }
  }
  
  cat("Best orientation:", best_orientation, 
      "| Rotation offset:", best_offset, 
      "| Match score:", best_score, 
      "| Rotated window:", substr(best_rotated, 1, window_size),
      "| Aligned_to_start:",substr(best_rotated_aligned_to_start, 1, window_size),
      "\n")
  
  

  
  
  return(DNAString(best_rotated_aligned_to_start))
}

# Function to extract ORF and translate

extract_orf <- function(seq, start, id = "current") {
  
  orf_ignore_stop <- subseq(seq, start = start)
  len <- length(orf_ignore_stop)

  orf_ignore_stop_trimmed <- subseq(orf_ignore_stop, start = 1, width = len - (len %% 3))
  aa <- Biostrings::translate(orf_ignore_stop_trimmed)
  # Truncate at first stop codon
  stop_pos <- regexpr("\\*", as.character(aa))[1]
  if (stop_pos > 0) {
    aa <- subseq(aa, start = 1, end = stop_pos)
  }
  orf<- subseq(orf_ignore_stop, start = 1, width = length(aa)*3)
  
  # Report lengths
  cat("SequenceID:", id, "| ORF length (nt):", length(orf), "| AA length:", length(aa), "\n")
  
  return(aa)
}




#write clustal 
write_clustal <- function(aln, file) {
  con <- file(file, "w")
  on.exit(close(con))
  
  writeLines("CLUSTAL W (1.82) multiple sequence alignment\n", con)
  
  max_name_len <- max(nchar(aln$nam))
  block_size <- 60
  nseq <- aln$nb
  seqlen <- nchar(aln$seq[[1]])
  
  for (start in seq(1, seqlen, by = block_size)) {
    end <- min(start + block_size - 1, seqlen)
    for (i in 1:nseq) {
      name <- format(aln$nam[i], width = max_name_len)
      seq_block <- substr(aln$seq[[i]], start, end)
      writeLines(paste(name, seq_block), con)
    }
    writeLines("", con)
  }
}



#write clustal alignment with numbers
write_custom_alignment <- function(alignment, file) {
  seqs <- alignment$seq
  names <- alignment$nam
  nseqs <- length(seqs)
  seqlen <- nchar(seqs[[1]])
  block_size <- 60
  
  # Determine max name width for alignment
  name_width <- max(nchar(names)) + 1  # +1 for spacing
  
  # Track cumulative non-gap positions
  cum_pos <- rep(0, nseqs)
  
  # Open file connection
  con <- file(file, "w")
  on.exit(close(con))
  
  # Write header
  writeLines("CLUSTAL W (custom format)\n", con)
  
  for (start in seq(1, seqlen, by = block_size)) {
    end <- min(start + block_size - 1, seqlen)
    
    # Extract block
    block <- sapply(seqs, function(s) substr(s, start, end))
    
    # Write each sequence line
    for (i in seq_along(block)) {
      seq_line <- block[i]
      nongap_count <- nchar(gsub("-", "", seq_line))
      cum_pos[i] <- cum_pos[i] + nongap_count
      line <- sprintf(paste0("%-", name_width, "s%s %d"), names[i], seq_line, cum_pos[i])
      writeLines(line, con)
    }
    
    # Identity line
    #chars <- do.call(rbind, strsplit(block, "")) - hits a 10000 char limit as it names the columns!
    chars <- matrix(unlist(strsplit(block, "")), nrow = length(block), byrow = TRUE)
    identity <- apply(chars, 2, function(col) if (length(unique(col)) == 1) "*" else " ")
    identity_line <- sprintf(paste0("%-", name_width, "s%s"), "", paste(identity, collapse = ""))
    writeLines(identity_line, con)
    
    writeLines("", con)  # Blank line between blocks
  }
}

generate_variant_table <- function(ref, seqs, orf_start=1, ORF_lengths, type = "nuc") { #ORF_lengths only for nuc
  
  report_table <- data.frame(name = character(), variants = character(), summary = character(), stringsAsFactors = FALSE)
  if (type == "nuc") {
    alignment_orf_start <- which(test_seq != "-")[orf_start] # it is an alignment this is true for all
    ref_orf <- ref_seq[alignment_orf_start:(alignment_orf_start+3*ORF_lengths[1])]
  }
  
  for (name in names(seqs)) {
    if (name == names(ref)) next
    
    ref_seq <- unlist(strsplit(as.character(ref[[1]]), ""))
    test_seq <- unlist(strsplit(as.character(seqs[[name]]), ""))
    
    variants <- c()
    i <- 1
    pos <- 1  # 1-based position for c. or p. notation
    
    while (i <= length(ref_seq) && i <= length(test_seq)) {
      # Match
      if (ref_seq[i] == test_seq[i]) {
        i <- i + 1
        pos <- pos + 1
        next
      }
      
      # --- Deletion ---
      if (test_seq[i] == "-") {
        del_start <- pos
        while (i <= length(test_seq) && test_seq[i] == "-") {
          i <- i + 1
          pos <- pos + 1
        }
        if ((pos - del_start) == 1) {
          variants <- c(variants, sprintf("%s%ddel (del_start: %d)", 
                                          ifelse(type == "nuc", "c.", "p."), 
                                          del_start-orf_start+1,
                                          del_start))
        } else {
          variants <- c(variants, sprintf("%s%d_%ddel (del_start: %d)", 
                                          ifelse(type == "nuc", "c.", "p."), 
                                          del_start-orf_start+1, 
                                          pos-orf_start+1 - 1,
                                          del_start))
        }
        next
      }
      
      # --- Insertion ---
      if (ref_seq[i] == "-") {
        ins_start <- pos - 1
        inserted <- ""
        while (i <= length(ref_seq) && ref_seq[i] == "-") {
          inserted <- paste0(inserted, test_seq[i])
          i <- i + 1
        }
        variants <- c(variants, sprintf("%s%d_%dins%s (ins_start: %d)", 
                                        ifelse(type == "nuc", "c.", "p."), 
                                        ins_start-orf_start+1, 
                                        ins_start-orf_start+1+ 1, 
                                        inserted,
                                        ins_start))
        next
      }
      
      # --- Mismatch ---
      variants <- c(variants, sprintf("%s%d%s>%s (pos: %d)", 
                                      ifelse(type == "nuc", "c.", "p."), 
                                      pos-orf_start+1, 
                                      ref_seq[i], 
                                      test_seq[i],
                                      pos)) #)
      i <- i + 1
      pos <- pos + 1
    }
    
    # --- Handle trailing insertions or deletions ---
    while (i <= length(ref_seq)) {
      if (ref_seq[i] != "-") {
        variants <- c(variants, sprintf("%s%ddel (pos: %d)", 
                                        ifelse(type == "nuc", "c.", "p."), 
                                        pos-orf_start+1,
                                        pos))
        pos <- pos + 1
      }
      i <- i + 1
    }
    while (i <= length(test_seq)) {
      if (test_seq[i] != "-") {
        variants <- c(variants, sprintf("%s%d_%dins%s (pos-1: %d)", 
                                        ifelse(type == "nuc", "c.", "p."), 
                                        pos -orf_start+1- 1, 
                                        pos-orf_start+1, 
                                        test_seq[i],
                                        pos-1))
      }
      i <- i + 1
    }
    
    # --- Finalize table ---
    if (length(variants) == 0) {
      report_table <- rbind(report_table, data.frame(SequenceID = name, Variant = "ORF OK"))
    } else {
      report_table <- rbind(report_table, data.frame(SequenceID = name, Variant = paste(variants[1:min(100, length(variants))],
                                                                                        collapse = ",")))
    }
    
    #  ADDITIONAL LINE: TRUNCATION (for aa mode)
    if (type == "aa") {
      stop_pos <- which(test_seq == "*")[1]
      #stop_pos <- regexpr("\\*", test_seq)[1]
      
      
      ref_seq_no_gaps <- gsub("-", "", paste(ref_seq, collapse = ""))
      ref_length <- nchar(ref_seq_no_gaps)
      #Fails if some entries longer: ref_length <- length(ref_seq) #gsub("-", "", ref_seq))  # assuming ref_seq is also a string
      
      if (stop_pos > 0 && stop_pos < ref_length) {
        report_table <- rbind(report_table, data.frame(SequenceID = name, Variant = paste("TRUNCATION at residue", stop_pos-orf_start+1, "of", ref_length-orf_start)))
      }
      
    }
    
    #  ADDITIONAL LINE: FRAME SHIFT (in ORFfor nuc mode)
    if (type == "nuc") {
      query_orf <- test_seq[alignment_orf_start:(alignment_orf_start+3*ORF_lengths[name])]
      gap_run <- 0
      for (pos in seq_along(query_orf)) {
        if (query_orf[pos] == "-" && ref_orf[pos] != "-") { #mismatch to ref
          gap_run <- gap_run + 1
        } else {
          if (gap_run > 0) {
            if (gap_run %% 3 != 0) {
              #message(paste("Frame shift at ORF position", pos - gap_run))
              report_table <- rbind(report_table, data.frame(SequenceID = name, Variant = paste("FRAME SHIFT at ORF nt",  pos - gap_run)))
              # Optionally: store in report_table
            }
            
            gap_run <- 0
          }
        }
      }
    }
    
    
    
    # r <- rle(test_seq[c(alignment_orf_start:(alignment_orf_start+3*ORF_lengths[name]))]) #length(test_seq))]) # check only the ORF not the rest
    # pos <- cumsum(r$lengths) - r$lengths + 1  # Start positions of each run
    # 
    # for (j in seq_along(r$values)) {
    #   if (r$values[j] == "-") {
    #     len <- r$lengths[j]
    #     if (len %% 3 != 0) {
    #       shift_pos <- pos[j]
    #       report_table <- rbind(report_table, data.frame(SequenceID = name, Variant = paste("FRAME SHIFT at ORF nt", shift_pos)))
    #     }
    #   }
    # }
    
    #indels <- gregexpr("[-]+", test_seq)[[1]]
    # if (indels[1] != -1) {
    #   for (j in seq_along(indels)) {
    #     len <- attr(indels, "match.length")[j]
    #     if (len %% 3 != 0) {
    #       shift_pos <- indels[j]
    #       report_table <- rbind(report_table, data.frame(SequenceID = name, Variant = paste("FRAME SHIFT at nt", shift_pos)))
    #       break
    #     }}}
  }
  
  return(report_table)
}



# End of Functions
########################

#*
#*
#*
#*
#*
#*
#*
#*
#*
#*
#*
#*
#*
#*
#*
#*
#*


# Load csv  file (fsv not robust if users generate in notepad)
seq_list<- read.csv(file.path(input_path, SampleListName), header = TRUE, stringsAsFactors = FALSE)
#tsv <- read.table(paste0(input_path,"SeqToAnalyse_Nluc.txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
#colnames(seq_list) <- c("SequenceID", "Path", "Filename") - expects tehse

# Apply trimws only to character columns and keep as df
seq_list <- as.data.frame(lapply(seq_list, function(x) if (is.character(x)) trimws(x) else x))



# Read the reference sequence and clean out numbers spaces and \r\n


# Load reference sequence
# 1. Read raw template
raw_lines <- readLines(file.path(seq_list$Path[1],seq_list$Filename[1]))
# 2. Remove header and clean sequence
clean_seq <- gsub("[^ACGTacgt]", "", paste(raw_lines[-1], collapse = ""))
ref_seq <- DNAString(toupper(clean_seq))
#ref_seq <- readDNAStringSet("reference_cleaned.fasta")[[1]]



# Initialize alignments
nuc_alignments <- list()
aa_alignments <- list()
nuc_alignments[[seq_list$SequenceID[1]]] <- ref_seq
aa_alignments[[seq_list$SequenceID[1]]] <- extract_orf(seq=ref_seq, start=orf_start, id=seq_list$SequenceID[1])


# Process each plasmid
for (i in 2:nrow(seq_list)) {
  file_path <- file.path(seq_list$Path[i], seq_list$Filename[i])
  seq <- readDNAStringSet(file_path, format="fasta")[[1]]
  seq_rot <- rotate_to_match(seq, ref_seq, window_start = orf_start) # get match around start site otherwise the sample read will be out with an indel
  if (seq_list$MacrogenCorrection[i]=="yes")
    {
    # Macrogen_insert_pos is relative to SynGAP read, so correct for plasmid position and tag length (A1 start-compensated) 
    insert_pos <- Macrogen_insert_pos + orf_start-1 + tag_offset*3
    # Insert base
    seq_rot <- DNAString(paste0(subseq(seq_rot, 1, insert_pos - 1), Macrogen_base_to_insert, 
                                subseq(seq_rot, insert_pos, length(seq_rot))))
  }
  
  nuc_alignments[[seq_list$SequenceID[i]]] <- seq_rot
    # Save rotated sequence to individual FASTA file
    rotated_file <- file.path(input_path, paste0(seq_list$SequenceID[i], "_rotated.fasta"))
    writeXStringSet(DNAStringSet(seq_rot), filepath = rotated_file)
    
  aa_alignments[[seq_list$SequenceID[i]]] <- extract_orf(seq=seq_rot, start=orf_start,id=seq_list$SequenceID[i])


}

# THIS SHOULD BE A SINGLE RUN NOT IN THE LOOP!
# Write nucleotide alignment
# Convert list to DNAStringSet
all_seqs <- DNAStringSet(nuc_alignments)
# Perform multiple sequence alignment

alignment <- msa(all_seqs, method = "ClustalOmega", order= "input")  # or "Muscle", "ClustalW"
# Export aligned sequences
# Assuming `alignment` is the result of msa()
aligned_set <- as(alignment, "DNAStringSet")#[match(seq_list$SequenceID, names(as(alignment, "DNAStringSet")))]
# Convert to seqinr::alignment format
alignment_seqinr <- msaConvert(alignment, type = "seqinr::alignment")
# Write to Clustal format - first does not have numbers or identity line
#write_clustal(alignment_seqinr, file = paste0(input_path, "clustal_alignment.aln"))
align_filename <-file.path(input_path, "custom_clustal_alignment.aln")
write_custom_alignment(alignment=alignment_seqinr, file = align_filename)
custom_clustal_alignment <- readLines(align_filename)# read back to add to pdw

# Write amino acid alignment
# Assuming `aa_alignments` is a list or vector of amino acid sequences
all_aa_seqs <- AAStringSet(aa_alignments)


# Perform multiple sequence alignment and align according to the original
aa_alignment <- msa(all_aa_seqs, method = "ClustalOmega",order = "input")  # or "Muscle", "ClustalW"
# Convert to AAStringSet and make sure the order is as the input
aligned_aa_set <- as(aa_alignment, "AAStringSet") # if need to reorder [match(seq_list$SequenceID, names(as(aa_alignment, "AAStringSet")))]
# If still need to reorder,  aa_alignment <- msa(aligned_aa_set, method = "ClustalOmega", order = "input")
# Convert to seqinr::alignment format
alignment_seqinr_aa <- msaConvert(aa_alignment, type = "seqinr::alignment")
# Write to Clustal format - firs does not have numbers or identity line
#write_clustal(alignment_seqinr_aa, file = paste0(input_path, "clustal_aa_alignment.aln"))
align_aa_filename <-paste0(input_path, "custom_clustal_alignment_aa.aln")
write_custom_alignment(alignment_seqinr_aa, file = align_aa_filename)
custom_clustal_alignment_aa <- readLines(align_aa_filename) # read back to add to pdw

#Summarise differences
ORF_lengths=sapply(aa_alignments,length)
nuc_table <- generate_variant_table(ref=aligned_set[1], seqs=aligned_set, orf_start=orf_start+tag_offset*3,ORF_lengths,type="nuc")
#nuc_table <- generate_variant_table(ref=nuc_alignments[1], seqs=nuc_alignments, "nuc")
aa_table <- generate_variant_table(ref=aligned_aa_set[1], seqs=aligned_aa_set, orf_start=tag_offset, ORF_lengths,type="aa")


write.table(nuc_table, file=(file.path(input_path, "nucleotide_variants.tsv")), sep="\t", row.names=FALSE, quote=FALSE)
write.table(aa_table, file=(file.path(input_path, "amino_acid_variants.tsv")), sep="\t", row.names=FALSE, quote=FALSE)


# Save pdw format

#Load pdw template
template_file <-file.path(input_path, "Template.pdt")
template_lines <- readLines(template_file)# <- readChar(template_file, nchars = file.info(template_file)$size)
#convert the original seq file list
lines <- apply(seq_list, 1, function(row) paste(row, collapse = ","))
header_line <- paste(colnames(seq_list), collapse = ",")

for (i in 2:nrow(seq_list)) {
sample_name<-seq_list$Filename[i]
p_variant_type <- aa_table[aa_table$SequenceID==seq_list$SequenceID[i],2][1]
c_variant_type <- nuc_table[nuc_table$SequenceID==seq_list$SequenceID[i],2]

pdw<-c(template_lines[1:2], paste0(template_lines[3],
              sample_name, ":", p_variant_type),
       template_lines[-c(1:3, length(template_lines))])
# indicate MacrogenCorrection applied
if (seq_list$MacrogenCorrection[i]=="yes")
  pdw<-c(pdw, paste0("HEADER ", # every header line starts with this 
                     c("", "Macrogen correction applied - adding one G to the a run of 6G (read as 5) at SynGAP base 330",
                    "NOTE - confirm this with Sanger Sequencing!",   
                     "")))
# sample name and result
pdw<-c(pdw, paste0("HEADER ", # every header line starts with this 
              c("", paste0("Sample name: ", sample_name),
              paste0("protein variant type: [", p_variant_type,"]"), 
              paste0("nucleotide variant type: [", c_variant_type,"]"), 
              paste0("(nuc match aligned to target cds, for location in plasmid add", orf_start+tag_offset*3," bases)")
              )))


# add the alignments 
pdw<-c(pdw, rep("HEADER ", 5),paste0("HEADER amino acid alignment"),paste0("HEADER ", custom_clustal_alignment_aa), "HEADER ")
pdw<-c(pdw, rep("HEADER ", 5),paste0("HEADER nucleotide alignment"),paste0("HEADER ", custom_clustal_alignment), "HEADER ")

# add the original file defining what exactly was aligned
pdw<-c(pdw, paste0("HEADER ", # every header line starts with this 
                   c(rep("", 5), 
                    "List of Sequences Analysed","",
                    header_line,
                    lines)
))
#add the nucleotide sequence
pdw<-c(pdw, "Sequence ..",as.character( nuc_alignments[[seq_list$SequenceID[i]]]))
# write to the analysis folder AND to the original folder
 writeLines(pdw, file.path(input_path, paste0(seq_list$SequenceID[i],"_",seq_list$Filename[i], "_rotated.pdw")))
 writeLines(pdw, file.path(seq_list$Path[i], paste0(seq_list$SequenceID[i],"_",seq_list$Filename[i], "_autorotated.pdw")))
 

}





