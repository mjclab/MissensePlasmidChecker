# Missense Plasmid Checker

Purpose

R script developed for checking the sequences and reporting on panels of plasmids encoding missense coding sequences preceded by a tag of a specified sequence. 
Adapted for projects on SynGAP1 and works for that. Should be broadly applicable, but current implementation may fail where assumptions are not met (e.g. search sequence should not overlap end of file).

See R script for details. 

Here is a summary

Input files
i) reference .fasta  
ii) 1 or more .fasta to compare against reference
iii) csv ascii text file 
- line 1 is header with 4 column titles (SequenceID,MacrogenCorrection,Filename,Path) describing the 4 items rewuired for each subsequent line: an ID, yes/no to apply a full-length sequencing error correction, the fasta file name including extension, the path to the fasta file
- line 2 is the reference
- subsequent lines list the samples files
iv) a Template.pdt file, this is a template to generate a pdraw sequence file (acaclone.com) but could be adapted to other viewers
(iii) and (iv) should be in the folder for a specific analysis series
(i) and (ii) can be anywhere as defined in (iii). I usually copy the reference into the analysis folder 

Output files
i)  rotated .fasta files to align the ORF start site
ii) a summary of cds and nucleotide differences in tsv format
ii) pdw plasmid viewer versions each containing, as a header
- the full cds and nucleotide alignments
- summary of cds and nucleotide differences
  

Other inputs required
at the top of the R script
i) the path to the csv ascii file
ii) the analysis folder path as an output path
iii) the offset in nucleotides to the ORF-of-interest in the reference, where the N-terminal tag starts of there is one
iv) the offset in codons (!) from the above to the start of the coding sequence (this is so the reported missense locations refer to the coding sequence, regardless of tag length. If there is a mutation in the tag a correspondingnegative position is reported)
v) the position (nucleotides) of the optional insert to correct for a misread in data produced by some full plasmid sequencing service providors
vi) the identity of the base that was lost at the above position and needs to be readded in data produced by some full plasmid sequencing service providors
(v)-(vi) are applied only if "yes" specified in the csv file second column. The correction has to be validated manually by a separate sanger sequencing run


Developed by: Michael J Courtney, Turku Bioscience Centre, University of Turku
Tested by: as above, add your name here 
Funding : SynGAP Research Fund, Research Council of Finland, EU-IMPULSE/EU Horizon Europe

R version requirements in renv.lock

