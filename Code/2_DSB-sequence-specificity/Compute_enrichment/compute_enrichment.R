#autor: paulina dziubanska-kusibab, contribution: hilmar berger
library(Biostrings)
library(stringr)

####################################################################################################
compute_enrichment <- function(testfile, controlfile, oligo_size, outputname, output_folder) {
  sample <- readDNAStringSet(testfile)
  background <- readDNAStringSet(controlfile)
  
  oligo_size <- (oligo_size)
  #matrix: length(input) x 4^oligo_size
  ok_matrix_size <- 5 * 10^5 * 4^6
  
  Cnt_oligo <- function(input)
  {
    
    len <- length(input)
    part_width <- floor((ok_matrix_size)/(4^oligo_size))
    number_of_parts <- ceiling(len / part_width)
    
    oligo_counts <- rep(0,4^oligo_size)
    
    print(paste0(deparse(substitute(input)),', length: ', len))
    for (i in 1:number_of_parts)
    {
      start <- 1+(i-1)*part_width
      end <- min(len, i*part_width)
      print(paste0('start: ', start,', end: ',end))
      # make sure that each motif is counted at most once per input sequence
      oligo_counts <- oligo_counts + apply( oligonucleotideFrequency(input[start:end],oligo_size)>0 ,2,sum)
    }
    return(oligo_counts)
  }
  
  # oligo counts define the number of input sequences containing each oligo
  cnt_sample <- Cnt_oligo(sample)
  cnt_background <- Cnt_oligo(background)
  # ratios define the proportion of input sequences containing each oligo
  ratio_sample <- cnt_sample/length(sample)
  ratio_background <- cnt_background/length(background)
  ratio <- ratio_sample/ratio_background
  log_ratio <- log2(ratio)
  scale <- ifelse(log_ratio <0, "below", "above")
  nt_content <- data.frame('Nt' = names(cnt_sample), 'count_sample' = cnt_sample, 'ratio_sample' = ratio_sample,
                           'count_background' = cnt_background, 'ratio_background' = ratio_background, 'ratio' = ratio, 
                           size_sample = length(sample), size_background = length(background), 
                           'logRatio' = log_ratio, 'scale' = scale)
  
  save(nt_content, file=file.path(output_folder, paste0(outputname,"_",oligo_size,"nt.RData")) )
  write.table(nt_content, file.path(output_folder, paste0(outputname,"_",oligo_size,"nt.txt")), sep="\t", dec=".", quote=F, col.names=NA)
  
}
