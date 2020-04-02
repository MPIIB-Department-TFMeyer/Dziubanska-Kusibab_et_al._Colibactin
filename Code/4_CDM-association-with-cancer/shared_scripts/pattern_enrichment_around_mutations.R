require(GenomicRanges)
require(data.table)

# Compute enrichment of all patterns in pattern_positions among mutations in mut_ranges (filtered to those contained in target regions) 
# pattern_postions: GRanges containing genomic positions of selected kmer patterns (should be on both strands). Must contain a metadata column named "pattern" which describes the pattern sequence
# mutation ranges: GRanges containing genomic positions of called variants. Must contain a column defining the entity (primary grouping of cases) to which the sample belongs a an sample ID column. Mutations should be unique within each entity. Each variant (i.e. unique position + nucleotide change ) must be labeled with a unique variant ID.
# : size (bp) of the genomic regions considered for pattern overlap (e.g. exon ranges or all chromosomes)
# oligo_freqs_in_target: named integer vector containing the frequency of each pattern within target regions
# pattern_bp_covered: named integer vector containing the number of bp covered by each pattern in pattern_positions within  target regions
# entity_column: column in metadata of mutation_ranges object containing the subclass each sample/mutation belongs to (e.g. cancer type/site)
# sample_column: column in metadata of mutation_ranges object containing the sample ID
# variant_id_column: column in metadata of mutation_ranges object containing the unique (chromosome, position, nucĺeotide change) variant ID
compute_pattern_enrichment_targeted <- function(pattern_positions, mutation_ranges, target_size, oligo_freqs_in_target, pattern_bp_covered, entity_column="entity", sample_id_column="sampleID", variant_id_column = "Mutation.ID") {
  
  # All ranges covered by the seleced pattern on any of both strands
  #tmp = reduce(pattern_positions[pattern_positions$index==pi], ignore.strand=T)
  # identify pattern covered positions that overlap mutated sites
  po = findOverlaps(mutation_ranges, pattern_positions, ignore.strand=T)
  
  # mutations overlapping the same pattern on both + and - strand should only be counted once
  tmp = data.table(mutID = eval(parse(text=paste0("mutation_ranges[queryHits(po)]$",variant_id_column))),
                   entity = eval(parse(text=paste0("mutation_ranges[queryHits(po)]$",entity_column))), 
                   pattern =as.character(pattern_positions[subjectHits(po)]$pattern), 
                   sampleID = eval(parse(text=paste0("mutation_ranges[queryHits(po)]$",sample_id_column))),
                   stringsAsFactors = F)
  
  # count unique mutation overlaps for each pattern
  #tmp2 = with(unique(tmp[, c("mutID", "entity", "p")]), tapply(entity, p, length))
  
  
  by_entity_mutation_counts = tapply(eval(parse(text=paste0("mutation_ranges[queryHits(po)]$",variant_id_column))), eval(parse(text=paste0("mutation_ranges[queryHits(po)]$",entity_column))), function(x) length(unique(x)))

  by_entity_sample_counts = tapply(eval(parse(text=paste0("mutation_ranges[queryHits(po)]$",sample_id_column))), eval(parse(text=paste0("mutation_ranges[queryHits(po)]$",entity_column))), function(x) length(unique(x)))
  
  # mutation hits is number of unique mutation sites overlapped by a given pattern
  #d = data.frame(pattern = names(tmp2), mutation_hit = tmp2, stringsAsFactors = F)
  d = unique(tmp[, c("mutID", "entity", "pattern")])[, .(mutation_hit=.N), by=c("pattern", "entity")]
  rm(tmp)
  # total count of pattern sites in target regions
  d$total_pattern_sites = oligo_freqs_in_target[d$pattern]
  # total bases covered by given pattern in target regions
  d$pattern_bp = pattern_bp_covered[d$pattern]
  # Proportion of total pattern sites overlapping a mutation
  d$hit_prop = d$mutation_hit / d$total_pattern_sites
  d$total_mutation_counts = by_entity_mutation_counts[as.character(d$entity)]
  d$sample_counts = by_entity_sample_counts[as.character(d$entity)]
  # Number of mutations not overlapping a given pattern. Those necessarily overlap all other patterns.
  d$mutation_no_hit = d$total_mutation_counts-d$mutation_hit
  # total bases not covered by given pattern as remainder from exonic ranges minus pattern covered sequence
  d$non_pattern_bp = target_size-d$pattern_bp
  # Proportions of bp covered by pattern or not affected by mutations. 
  # Note that this assumes each mutation covers a single bp - this could be improved if necessary
  d$hit_per_bp_in_pattern = d$mutation_hit/d$pattern_bp
  d$hit_per_bp_in_non_pattern = d$mutation_no_hit/d$non_pattern_bp
  # ratio of proportions wiith vs. without pattern
  d$ratio_pattern_vs_rest = d$hit_per_bp_in_pattern/d$hit_per_bp_in_non_pattern
  
  ds = d[order(-d$ratio_pattern_vs_rest),]
  
  return(ds)
}

