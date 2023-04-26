###############################
#         CMD Options         #
###############################

suppressPackageStartupMessages( library( optparse ) )
#suppressPackageStartupMessages( library( futile.logger ) )

option_list = list(
  make_option( c( "-c", "--counts" ), type="character", default='.',
               help="Count file", metavar="character" ),
  make_option( c( "-e", "--ext" ), type="character", default='.alleles_frequency_table.tsv',
               help="Count file extension [Default = .alleles_frequency_table.tsv]", metavar="character" ),
  make_option( c( "-a", "--amplicon" ), type="character",
               help="Amplicon sequence", metavar="character" ),
  make_option( c( "-o", "--out" ), type="character", default='.',
               help="Output directory [Default = .]", metavar="character" )
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

###############################
#     Option validation       #
###############################

#futile.logger::flog.info( "Validating inputs..." )

if ( is.null( opt$counts) ) {
#  flog.fatal( "Please provide a count file" )
  stop( "Please provide a count directory" )
}

if ( ! file.exists( opt$counts ) ) {
#  flog.fatal( "Count file does not exist." )
  stop( "Count file does not exist." )
}

if ( is.null( opt$amplicon ) ) {
#  flog.fatal( "Please provide an amplicon sequence" )
  stop( "Please provide an amplicon sequence" )
}

if ( !dir.exists( opt$out ) ) {
#  flog.fatal( "Output directory does not exist." )
  stop( "Output directory does not exist." )
}

countFile = opt$counts
countFileExtension = opt$ext
amplicon = opt$amplicon
outputDirectory = opt$out

###############################
#         Libraries           #
###############################

#futile.logger::flog.info( "Loading libraries..." )

suppressPackageStartupMessages( library( tidyverse ) )

###############################
#         Functions           #
###############################

get_cigar <- function( ref, alt ) {
  
  # Split string into character vector 
  alt.split <- str_split( alt, "") %>% unlist()
  ref.split <- str_split( ref, "") %>% unlist()
  
  aln.length <-  length( ref.split )
  
  if( length( alt.split) != aln.length ) {
    stop( "Reference and aligned strings are not the same length" )
  }
  
  # Get mapping codes at each position in the alignment
  # M = match/mismatch, D = deletion, I = indel 
  aln <- vector()
  for ( i in 1:aln.length ) {
    aln[i] = case_when(
      alt.split[i] == ref.split[i] ~ "M",
      alt.split[i] != ref.split[i] && alt.split[i] != '-' && ref.split[i] != '-' ~ "M",
      alt.split[i] == '-' && ref.split[i] != '-' ~ "D",
      ref.split[i] == '-' && alt.split[i] != '-' ~ "I",
      TRUE ~ 'X'
    )
  }
  
  # Convert alignment codes into cigar string
  current.code = ''
  cigar = ''
  count = 0
  for ( i in 1:aln.length ) {
    if ( current.code == '' ) {
      current.code = aln[i]
      count = 1
    } else {
      if ( current.code == aln[i] ) {
        count <- count + 1
      } else {
        cigar <- paste0( cigar, count, current.code )
        current.code = aln[i]
        count = 1
      }
    }
    if ( i == aln.length ) {
      cigar <- paste0( cigar, count, current.code )
    }
  }
  
  return( cigar )
}


###############################
#    Get count filenames      #
###############################

#futile.logger::flog.info( "Getting alignments..." )

prefix <- str_replace( basename( countFile ), countFileExtension, '' )

crispresso <- read.delim( file.path( countFile ), 
                          header = T, stringsAsFactors = FALSE,
                          col.names = c( 'aln_seq', 'ref_seq', 'ref_name', 'status', 'del', 'ins', 'mut', 'n_reads', 'pct_reads' ) )

#futile.logger::flog.info( "Building CIGAR strings..." )

crispresso.cigar <- crispresso %>% 
                        rowwise() %>% 
                        mutate( ref_start = str_locate( ref_seq, amplicon )[1],
                                ref_end = str_locate( ref_seq, amplicon )[2],
                                mod_ref_seq = ifelse( str_detect( ref_seq, amplicon ), str_sub( ref_seq, ref_start, ref_end), ref_seq ),
                                mod_aln_seq = ifelse( str_detect( ref_seq, amplicon ), str_sub( aln_seq, ref_start, ref_end ), aln_seq ), 
                                cigar = get_cigar( mod_ref_seq, mod_aln_seq ) ) 

#futile.logger::flog.info( "Writing alignments..." )

write.table( crispresso.cigar, file.path( outputDirectory, paste0( prefix, '.cigar_alignments.tsv' ) ), row.names = F, sep = "\t", quote = F )

#futile.logger::flog.info( "Getting CIGAR counts..." )

crispresso.cigar.counts <- crispresso.cigar %>%
                            group_by( cigar ) %>% 
                            summarise( n = sum( n_reads ), pct = sum( pct_reads ) ) %>%
                            arrange( -pct ) 

#futile.logger::flog.info( "Writing CIGAR counts..." )  
write.table( crispresso.cigar.counts, file.path( outputDirectory, paste0( prefix, '.cigar_counts.tsv' ) ), row.names = F, sep = "\t", quote = F )
