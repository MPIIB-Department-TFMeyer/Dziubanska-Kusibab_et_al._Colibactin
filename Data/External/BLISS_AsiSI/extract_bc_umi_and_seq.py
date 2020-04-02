import sys, random, os, re, gzip
import HTSeq
from string import maketrans
from operator import itemgetter
#from fuzzywuzzy import fuzz
#from fuzzywuzzy import process

nuc_trans = maketrans("ATCG", "TAGC")

def extract_bc_umi_seq(input_files, barcodes):

    seq_ofile_handles = {}
    index_ofile_handles = {}
    for bc in barcodes.keys():
        sample = barcodes[bc]
        seq_ofile_handles[bc] = gzip.GzipFile("%s_genomic.fastq.gz" % sample, "w")
        index_ofile_handles[bc] = gzip.GzipFile("%s_BC_UMI.fastq.gz" % sample,"w")

    excluded_ofile = gzip.GzipFile("undetermined.fastq.gz", "w")
        
    for ifile_name in input_files:
        #print ifile_name

        if ifile_name.endswith(".gz"):
            ifile=gzip.GzipFile(ifile_name,"r")
        else:
            ifile = ifile_name
        
        total_read_cnt = 0
        accepted_read_cnt = 0
        all_bcs = {}
        for r in HTSeq.FastqReader(ifile ):
            total_read_cnt += 1
            #umi = r.seq[0:8]
            bc = r.seq[8:16]
            
            if bc in all_bcs:
                all_bcs[bc] += 1
            else:
                all_bcs[bc] = 1

            if bc in barcodes:
                accepted_read_cnt += 1
                new_genomic_item = HTSeq.SequenceWithQualities(r.seq[16:(len(r.seq)+1)], r.name, r.qualstr[16:(len(r.seq)+1)] )
                new_index_item = HTSeq.SequenceWithQualities(r.seq[0:16], r.name, r.qualstr[0:16] )

                new_genomic_item.write_to_fastq_file(seq_ofile_handles[bc])
                new_index_item.write_to_fastq_file(index_ofile_handles[bc])
            else:
                r.write_to_fastq_file(excluded_ofile)

    print "\nTotal reads seen: %d\tAssigned to known barcodes: %d (%g %%)\tUndetermined: %d (%g %%)" % (total_read_cnt, accepted_read_cnt, round(float(accepted_read_cnt)/total_read_cnt * 100,2), total_read_cnt - accepted_read_cnt, round(float(total_read_cnt - accepted_read_cnt)/total_read_cnt * 100,2) )

    print "\nTop 100 observed barcodes with frequencies:"
    for k,v in sorted(all_bcs.items(), key=itemgetter(1), reverse=True)[0:100]:
        found = "?" if not k in barcodes else barcodes[k]
        print "%s\t%d\t%s" % (k, v, found )
            
    for bc in barcodes.keys():
        seq_ofile_handles[bc].close()
        index_ofile_handles[bc].close()

    excluded_ofile.close()

###############################################################################    
if __name__ == '__main__':

    import argparse
    import time
    parser = argparse.ArgumentParser(description='Extract barcopdes, UMI and sequences for alignment from FASTQ. Split output FASTQ and indices by BC')
    parser.add_argument('fastq_r1', type=str, nargs='+', help='FASTQ read 1')

    parser.add_argument('-b','--barcode-table', action='store', type=str, default=None,
                       help='File name of tab-separated text files with sample name in first column and barcode in second column. First row will be ignored.')

    args = parser.parse_args()

    input_files = args.fastq_r1
    print "Input Files:"
    for ff in input_files:
        print ff
    
    barcodes = {}
    if args.barcode_table is None:
        sys.exit("No barcode to sample translation has been provided")
    else:
        with open(args.barcode_table, "r") as bcfile:
            linecnt = 0
            for line in bcfile.readlines():
                linecnt += 1
                if linecnt == 1:
                    continue
        
                fields = line.strip("\n\r").split("\t", -1)
                barcodes[fields[1]] = fields[0]

    print "\nKnown barcodes:"
    for k in barcodes.keys():
        print "%s\t%s" % (k, barcodes[k])
        
    extract_bc_umi_seq(input_files, barcodes)
