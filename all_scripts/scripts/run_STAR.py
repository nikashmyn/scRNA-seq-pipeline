#!/usr/bin/env python
# Original Author: Francois Aguet
# Modifications: Etai Jacob
# Last update: Dec. 13 2017

#When running STARlong - make sure the indexed file was made using STARlong and not STAR

#Not good for splice variants

import argparse
import os, sys
import subprocess
import gzip
import shutil
from datetime import datetime

parser = argparse.ArgumentParser(description='Run STAR')
parser.add_argument('index', help='e.g. /xchip/gtex/resources/STAR_genomes/STAR_genome_v19')
parser.add_argument('annotation_gtf', help='e.g. /xchip/gtex/resources/gencode.v19.transcripts.patched_contigs.gtf')
parser.add_argument('fastq_fnames', help='Paths to fastq or fastq.gz files (at most 2 for paired ends separeted by ,')
parser.add_argument('prefix', help='Prefix for output file names')
parser.add_argument('--twopassMode', default='Basic', help='1 or 2-pass mapping mode.')
parser.add_argument('--isPairedEnd', help='If one file is provided with R1, with Illumina R1/R2 naming, then R2 will be substituted')
parser.add_argument('--hasMultipleFastqs', default="0", help='If there is more than a single fastq file (paired or not)')
parser.add_argument('--outFilterMultimapNmax', default='20')
parser.add_argument('--alignSJoverhangMin', default='8')
parser.add_argument('--alignSJDBoverhangMin', default='2')
parser.add_argument('--outFilterMismatchNmax', default='999')
parser.add_argument('--outFilterMismatchNoverLmax', default='0.1')
parser.add_argument('--alignIntronMin', default='20')
parser.add_argument('--alignIntronMax', default='1000000')
parser.add_argument('--alignMatesGapMax', default='1000000')
parser.add_argument('--outFilterType', default='BySJout')
parser.add_argument('--outFilterScoreMinOverLread', default='0.33')
parser.add_argument('--outFilterMatchNminOverLread', default='0.33')
parser.add_argument('--limitSjdbInsertNsj', default='1200000')
parser.add_argument('--outSAMstrandField', default='intronMotif')
parser.add_argument('--outFilterIntronMotifs', default='None', help="Use 'RemoveNoncanonical' for Cufflinks compatibility")
parser.add_argument('--alignSoftClipAtReferenceEnds', default='No')
parser.add_argument('--quantMode', default='TranscriptomeSAM GeneCounts', help='Outputs read counts, and a BAM with reads in transcriptome coordinates')
parser.add_argument('--outSAMtype', default='BAM SortedByCoordinate')
parser.add_argument('--outSAMunmapped', default='Within', help='Keep unmapped reads in output BAM')
parser.add_argument('--outSAMattrRGline', default='ID:rg1 SM:sm1', help='Adds read group line to BAM header; required by GATK')
parser.add_argument('--outSAMattributes', default='NH HI AS nM NM')
parser.add_argument('--chimSegmentMin', default='15', help='Minimum chimeric segment length; switches on detection of chimeric (fusion) alignments')
parser.add_argument('--chimJunctionOverhangMin', default='15', help='Minimum overhang for a chimeric junction')
parser.add_argument('--genomeLoad', default='NoSharedMemory')
parser.add_argument('--limitBAMsortRAM', default=str(1024**3*40*4))  # 95GB required for some long read samples
parser.add_argument('-t', '--threads', default='32', help='Number of threads')
parser.add_argument('--reSortBAM', dest='reSortBAM', action='store_true', help='Sort transcriptome BAM by read name and into adjacent pairs')
parser.set_defaults(reSortBAM=False)
args = parser.parse_args()

# ENCODE options (from STAR manual)
# --outFilterType BySJout
# --outFilterMultimapNmax 20 # same in ICGC
# --alignSJoverhangMin 8 (STAR default 5, not set in ICGC)
# --alignSJDBoverhangMin 1 # same in ICGC (STAR default 3)
# --outFilterMismatchNmax 999 # 10 in ICGC -> problem for long reads?
# --outFilterMismatchNoverLmax 0.04 (0.3 STAR default, unchanged ICGC)
# --alignIntronMin 20
# --alignIntronMax 1000000 # 500000 in ICGC
# --alignMatesGapMax 1000000 # same in ICGC
# ICGC defaults:
# --outFilterScoreMinOverLread 0.33
# --outFilterMatchNminOverLread 0.33

# --outFilterMismatchNoverLmax 0.04 corresponds to: 0.04*152=6; 0.04*500=20 (STAR default: 0.3)
# --outFilterMismatchNoverLmax 0.1 corresponds to: 0.1*152=15; 0.1*500=50 (STAR default: 0.3)

# Additional options:
#  --outSAMstrandField intronMotif: required for Cufflinks (for unstranded data)
#  --outFilterIntronMotifs RemoveNoncanonical (recommended for Cufflinks)

# ICGC: --outSAMattributes", default=["NH", "HI", "NM", "MD", "AS", "XS"]
# STAR default :NH HI AS nM


# with open(args.fastq_list) as fqlist:
#     fastq1 = fqlist.readline().strip()
#     fastq2 = fqlist.readline().strip()
if args.hasMultipleFastqs == "0":
  print "Going for single fastq mode.."
  fnames = args.fastq_fnames
  fastq1 = fnames
  fastq2 = ""
  if args.isPairedEnd == "1":
    fastq2 = fastq1
    fastq2 = fastq1.replace("_R1_", "_R2_")
    if(id(fastq2) == id(fastq1)):
      print "No change (_R1_) - going for next format."
      fastq2 = fastq1.replace(".1.fastq", ".2.fastq")
      if(id(fastq2) == id(fastq1)):
        print "No change (.1.fastq) - going for next format."
        fastq2 = fastq1.replace(".R1.fastq", ".R2.fastq") 
        if(id(fastq2) == id(fastq1)):
          print "No change (.R1.fastq) - going for next format."
          fastq2 = fastq1.replace("_1.fastq", "_2.fastq")
          if(id(fastq2) == id(fastq1)):
            print "No change (_1.fastq) - exiting program since no other format exists."
            sys.exit()
  
  tmpFastq = args.fastq_fnames.split(",")[0]
  if tmpFastq.endswith(".gz"): 
    with gzip.open(tmpFastq) as f:
      f.readline()
      seq = f.readline().strip()
      read_length = len(seq)
  else:
    with open(tmpFastq) as f:
      f.readline()
      seq = f.readline().strip()
      read_length = len(seq)  
    
else:
  print "Going for multiple fastq mode.."
  fnames = args.fastq_fnames.split(",")
  fastq1 = fnames[0]
  fastq2 = ""
  if len(fnames) > 1:
    fastq2 = fnames[1]
  elif args.isPairedEnd == "1":
    fastq2 = fastq1.replace("_R1_", "_R2_")
    if(id(fastq2) == id(fastq1)):
      print "No change (_R1_) - going for next format."
      fastq2 = fastq1.replace(".1.fastq", ".2.fastq")
      if(id(fastq2) == id(fastq1)):
        print "No change (.1.fastq) - going for next format."
        fastq2 = fastq1.replace(".R1.fastq", ".R2.fastq") 
        if(id(fastq2) == id(fastq1)):
          print "No change (.R1.fastq) - exiting program."
          sys.exit()
  
  if fastq1.endswith(".gz"): 
    with gzip.open(fastq1) as f:
      f.readline()
      seq = f.readline().strip()
      read_length = len(seq)
  else:
    with open(fastq1) as f:
      f.readline()
      seq = f.readline().strip()
      read_length = len(seq)
      

print "\nfastq1 is "+fastq1
print "fastq2 is "+fastq2
print "read length is ", read_length,"\n"

#sys.exit()



overhang = int((read_length + 5)/10.0)*10

if read_length<=100:
    starcmd='STAR'
else:
    starcmd='STAR' #'STARlong'

# set up command
cmd = starcmd+' --runMode alignReads --runThreadN '+args.threads
print "looking for "+args.index+'.readLength'+str(overhang-1)
if os.path.exists(args.index+'.readLength'+str(overhang-1)):
    cmd += ' --genomeDir '+args.index+'.readLength'+str(overhang-1)+ ' --sjdbOverhang '+str(overhang-1)
    print "Found index: "+ args.index+'.readLength'+str(overhang-1)
else:
    cmd += ' --genomeDir '+args.index +' --sjdbGTFfile '+args.annotation_gtf

cmd += ' --twopassMode '+args.twopassMode\
    +' --outFilterMultimapNmax '+args.outFilterMultimapNmax\
    +' --alignSJoverhangMin '+args.alignSJoverhangMin+' --alignSJDBoverhangMin '+args.alignSJDBoverhangMin\
    +' --outFilterMismatchNmax '+args.outFilterMismatchNmax+' --outFilterMismatchNoverLmax '+args.outFilterMismatchNoverLmax\
    +' --alignIntronMin '+args.alignIntronMin+' --alignIntronMax '+args.alignIntronMax+' --alignMatesGapMax '+args.alignMatesGapMax\
    +' --outFilterType '+args.outFilterType\
    +' --outFilterScoreMinOverLread '+args.outFilterScoreMinOverLread+' --outFilterMatchNminOverLread '+args.outFilterMatchNminOverLread\
    +' --limitSjdbInsertNsj '+args.limitSjdbInsertNsj\
    +' --readFilesIn '+fastq1+' '+fastq2
    #+' --sjdbOverhang '+str(overhang-1) \
if os.path.splitext(fastq1)[1]=='.gz':
    cmd += ' --readFilesCommand zcat'
cmd += ' --outFileNamePrefix '+args.prefix+'.'\
    +' --outSAMstrandField '+args.outSAMstrandField+' --outFilterIntronMotifs '+args.outFilterIntronMotifs\
    +' --alignSoftClipAtReferenceEnds '+args.alignSoftClipAtReferenceEnds+' --quantMode '+args.quantMode\
    +' --outSAMtype '+args.outSAMtype+' --limitBAMsortRAM '+args.limitBAMsortRAM\
    +' --outSAMunmapped '+args.outSAMunmapped+' --genomeLoad '+args.genomeLoad
if int(args.chimSegmentMin)>0:
    cmd += ' --chimSegmentMin '+args.chimSegmentMin+' --chimJunctionOverhangMin '+args.chimJunctionOverhangMin
cmd += ' --outSAMattributes '+args.outSAMattributes+' --outSAMattrRGline '+args.outSAMattrRGline

print cmd



# run STAR
r = subprocess.call(cmd, shell=True, executable='/bin/bash')
if r != 0:
    exit(1)

# set permissions
for r,d,f in os.walk(args.prefix+'._STARpass1'):
    os.chmod(r, 0o755)

# delete unneeded files
shutil.rmtree(args.prefix+'._STARgenome')
shutil.rmtree(args.prefix+'._STARtmp')

# index BAMs
print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Indexing BAM')
sys.stdout.flush()

cmd = 'samtools index '+args.prefix+'.Aligned.sortedByCoord.out.bam'
subprocess.call(cmd, shell=True, executable='/bin/bash')
print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finished indexing BAM')
sys.stdout.flush()

if 'TranscriptomeSAM' in args.quantMode:
    if not args.reSortBAM:
        # Transcriptome BAM is correctly sorted since STAR 2.4.1d; the following is no longer needed. Rename output for compatibility
        os.rename(args.prefix+'.Aligned.toTranscriptome.out.bam', args.prefix+'.Aligned.toTranscriptome.sorted.bam')
    else:
        # sort transcriptome BAM by read name and into adjacent pairs (for multi-mapped reads)
        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Sorting transcriptome BAM')
        sys.stdout.flush()
        cmd = 'samtools view '+args.prefix+'.Aligned.toTranscriptome.out.bam | sort -S 75G -s -k 1,1 -k 13,13 |\
            cat <(samtools view -H '+args.prefix+'.Aligned.toTranscriptome.out.bam) - | samtools view -@ '+args.threads+' -bS - > \
                '+args.prefix+'.Aligned.toTranscriptome.sorted.bam'
        # -T /broad/hptmp
        subprocess.call(cmd, shell=True, executable='/bin/bash')
        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finished sorting transcriptome BAM')
        sys.stdout.flush()

        # remove unsorted BAM
        os.remove(args.prefix+'.Aligned.toTranscriptome.out.bam')


if int(args.chimSegmentMin)>0:
#in cmd made the " - " a " -o " in the sort cmd changed "...sorted" to "...sorted.bam" | could have just made it temp with -T and kept the prefix argument "...sorted"
    cmd = 'samtools view -bS '+args.prefix+'.Chimeric.out.sam | samtools sort -o '+args.prefix+'.Chimeric.out.sorted.bam' 
    subprocess.call(cmd, shell=True, executable='/bin/bash')
    cmd = 'samtools index '+args.prefix+'.Chimeric.out.sorted.bam'
    subprocess.call(cmd, shell=True, executable='/bin/bash')
    os.remove(args.prefix+'.Chimeric.out.sam')
