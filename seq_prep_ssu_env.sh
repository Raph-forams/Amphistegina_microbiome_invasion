##########################################################
#                                                        #
#  SSU envinronmental samples eukaryotes + prokaryotes   #
#    + re-done libraries 7 & 8                           #
#                                                        #
##########################################################

# Libraries analyzed in this run
#   Lib n° - concentration -      target     - sample type
#   Lib 7  -     10%       -      proks      - single forams
#   Lib 8  -     10%       -      proks      - single forams
#   Lib 13 -     20%       -      euks       - environmental (sediment & filter)
#   Lib 14 -     20%       -      euks       - environmental (sediment & filter)
#   Lib 15 -     20%       -      proks      - environmental (sediment & filter)
#   Lib 16 -     20%       -      proks      - environmental (sediment & filter)


# set working directory on server
WDIR="/storage/hdd1/chh/Debora_amplicons/Environmental" 
cd $WDIR

# activate conda environment
conda activate trimming
# cutadapt 3.1

# create directory for logfiles
mkdir Logfiles

# CAUTION: It is not recommended to combine the output from several sequencer runs. Some argue, lanes should also be processed independently.
# However, as each library only contains a limited number of samples (with a sequencing depth <100K per sample), we will process all libraries in a project together regardless of sequencing run.
# Otherwise, the number of sequences per library may not be sufficient for good denoising

# prepare index files
# split index information by library
# and prepare index fasta for demultiplexing
# make sure that there are no spaces in the table with the index and primer information
for library in $(sed '1d' Environmental_and_repeat_library_info.txt | cut -f6 | sort | uniq)
do
  sed '1d' Environmental_and_repeat_library_info.txt | awk -v FS="\t" -v OFS="\t" -v library=${library} '$6 == library' > "Library_"${library}"_index_primer.txt"
done


### Demultiplex using cutadapt version 3.0

# create directory for demultiplexed files
mkdir Demux

# as illumina adapters were ligated regardless of orientation, demultiplex both fwd-rev and rev-fwd
# treat them separately until after denoising and merging
for i in $(ls -1 Library_*)
do
  LIB=$(echo "${i}" | sed -e 's/^Library_//' -e 's/_index_primer.txt$//')
  while read line
  do
    SID=$(echo "${line}" | cut -f1)
    FWD=$(echo "${line}" | cut -f2 | sed 's/^/\^/')
    REV=$(echo "${line}" | cut -f4 | sed 's/^/\^/')
    OLP=7 # index length minus 1
    
    cutadapt -j 2 -O ${OLP} --no-indels -e 0 -g ${FWD} -G ${REV} --discard-untrimmed -o Demux/${SID}"_fr_R1.fastq.gz" -p Demux/${SID}"_fr_R2.fastq.gz" Libraries/Library${LIB}_S*_L001_R1_001.fastq.gz Libraries/Library${LIB}_S*_L001_R2_001.fastq.gz > Logfiles/${SID}".demux_fr.log" 2>&1
    
    cutadapt -j 2 -O ${OLP} --no-indels -e 0 -g ${REV} -G ${FWD} --discard-untrimmed -o Demux/${SID}"_rf_R1.fastq.gz" -p Demux/${SID}"_rf_R2.fastq.gz" Libraries/Library${LIB}_S*_L001_R1_001.fastq.gz Libraries/Library${LIB}_S*_L001_R2_001.fastq.gz > Logfiles/${SID}".demux_rf.log" 2>&1 
  done < ${i}
done


# further processing steps are easier and faster without gz
gunzip Demux/*.gz

# from here on work with three different approaches:
# 1) Environmental eukaryotes: Library 13, 14
# 2) Environmental prokaryotes: Library 15, 16
# 3) Re-done single forams prokaryotes: Library 7, 8


### Primer clipping using cutadapt version 3.0

# select only the libraries we want 
cat Library_13_index_primer.txt Library_14_index_primer.txt | cut -f1 > sample_names_13_14.txt
cat Library_15_index_primer.txt Library_16_index_primer.txt | cut -f1 > sample_names_15_16.txt
cat Library_7_index_primer.txt Library_8_index_primer.txt | cut -f1 > sample_names_7_8.txt

# create directory for primer-clipped files
mkdir Clipped

# Primer clipping 1) Environmental eukaryotes: Library 13, 14

# set clipping parameters (V4 - eukaryotes)
FWD="^CCAGCASCYGCGGTAATTCC"
REV="^ACTTTCGTTCTTGATYRA"
OFWD=$(expr ${#FWD} - 2)
OREV=$(expr ${#REV} - 2)
ERROR=0.16

# remove primer from each sample
while read line
do
  # process fwd-rev orientation
  cutadapt -j 2 --no-indels -e ${ERROR} -g "${FWD};o=${OFWD}" -G "${REV};o=${OREV}" -m 50 --discard-untrimmed -o Clipped/${line}"_clip_fr_R1.fastq" -p Clipped/${line}"_clip_fr_R2.fastq" Demux/${line}"_fr_R1.fastq" Demux/${line}"_fr_R2.fastq" > Logfiles/${line}".clip_fr.log" 2>&1
  # process rev-fwd orientation
  cutadapt -j 2 --no-indels -e ${ERROR} -g "${REV};o=${OREV}" -G "${FWD};o=${OFWD}" -m 50 --discard-untrimmed -o Clipped/${line}"_clip_rf_R1.fastq" -p Clipped/${line}"_clip_rf_R2.fastq" Demux/${line}"_rf_R1.fastq" Demux/${line}"_rf_R2.fastq" > Logfiles/${line}".clip_rf.log" 2>&1
done < sample_names_13_14.txt

# Primer clipping 2) Environmental prokaryotes: Library 15, 16 and 3) Re-done single forams prokaryotes: Library 7, 8

# set clipping parameters (prokaryotes)
FWD="^CCTACGGGNGGCWGCAG"
REV="^GACTACHVGGGTATCTAATCC"
OFWD=$(expr ${#FWD} - 2)
OREV=$(expr ${#REV} - 2)
ERROR=0.16

# remove primer from each sample
cat sample_names_15_16.txt sample_names_7_8.txt | while read line
do
  # process fwd-rev orientation
  cutadapt -j 2 --no-indels -e ${ERROR} -g "${FWD};o=${OFWD}" -G "${REV};o=${OREV}" -m 50 --discard-untrimmed -o Clipped/${line}"_clip_fr_R1.fastq" -p Clipped/${line}"_clip_fr_R2.fastq" Demux/${line}"_fr_R1.fastq" Demux/${line}"_fr_R2.fastq" > Logfiles/${line}".clip_fr.log" 2>&1
  # process rev-fwd orientation
  cutadapt -j 2 --no-indels -e ${ERROR} -g "${REV};o=${OREV}" -G "${FWD};o=${OFWD}" -m 50 --discard-untrimmed -o Clipped/${line}"_clip_rf_R1.fastq" -p Clipped/${line}"_clip_rf_R2.fastq" Demux/${line}"_rf_R1.fastq" Demux/${line}"_rf_R2.fastq" > Logfiles/${line}".clip_rf.log" 2>&1
done 


### Count sequences for each step and switch to dada2 in R
ls -1v Demux/*_fr_R1.fastq | xargs wc -l | grep -v "total" | awk '{print $1/4}' | paste <(ls -1v Demux/*_fr_R1.fastq | xargs -n1 basename | sed 's/_fr_R1\.fastq//') - > tmp1
ls -1v Clipped/*_clip_fr_R1.fastq | xargs wc -l | grep -v "total" | awk '{print $1/4}' | paste tmp1 - > tmp2
echo -e 'SID\tDemux\tClipped' | cat - tmp2 > nSeqs_fr.txt
rm tmp*
ls -1v Demux/*_rf_R1.fastq | xargs wc -l | grep -v "total" | awk '{print $1/4}' | paste <(ls -1v Demux/*_rf_R1.fastq | xargs -n1 basename | sed 's/_rf_R1\.fastq//') - > tmp1
ls -1v Clipped/*_clip_rf_R1.fastq | xargs wc -l | grep -v "total" | awk '{print $1/4}' | paste tmp1 - > tmp2
echo -e 'SID\tDemux\tClipped' | cat - tmp2 > nSeqs_rf.txt
rm tmp*

