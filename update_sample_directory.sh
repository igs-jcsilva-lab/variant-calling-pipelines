#! /bin/bash
#$ -S /bin/bash
#$ -P jcsilva-gcid-proj4a-malaria
#$ -cwd

#Updates: can now run on rhel6 machines
#	  Becuase of this, now running updated versions of:
#	  java, gatk, bowtie2
			
###############################################################################

source $SGE_ROOT/igs/common/settings.sh

###############################################################################

set 1>&2

set -o errexit
set -o nounset
set -o pipefail
set -o xtrace

###############################################################################

export PATH=/home/edrabek/bin:$PATH
export PYTHONPATH=/home/edrabek/lib/python
export LD_LIBRARY_PATH=/home/edrabek/lib
export LC_ALL=C

#export PATH=/usr/local/packages/tabix:$PATH

###############################################################################

temp_d=`mktemp -d`

function clean_up {
  echo 1>&2 "? = $?"
  rm -rf $temp_d
}

trap clean_up 0

###############################################################################
# command line parameters

# note the assumption that the second-to-last element of the directory name is
# the country and the last element of the directory name is the sample ID, and
# that this latter is unique (i.e., you don't have two samples from different
# countries with the same ID)

sample=`basename $1`
country=`basename $(dirname $1)`
cd $1
shift
scripts_d=$1
shift
auxiliary_files_d=$1
shift

###############################################################################
# make a list of sequencing runs and compile the correct files to use for each one, always preferring the trimmed
# version when present

f=sequencing_runs
if [[ ! -e $f ]]; then
  {
    pushd ILLUMINA_DATA > /dev/null
    ls *.fastq.gz \
      | perl -pe 's/(.*)_R[12](_trimmed)?.fastq.gz/$1/ or die $_' \
      | sort -u
    popd > /dev/null
  } > $f
fi

mkdir -p ILLUMINA_DATA.correct_files_to_use

pushd ILLUMINA_DATA.correct_files_to_use
for sequencing_run in `cat ../sequencing_runs`; do
  for mate in 1 2; do
    for suffix in .fastq.gz _stats.txt; do
      f=${sequencing_run}_R$mate$suffix
      if [[ ! -e $f ]]; then
        orig_f=../ILLUMINA_DATA/$f
        trimmed_file=`echo $orig_f | perl -pe 's/(_R[12])(.fastq.gz|_stats.txt)/$1_trimmed$2/'`
        if [[ -e $trimmed_file ]]; then
          orig_f=$trimmed_file
        fi

        if [[ -e $orig_f ]]; then
          ln -s $orig_f $f
        fi
      fi
    done


    # fill in the stats file when it's not present
    if [[ ! -e ${sequencing_run}_R${mate}_stats.txt ]]; then
      zcat ${sequencing_run}_R${mate}.fastq.gz \
        | $scripts_d/get_illumina_stats_rhel7 \
        > $f
    fi

    # compute mean quality
    f=${sequencing_run}_R$mate.base_count-mean_quality
    if [[ ! -e $f ]]; then
      zcat ${sequencing_run}_R${mate}.fastq.gz \
        | $scripts_d/compute_mean_quality_rhel7 \
        > $f
    fi
  done
done
popd

###############################################################################
# Align individual sequencing runs
# Adapted from snpcalling_prep_for_IGS.sh

ref=$auxiliary_files_d/reference.fa

java="/usr/bin/java"
bowtie2="/usr/local/packages/bowtie2-2.2.9/bowtie2"

mkdir -p alignments_v24

pushd alignments_v24
for sequencing_run in `cat ../sequencing_runs`; do
  mkdir -p $sequencing_run
  pushd $sequencing_run

  mate_1_f=../../ILLUMINA_DATA.correct_files_to_use/${sequencing_run}_R1.fastq.gz
  mate_2_f=../../ILLUMINA_DATA.correct_files_to_use/${sequencing_run}_R2.fastq.gz

  # If there are exactly five underscores in the sequencing run ID, assume that it encodes the library, barcode, and lane
  if [[ _____ == "`echo $sequencing_run | perl -pe 's/[^_\n]+//g'`" ]]; then
    library=$(echo $sequencing_run | cut -f4 -d'_')
    barcode=$(echo $sequencing_run | cut -f5 -d'_')
    lane=$(echo $sequencing_run | cut -f6 -d'_')
  else
    library=unknown
    barcode=unknown
    lane=unknown
  fi

  if [[ ! -e recalibrated.bam ]]; then
    # Alignment with bowtie2, creating bam file
    $bowtie2 --no-unal -x $auxiliary_files_d/reference_index -1 $mate_1_f -2 $mate_2_f \
      | samtools view - -S -q 0 -b \
      > raw.bam

    # Fixing read groups/header/sorting
    ${java} -Xmx16g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups I=raw.bam o=fixed.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT RGID=${sequencing_run} RGLB=${library} RGPL=Illumina RGPU=${barcode} RGSM="$country-$sample"

    rm raw.bam

    # Sort bam file (for deduping step
    ${java} -jar /usr/local/packages/picard-tools-2.5.0/picard.jar SortSam INPUT=fixed.bam OUTPUT=sorted.bam SORT_ORDER=coordinate

    rm fixed.bam

    # Mark duplicates
    ${java} -Xmx16g -jar /usr/local/packages/picard-tools-2.5.0/picard.jar MarkDuplicates INPUT=sorted.bam OUTPUT=dedup.bam METRICS_FILE=dedup_metrics.txt

    rm sorted.bam dedup_metrics.txt

    # Re-index with SAM tools
    samtools index dedup.bam

    # Indels- Could skip this if we don't care about SNPs in indels
    ${java} -Xmx16g -jar /usr/local/packages/gatk-3.5/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I dedup.bam -o indel.intervals

    # Interval + ref + BAM:
    ${java} -Xmx16g -jar /usr/local/packages/gatk-3.5/GenomeAnalysisTK.jar -T IndelRealigner -targetIntervals indel.intervals -o realigned.bam -I dedup.bam -R $ref

    rm dedup.ba*
    rm indel.intervals

    # Base-recalibration:
    ${java} -Xmx16g -jar /usr/local/packages/gatk-3.5/GenomeAnalysisTK.jar -T BaseRecalibrator -R $ref -I realigned.bam -knownSites $auxiliary_files_d/known_snps.vcf -o recal.grp
    ${java} -Xmx16g -jar /usr/local/packages/gatk-3.5/GenomeAnalysisTK.jar -T PrintReads -R $ref -I realigned.bam -BQSR recal.grp -o recalibrated.bam

    rm realigned.ba*
    rm recal.grp
  fi

  # Generating alignment and coverage information
  f=recalibrated.coverage.gz
  if [[ ! -e $f ]]; then
    bedtools genomecov -ibam recalibrated.bam -g $ref \
      | gzip -v \
      > $f
  fi

  f=recalibrated.flagstat
  if [[ ! -e $f ]]; then
    samtools flagstat recalibrated.bam \
      > recalibrated.flagstat
  fi

  popd
done
popd

###############################################################################
# SNP calling
# Adapted from snpcalling_submission_unknown_MPnum.sh and snpcalling_for_IGS.sh
#adding preliminary steps for transitioning to Haplotype Caller on July 18, 2016 - kmoser

java="/usr/local/packages/jdk-8u151/bin/java"

gatk="/usr/local/packages/gatk-4.0.4.0/gatk"

mkdir -p snpcalls_v24
pushd snpcalls_v24

f=bam.list
if [[ ! -e $f ]]; then
  ls ../alignments_v24/*/recalibrated.bam > bam.list
fi

f=diploid.g.vcf
if [[ ! -e $f ]]; then
  $gatk --java-options "-Xmx8G" HaplotypeCaller \
    -R $auxiliary_files_d/reference.fa \
    -I bam.list \
    -ERC GVCF \
    -O $f
fi

popd

###############################################################################
# Merge individual bam files and get genome coverage and flagstats

mkdir -p alignments_v24.merged
pushd alignments_v24.merged

prefix=all.sorted_and_recalibrated

f=$prefix.bam
if [[ ! -e $f ]]; then
  mkdir -p $(dirname $f)


  recalibrated_count=`ls ../alignments_v24/* | grep -c recalibrated.bam || true`
  if [[ $recalibrated_count == 0 ]]; then
    echo 1>&2 Error: No recalibrated files in `pwd`/../alignments_v24
    exit 1
  elif [[ $recalibrated_count == 1 ]]; then
    ln -s ../alignments_v24/*/recalibrated.bam $f
  else
    samtools merge $f ../alignments_v24/*/recalibrated.bam
  fi
fi

f=$prefix.coverage.gz
if [[ ! -e $f ]]; then
  bedtools genomecov -ibam $prefix.bam -g $auxiliary_files_d/reference.contig-length \
    | gzip -v \
    > $f
fi

f=$prefix.position-coverage.gz
if [[ ! -e $f ]]; then
  bedtools genomecov -d -ibam $prefix.bam -g $auxiliary_files_d/reference.contig-length \
    | gzip -v \
    > $f
fi

f=$prefix.gene-coverage.gz
if [[ ! -e $f ]]; then
  bedtools coverage -abam $prefix.bam -b $auxiliary_files_d/reference.gff \
    | gzip -v \
    > $f
fi

f=$prefix.flagstat
if [[ ! -e $f ]]; then
  samtools flagstat $prefix.bam \
  > $f
fi

popd

###############################################################################
# TODO: UPDATE ALL GATK TO 4
