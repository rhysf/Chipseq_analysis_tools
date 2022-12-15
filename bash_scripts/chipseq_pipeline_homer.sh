#!/bin/sh -l
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source $DIR/Bash_functions/File_input_output.sh

### Michele Busby, Jean Chang, Catherine Li, Tal Goodisman & Rhys Farrer

# Exit with message
diemsg() { 
	echo "Usage:    $0 -g <genome.fasta> -b <folder with BAM replicates> -c <control BAM merged>"
	echo "Optional: -o Output directory [./]"
	echo "          -p Peak style (factor/histone) [histone]"
	echo "          -a Alignment style (paired=TRUE, unpaired=FALSE) [TRUE]"
	echo "          -f fold enrichment over input tag count [4]"
	echo "          -q poisson p-value threshold relative to input tag count [0.0001]"
	echo "          -d False discovery rate [0.001]"
	echo ""
	echo "Notes: requires homer peak caller installed: http://homer.ucsd.edu/homer/download.html"
	echo "       requires igvtools: https://software.broadinstitute.org/software/igv/download"
	echo "       requires bedtools and samtools."
	die  "Notes: merge: samtools merge -o merged_Chip_input.sorted.bam -b merged_Chip_input.sorted.bam.txt"
}
[ $# -gt 1 ] || diemsg

# Gather input files and settings
declare -A opts
opts[o]="./"
opts[p]="histone"
opts[a]="TRUE"
opts[f]=4
opts[q]=0.0001
opts[d]=0.001
getOpts $@

# Check input files are readable
for file in ${opts[g]} ${opts[c]} ; do
	[[ -r $file ]] || die "File: $file does not appear to be valid"
done

# Check input folders are present
for file in ${opts[b]} ; do
	[[ -d $file ]] || die "Folder: $file does not appear to be valid"
done

# Warning
echo "Running $0 -g ${opts[g]} -b ${opts[b]} -c ${opts[c]} -o ${opts[o]} -p ${opts[p]} -a ${opts[a]}"

### Functions

Make_tag_directory_for_control_bam() {
	CONTROL=$1
	OUTPUT_DIR=$2
	GENOME_FILE=$3
	PAIRED_END=$4
	OPTIONAL_EXTRA_PARAMS=$5

	if [ "$CONTROL"  != "NONE" ]  && [ -f  "$CONTROL" ]; then
		echo "Control is a file, assuming bam"
		CONTROL_STEM=$(basename "${CONTROL}" .bam)
		CONTROL_BAM=$CONTROL;
		CONTROL_DIR="$OUTPUT_DIR/Tags/$CONTROL_STEM"

		if [ "$PAIRED_END" == "TRUE" ]; then
			echo "Paired end"
			CMD1="makeTagDirectory $CONTROL_DIR $CONTROL_BAM -genome $GENOME_FILE  -illuminaPE -checkGC $OPTIONAL_EXTRA_PARAMS"
      			echo "CMD1=$CMD1"
      			eval $CMD1
		else
			echo "Not paired end"
			CMD1="makeTagDirectory $CONTROL_DIR $CONTROL_BAM -genome $GENOME_FILE -checkGC"
      			echo "CMD1=$CMD1"
      			eval $CMD1
		fi
	elif [ "$CONTROL"  != "NONE" ]  &&  [ -d  "$CONTROL" ]; then
		echo "Control is a directory, assuming processed homer peaks"
		CONTROL_DIR="$CONTROL";
		echo "Control dir is $CONTROL_DIR "
	else
		echo "Nothing found for control $CONTROL"
	fi
}

prepareoutdirs() {
	INPUT_BAM=$1
	OUTPUT_DIR=$2

	echo "prepareoutdirs: Input bam: $INPUT_BAM"
	echo "prepareoutdirs: Writing to: $OUTPUT_DIR"

	#If output directory does not exist make it
	if [ -d $OUTPUT_DIR ] ; then
   		echo "Output directory $OUTPUT_DIR exists."
	else
   		echo "Output directory $OUTPUT_DIR does not exists. Creating directory."
   		mkdir $OUTPUT_DIR
	fi

	if [ -d $OUTPUT_DIR/Tags ] ; then
    	echo "Output directory $OUTPUT_DIR exists."
	else
   		echo "Output directory $OUTPUT_DIR does not exists. Creating directory."
   		mkdir $OUTPUT_DIR/Tags
	fi

	if [ -d $OUTPUT_DIR/Peaks ] ; then
   		echo "Output directory $OUTPUT_DIR/Peaks exists."
	else
   		echo "Output directory $OUTPUT_DIR/Peaks  does not exists. Creating directory."
   		mkdir $OUTPUT_DIR/Peaks
	fi

	if [ -d $OUTPUT_DIR/Peaks/bed ] ; then
   		echo "Output directory $OUTPUT_DIR/Peaks/bed exists."
	else
   		echo "Output directory $OUTPUT_DIR/Peaks/bed  does not exists. Creating directory."
   		mkdir $OUTPUT_DIR/Peaks/bed
	fi

	if [ -d $OUTPUT_DIR/tdfs ] ; then
   		echo "Output directory $OUTPUT_DIR/tdfs exists."
	else
   		echo "Output directory $OUTPUT_DIR/tdfs does not exists. Creating directory."
   		mkdir $OUTPUT_DIR/tdfs
	fi
}

make_bam_a_tdf() {
	INPUT_BAM=$1
	OUTPUT_DIR=$2
	GENOME_FILE=$3
	PAIRED_END=$4

	#Get stem
	STEM=$(basename "${INPUT_BAM}" .bam)

	#Make bam a tdf
	CMD2="igvtools count -e 100 $INPUT_BAM $OUTPUT_DIR/tdfs/$STEM.tdf $GENOME_FILE"
	echo "CMD2=$CMD2"
	eval $CMD2

	# Homer has two parts: makes a tag directory then does peak calling
	# Make tag directory for input bam
}

Make_tag_directory_for_input_bam() {
	INPUT_BAM=$1
	OUTPUT_DIR=$2
	GENOME_FILE=$3
	PAIRED_END=$4

	if [ "$PAIRED_END" == "TRUE" ]; then
		echo "Paired end"
		CMD3="makeTagDirectory $OUTPUT_DIR/Tags/$STEM $INPUT_BAM -genome $GENOME_FILE -checkGC -illuminaPE"
   		echo "CMD3=$CMD3"
   		eval $CMD3
	else
		echo "Not paired end"
		CMD3="makeTagDirectory $OUTPUT_DIR/Tags/$STEM $INPUT_BAM -genome $GENOME_FILE -checkGC"
   		echo "CMD3=$CMD3"
   		eval $CMD3
	fi
}

Find_peaks() {
	INPUT_BAM=$1
	OUTPUT_DIR=$2
	CONTROL_DIR=$3
	PEAK_STYLE=$4
	FOLD=$5
	PVALUE_RELATIVE_TO_INPUT=$6
	FDR=$7

	# Get stem
	STEM=$(basename "${INPUT_BAM}" .bam)

	# Find peaks (homer software e.g., /Users/rhysfarrer/Desktop/Programs/homer_peak_caller/bin/findPeaks)
	if [ "$CONTROL_DIR" == "NONE" ]; then
		CMD4="findPeaks $OUTPUT_DIR/Tags/$STEM -style $PEAK_STYLE -o $OUTPUT_DIR/Peaks/$STEM.calls -F $FOLD -P $PVALUE_RELATIVE_TO_INPUT -fdr $FDR"
   		echo "CMD4=$CMD4"
   		eval $CMD4
	else
		CMD4="findPeaks $OUTPUT_DIR/Tags/$STEM -style $PEAK_STYLE -o $OUTPUT_DIR/Peaks/$STEM.calls -i $CONTROL_DIR -F $FOLD -P $PVALUE_RELATIVE_TO_INPUT -fdr $FDR"
   		echo "CMD4=$CMD4"
   		eval $CMD4
	fi

	# Make calls file a bed
	CMD5="cut -f 2,3,4  $OUTPUT_DIR/Peaks/$STEM.calls> $OUTPUT_DIR/Peaks/bed/$STEM.bed"
	echo "CMD5=$CMD5"
	eval $CMD5

	CMD6="sed -i '' 's/^chr\t/#chr\t/g' $OUTPUT_DIR/Peaks/bed/$STEM.bed"
	echo "CMD6=$CMD6"
	eval $CMD6

	# Get the SPOT calls
	echo "Writing spot score to $OUTPUT_DIR/Peaks/spotscores.txt"
	readsInPeaks=$(bedtools intersect -abam $INPUT_BAM -b $OUTPUT_DIR/Peaks/bed/$STEM.bed -bed | wc -l)
	readsInAlignment=$(samtools view $INPUT_BAM | wc -l)
	echo "Reads in Alignment: $readsInAlignment"
	spotScore=$(echo "scale=4; $readsInPeaks/$readsInAlignment" |bc)
	echo "$STEM $spotScore" >> $OUTPUT_DIR/Peaks/spotscores.txt
}

# Get file names
fileNames=( $( ls ${opts[b]}/*.bam) )

# Make tag directory for control bams
Make_tag_directory_for_control_bam ${opts[c]} ${opts[o]} ${opts[g]} ${opts[a]} "-tbp 1"
echo Finished processing ${opts[c]}

for f in ${fileNames[@]}; do

	echo processing filename $f

	if [ $f != ${opts[c]} ]; then
		STEM=$(basename "${f}" .bam)
		echo "STEM is $STEM" # checking output

		# Run pipeline
		prepareoutdirs $f ${opts[o]}
		make_bam_a_tdf $f ${opts[o]} ${opts[g]}
		Make_tag_directory_for_input_bam $f ${opts[o]} ${opts[g]} ${opts[a]}
		EXTRA_PARAM=
		Make_tag_directory_for_control_bam $CONTROL_DIR ${opts[o]} ${opts[g]} ${opts[a]} $EXTRA_PARAM
		Find_peaks $f ${opts[o]} ${CONTROL_DIR} ${opts[p]} ${opts[f]} ${opts[q]} ${opts[d]}
	fi
done
