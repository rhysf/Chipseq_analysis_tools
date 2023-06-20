#!/bin/sh -l
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source $DIR/Bash_functions/File_input_output.sh

### Michele Busby, Jean Chang, Catherine Li, Tal Goodisman & Rhys Farrer

function echo_and_print_to_outfile {
	echo "$@" 
	[ -d ${opts[o]} ] || mkdir ${opts[o]}
	echo "$@" >> ${opts[o]}/commands.txt 2>&1 
}

# Exit with message
diemsg() { 
	echo "Usage:    $0 -g <genome.fasta> -b <folder with BAM replicates> -c <control BAM merged>"
	echo "Optional: -o Output directory [./]"
	echo "          -p Peak style (factor/histone) [histone]"
	echo "          -a Alignment style (paired=TRUE, unpaired=FALSE) [TRUE]"
	echo "          -f fold change [4]"
	echo "          -d FDR [0.001]"
	echo "          -v P-value [0.0001]"
	echo "          -e Effective genome size [2e9]"
	echo "          -s BAM replicates suffix [clean.bam]"
	echo ""
	echo "Notes: requires homer peak caller installed: http://homer.ucsd.edu/homer/download.html"
	echo "       requires igvtools: https://software.broadinstitute.org/software/igv/download"
	echo "       requires bedtools and samtools."
	die  "Notes2: merge: samtools merge -o merged_Chip_input.sorted.bam -b merged_Chip_input.sorted.bam.txt"
}
[ $# -gt 6 ] || diemsg

# Gather input files and settings
declare -A opts
opts[o]="./"
opts[p]="histone"
opts[a]="TRUE"
opts[f]=4
opts[d]=0.001
opts[v]=0.0001
opts[e]="2e9"
opts[s]="clean.bam"
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
echo_and_print_to_outfile "Running $0 -g ${opts[g]} -b ${opts[b]} -c ${opts[c]} -o ${opts[o]} -p ${opts[p]} -a ${opts[a]} -f ${opts[f]} -d ${opts[d]} -v ${opts[v]} -e ${opts[e]} -s ${opts[s]}"

# Get file names
fileNames=( $( ls ${opts[b]}/*clean.bam) )
echo_and_print_to_outfile "BAM replicates = ${fileNames}";

### Functions

Make_tag_directory_for_control_bam() {
	CONTROL=$1
	OUTPUT_DIR=$2
	GENOME_FILE=$3
	PAIRED_END=$4
	OPTIONAL_EXTRA_PARAMS=$5

	if [ "$CONTROL"  != "NONE" ]  && [ -f  "$CONTROL" ]; then
		echo_and_print_to_outfile "Make_tag_directory_for_control_bam: Control $CONTROL is a file, assuming bam"
		CONTROL_STEM=$(basename "${CONTROL}" .bam)
		CONTROL_BAM=$CONTROL
		CONTROL_DIR="$OUTPUT_DIR/Tags/$CONTROL_STEM"

		# cmd
		CMD1="makeTagDirectory $CONTROL_DIR $CONTROL_BAM -genome $GENOME_FILE -checkGC"

		# add paired end mode, and extra param
		if [ "$PAIRED_END" == "TRUE" ]; then
			CMD1="${CMD1} -illuminaPE $OPTIONAL_EXTRA_PARAMS" # why only optional extra params here? check
		fi

		echo_and_print_to_outfile "Make_tag_directory_for_control_bam: Paired end = $PAIRED_END"
		echo_and_print_to_outfile "CMD1=$CMD1"
		eval $CMD1 >> ${opts[o]}/commands.txt 2>&1 

	elif [ "$CONTROL"  != "NONE" ]  &&  [ -d  "$CONTROL" ]; then
		echo_and_print_to_outfile "Control is a directory, assuming processed homer peaks"
		CONTROL_DIR="$CONTROL";
		echo_and_print_to_outfile "Control dir is $CONTROL_DIR "
	else
		echo_and_print_to_outfile "Nothing found for control $CONTROL"
	fi
}

prepareoutdirs() {
	INPUT_BAM=$1
	OUTPUT_DIR=$2

	echo_and_print_to_outfile "prepareoutdirs: Input bam: $INPUT_BAM"
	echo_and_print_to_outfile "prepareoutdirs: Writing to: $OUTPUT_DIR"

	#If output directory does not exist make it
	if [ -d $OUTPUT_DIR ] ; then
   		echo_and_print_to_outfile "Output directory $OUTPUT_DIR exists."
	else
   		echo_and_print_to_outfile "Output directory $OUTPUT_DIR does not exists. Creating directory."
   		mkdir $OUTPUT_DIR
	fi

	if [ -d $OUTPUT_DIR/Tags ] ; then
		echo_and_print_to_outfile "Output directory $OUTPUT_DIR exists."
	else
   		echo_and_print_to_outfile "Output directory $OUTPUT_DIR does not exists. Creating directory."
   		mkdir $OUTPUT_DIR/Tags
	fi

	if [ -d $OUTPUT_DIR/Peaks ] ; then
   		echo_and_print_to_outfile "Output directory $OUTPUT_DIR/Peaks exists."
	else
   		echo_and_print_to_outfile "Output directory $OUTPUT_DIR/Peaks  does not exists. Creating directory."
   		mkdir $OUTPUT_DIR/Peaks
	fi

	if [ -d $OUTPUT_DIR/Peaks/bed ] ; then
   		echo_and_print_to_outfile "Output directory $OUTPUT_DIR/Peaks/bed exists."
	else
   		echo "Output directory $OUTPUT_DIR/Peaks/bed  does not exists. Creating directory."
   		mkdir $OUTPUT_DIR/Peaks/bed
	fi

	if [ -d $OUTPUT_DIR/tdfs ] ; then
   		echo_and_print_to_outfile "Output directory $OUTPUT_DIR/tdfs exists."
	else
   		echo_and_print_to_outfile "Output directory $OUTPUT_DIR/tdfs does not exists. Creating directory."
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
	echo_and_print_to_outfile "CMD2=$CMD2"
	eval $CMD2 >> ${opts[o]}/commands.txt 2>&1 

	# Homer has two parts: makes a tag directory then does peak calling
	# Make tag directory for input bam
}

Make_tag_directory_for_input_bam() {
	INPUT_BAM=$1
	OUTPUT_DIR=$2
	GENOME_FILE=$3
	PAIRED_END=$4

	CMD3="makeTagDirectory $OUTPUT_DIR/Tags/$STEM $INPUT_BAM -genome $GENOME_FILE -checkGC"

	# Paired end mode
	if [ "$PAIRED_END" == "TRUE" ]; then
		CMD3="${CMD3} -illuminaPE"
	fi

	echo_and_print_to_outfile "$0: Paired end $PAIRED_END"
	echo_and_print_to_outfile "CMD3=$CMD3"
   	eval $CMD3 >> ${opts[o]}/commands.txt 2>&1 
}

Find_peaks() {
	INPUT_BAM=$1
	OUTPUT_DIR=$2
	CONTROL_DIR=$3
	PEAK_STYLE=$4
	FOLD_CHANGE=$5
	FDR=$6
	P_VALUE=$7
	EFFECTIVE_GENOME_SIZE=$8

	# Get stem
	STEM=$(basename "${INPUT_BAM}" .bam)

	# Add control dir
	CMD4="findPeaks $OUTPUT_DIR/Tags/$STEM -style $PEAK_STYLE -o $OUTPUT_DIR/Peaks/$STEM.calls"
	if [ "$CONTROL_DIR" != "NONE" ]; then
		CMD4="${CMD4} -i $CONTROL_DIR"
	fi

	# Add settings
	SETTINGS="-F ${FOLD_CHANGE} -fdr ${FDR} -P ${P_VALUE} -gsize {$EFFECTIVE_GENOME_SIZE}"
	$CMD4="${CMD4} ${SETTINGS}"

	# Find peaks
   	echo_and_print_to_outfile "CMD4=$CMD4"
   	eval $CMD4 >> ${opts[o]}/commands.txt 2>&1 

	# Make calls file a bed
	CMD5="cut -f 2,3,4  $OUTPUT_DIR/Peaks/$STEM.calls> $OUTPUT_DIR/Peaks/bed/$STEM.bed"
	echo_and_print_to_outfile "CMD5=$CMD5"
	eval $CMD5 >> ${opts[o]}/commands.txt 2>&1 

	CMD6="sed -i '' 's/^chr\t/#chr\t/g' $OUTPUT_DIR/Peaks/bed/$STEM.bed"
	echo_and_print_to_outfile "CMD6=$CMD6"
	eval $CMD6 >> ${opts[o]}/commands.txt 2>&1 

	# Get the SPOT calls
	echo_and_print_to_outfile "Writing spot score to $OUTPUT_DIR/Peaks/spotscores.txt"
	readsInPeaks=$(bedtools intersect -abam $INPUT_BAM -b $OUTPUT_DIR/Peaks/bed/$STEM.bed -bed | wc -l)
	readsInAlignment=$(samtools view $INPUT_BAM | wc -l)
	echo_and_print_to_outfile "Reads in Alignment: $readsInAlignment"
	spotScore=$(echo "scale=4; $readsInPeaks/$readsInAlignment" |bc)
	echo "$STEM $spotScore" >> $OUTPUT_DIR/Peaks/spotscores.txt
}

# Make tag directory for control bams
#EXTRA_PARAM="-tbp 1"
EXTRA_PARAM=
Make_tag_directory_for_control_bam ${opts[c]} ${opts[o]} ${opts[g]} ${opts[a]} $EXTRA_PARAM

echo_and_print_to_outfile Finished processing ${opts[c]}

for f in ${fileNames[@]}; do

	echo_and_print_to_outfile "processing filename $f" 

	if [ $f != ${opts[c]} ]; then
		STEM=$(basename "${f}" .bam)
		echo_and_print_to_outfile "STEM is $STEM" # checking output

		# Run pipeline
		prepareoutdirs $f ${opts[o]}
		make_bam_a_tdf $f ${opts[o]} ${opts[g]}
		Make_tag_directory_for_input_bam $f ${opts[o]} ${opts[g]} ${opts[a]}
		EXTRA_PARAM2=
		Make_tag_directory_for_control_bam $CONTROL_DIR ${opts[o]} ${opts[g]} ${opts[a]} $EXTRA_PARAM2
		Find_peaks $f ${opts[o]} ${CONTROL_DIR} ${opts[p]} ${opts[f]} ${opts[d]} ${opts[v]} ${opts[e]}
	fi
done