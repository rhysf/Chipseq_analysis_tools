#!/bin/bash -l

### rfarrer@broadinstitute.org

# Exit with message
die() { echo "$@" ; exit 1; }

# Get command line parameters (opts global variable)
function getOpts {
	local input=$@

	# Gather input files and settings
	while [ ${#input[@]} -gt 0 ]
	do
		case "$1" in
        		-a) opts["a"]="$2";; 
        		-b) opts["b"]="$2";;
        		-c) opts["c"]="$2";;
        		-d) opts["d"]="$2";;
        		-e) opts["e"]="$2";;
        		-f) opts["f"]="$2";;
        		-g) opts["g"]="$2";;
        		-h) opts["h"]="$2";;
        		-i) opts["i"]="$2";;
        		-j) opts["j"]="$2";;
        		-k) opts["k"]="$2";;
        		-l) opts["l"]="$2";;
        		-m) opts["m"]="$2";;
        		-n) opts["n"]="$2";;
        		-o) opts["o"]="$2";;
        		-p) opts["p"]="$2";;
        		-q) opts["q"]="$2";;
        		-r) opts["r"]="$2";;
        		-s) opts["s"]="$2";;
        		-t) opts["t"]="$2";;
        		-u) opts["u"]="$2";;
        		-v) opts["v"]="$2";;
        		-w) opts["w"]="$2";;
        		-x) opts["x"]="$2";;
        		-y) opts["y"]="$2";;
        		-z) opts["z"]="$2";;
			-*) echo >&2 \
				die "Unusual option given in getOpts in $0";;
			*)  break;;	# terminate while loop
		esac
		shift
		shift
	done

	#echo "In File_input_output.sh and -a =" 
	#echo ${opts[a]}
	#echo "In File_input_output.sh and -b =" 
	#echo ${opts[b]}
}

# Folder with R1 and R2 only (for batch analysis) (R1 and R2 global variables)
function findFastqPairsInFolder {
	local input=$1
	echo "Checking for R1 and R2 fastq files..."
	yourfilenames=`ls $1/*`
	for eachfile in $yourfilenames
	do
		if [[ ${eachfile} =~ _R1.fastq$ ]] ; then
			echo "$eachfile = R1"
			R1=$eachfile
		elif [[ ${eachfile} =~ _R2.fastq$ ]] ; then
			echo "$eachfile = R2"
			R2=$eachfile
		
		elif [[ ${eachfile} =~ _R1.fq$ ]] ; then
			echo "$eachfile = R1"
			R1=$eachfile
		elif [[ ${eachfile} =~ _R2.fq$ ]] ; then
			echo "$eachfile = R2"
			R2=$eachfile
		
		elif [[ ${eachfile} =~ _1.fastq$ ]] ; then
			echo "$eachfile = R1"
			R1=$eachfile
		elif [[ ${eachfile} =~ _2.fastq$ ]] ; then
			echo "$eachfile = R2"
			R2=$eachfile
		
		elif [[ ${eachfile} =~ _1.fq$ ]] ; then
			echo "$eachfile = R1"
			R1=$eachfile
		elif [[ ${eachfile} =~ _2.fq$ ]] ; then
			echo "$eachfile = R2"
			R2=$eachfile
		
		elif [[ ${eachfile} =~ _R1_001.fastp.fastq$ ]] ; then
			echo "$eachfile = R1"
			R1=$eachfile
		elif [[ ${eachfile} =~ _R2_001.fastp.fastq$ ]] ; then
			echo "$eachfile = R2"
			R2=$eachfile
		
		else
			echo "Ignoring $eachfile"
		fi
	done
	#[[ -f $R1 ]] || die "Haven't found R1 in $2"
	#[[ -f $R2 ]] || die "Haven't found R2 in $2"
}

# Folder with FASTQ (for batch analysis) (R1 global variable)
function findFastqInFolder {
	local input=$1
	echo "Checking fastq file..."
	yourfilenames=`ls $1/*`
	for eachfile in $yourfilenames
	do
		if [[ ${eachfile} =~ .fastq$ ]] ; then
			echo "$eachfile = R1"
			R1=$eachfile
		elif [[ ${eachfile} =~ .fq$ ]] ; then
			echo "$eachfile = R1"
			R1=$eachfile
		else
			echo "Ignoring $eachfile"
		fi
	done
	#[[ -f $R1 ]] || die "Haven't found R1 in $2"
}
