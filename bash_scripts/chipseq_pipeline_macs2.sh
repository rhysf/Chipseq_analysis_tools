#!/bin/bash -l
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source $DIR/Bash_functions/File_input_output.sh

### rfarrer@broadinstitute.org

# Exit with message
diemsg() {
	echo "Usage: $0 -r <treated.BAM> -c <control.BAM>"
	echo ""
	echo "Optional: -e Experiment tag [treated_vs_untreated_MACS]"
	echo "          -i IDR directory [/Users/rhysfarrer/Desktop/Programs/idr-2.0.2/]"
	echo "          -n Effective genome size (Gemtools) (Default CNA2 with 150nt reads) [18627656]"
	echo "          -q Qvalue (For broad marks, you can try 0.05) [0.05]"
	echo "          -b Broad peaks (y/n) [y]"
	echo "          -a Cutoff-analysis (y/n) [n]"
	echo "          -m Min MFOLD range [5]"
	echo "          -x Max MFOLD range [50]"
	echo "          -t ExtSize (The arbitrary extension size in bp. When nomodel is true, MACS will use this value as fragment size to extend each read towards 3' end) [200]"
	echo "          -o NoModel (Whether or not to build the shifting model) [n]"
	echo ""
	echo "Notes:    >All replicates should be in current working directory, needed for merging files"
	echo "          >Prerequisites: Bedtools, Samtools and MACS3 are in PATH"
	echo "          >IDR installed from https://github.com/nboley/idr"
	echo "          >Cutoff-analysis will analyze number or total length of peaks that can be called by different p-value cutoff then output a summary table to help"
	echo "          user decide a better cutoff. The table will be saved in NAME_cutoff_analysis.txt file"
	echo "          >Effective genome size calculated with perl script. shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8)"
	echo "          >MFOLD range Tweaking it is not recommended. However, Too few paired peaks (0) so I can not build the model! Broader your MFOLD range parameter may erase this error. If it still can't build the model, we suggest to use --nomodel and --extsize 147 or other fixed number instead."
	echo "          >NoModel suggestion use both 'nomodel' and 'extsize' to make signal files from different datasets comparable."
	die  "Outputs: ./consistency/ ./mapped/ ./peaks/ and ./stats/" 
}
[ $# -gt 3 ] || diemsg

exists() {
	command -v "$1" >/dev/null 2>&1
}

# Gather input files and settings
declare -A opts
opts[a]="n"
opts[b]="y"
opts[e]="treated_vs_untreated_MACS"
opts[i]="/Users/rhysfarrer/Desktop/Programs/idr-2.0.2/"
opts[m]=5
opts[n]=18627656
opts[o]="n"
opts[q]="0.05"
opts[t]=200
opts[x]=50
getOpts $@

# Check input files are readable
[ -r ${opts[r]} ] || die "File: ${opts[r]} does not appear to be valid"
[ -r ${opts[c]} ] || die "File: ${opts[c]} does not appear to be valid"

# Programs required
#[ -x bamToBed ] || die "bamToBed in bedtools/bin not found. Add to PATH"
#[ -x samtools ] || die "samtools not found. Add to PATH"
#[ -x macs2 ] || die "macs2 not found. Add to PATH"
if exists bash; then
	echo "macs2 found in PATH"
else
	die "macs2 not found in PATH"
fi

# Broad peaks
SETTINGS=
[ ${opts[b]} == "n" ] || SETTINGS="--broad"
[ ${opts[a]} == "n" ] || SETTINGS="${SETTINGS} --cutoff-analysis"
[ ${opts[o]} == "n" ] || SETTINGS="${SETTINGS} --nomodel"

# make output directory
CMD0="mkdir ${opts[e]}"
echo "CMD0 : $CMD0"
eval $CMD0

# Call Peaks using MACS2 on individual replicates (using relaxed thresholds)
CMD1="macs2 callpeak --treatment ${opts[r]} --control ${opts[c]} -f BAM -g ${opts[n]} -n ${opts[e]} --outdir ${opts[e]} -q ${opts[q]} -m ${opts[m]} ${opts[x]} --extsize ${opts[t]} ${SETTINGS}"
echo "CMD1 : $CMD1"
echo "CMD1 : $CMD1" > ${opts[e]}/commands.txt 2>&1 
eval ${CMD1} >> ${opts[e]}/commands.txt 2>&1 

