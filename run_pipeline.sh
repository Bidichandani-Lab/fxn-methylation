#!/usr/bin/env bash
amp_values=()
output_dir=""

function help() {
    printf 'Usage: ./%s [OPTIONS] R1 R2 [R1 R2 ...]\n' "$(basename "$0")"
    echo ""
    echo "Description:"
    echo "  This script can do stuff"
    echo ""
    echo "Options:"
    echo "  -a amplicon_number"
    echo "    Amplicon numbers (1-5) that you would like to process (Default: 3)"
    echo ""
    echo "  -o output_directory"
    echo "    Directory where all artifacts and output is saved (Default: ./out/)"
    echo ""
    echo "  -n  selection_count"
    echo "    The number of cells to select for visualizations and statistics (Default: 300)"
    echo ""
    echo "General Options:"
    echo "  -h        Help message"
    echo ""
    echo "Arguments:"
    echo "  R1 R2     Read 1 and read 2 FASTQ(.gz) files"
    echo ""
    echo "Usage Examples:"
    echo "  - Run with all defaults. ONLY AMP3 WILL BE PROCESSED"
    printf '    ./%s Example_L001_R1_001.fastq.gz Example_L001_R2_001.fastq.gz\n' "$(basename "$0")"
    echo ""
    echo "  - Use multiple FASTQ(s). Order doesn't matter as long as each has a pair (R1 and R2). Pairing is made by sorting prefixes before R1/R2."
    printf '    ./%s Example1_L001_R1_001.fastq.gz Example1_L001_R2_001.fastq.gz Example2_L001_R1_001.fastq.gz Example2_L001_R2_001.fastq.gz\n' "$(basename "$0")"
    echo ""
    echo "  - Run using all amplicons (1-5)"
    printf '    ./%s -a1 -a2 -a3 -a4 -a5 Example_L001_R1_001.fastq.gz Example_L001_R2_001.fastq.gz\n' "$(basename "$0")"
}

function display() {
    printf "\033[1;94m%s\033[0m\n" "$1"
}

set -e

while getopts ":a:o:n:h" opt; do
    case $opt in
        a)
            amp_values+=("$OPTARG")
            ;;
        o)
            output_dir="$OPTARG"
            ;;
        n)
            selection_count="$OPTARG"
            ;;
        h)
            help
            exit
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

shift $((OPTIND - 1))
sorted_args=($(printf '%s\n' "$@" | sort))

# Validation
# Make sure output directory exists
if [[ -z $output_dir ]] ; then
    output_dir="out"
    if ! [[ -d "$output_dir" ]] ; then
        mkdir "$output_dir"
    fi
fi
if ! [[ -d "$output_dir" ]] ; then
    echo "output directory \"$output_dir\" does not exist" >&2
    exit 1
fi

# If no amplicons are selected, default to just amplicon 3
if [[ -z $amp_values ]] ; then
    amp_values=(3)
fi
for amp in "${amp_values[@]}" ; do
    if ! [[ "${amp}" =~ ^[1-5]$ ]] ; then
        echo "invalid amp \"${amp}\" must be between 1 and 5" >&2
        exit 1
    fi
    if [[ $amp == 1 ]] ; then
        has_amp1=1
    fi
done

# Validate selection count. All stats, and visualization is based on selection count.
# All else stays the same. Default is 300.
if [[ -z $selection_count ]] ; then
    selection_count=300
fi
if ! [[ $selection_count =~ ^[0-9]+$ ]] ; then
    echo "invalid selection count: $selection_count" >&2
    exit 1
fi

# Make sure each FASTQ is paired (Has R1 and R2). Right now we are just
# sorting and pairing based on lexicographical order.
# A better solution may be to match on a RegEx. (.*)R[12].fastq. That will
# be an exercise for the reader. Make a PR.
for ((i = 0; i < ${#sorted_args[@]}; i+=2)) ; do
    r1_path="${sorted_args[i]}"
    r2_path="${sorted_args[i+1]}"
    if ! [[ -f $r1_path ]] ; then
        echo "read 1 (R1) FASTQ \"$r1_path\" file does not exist " >&2
        exit 1
    fi
    if ! [[ -f $r2_path ]] ; then
        echo "read 2 (R2) FASTQ \"$r2_path\" file does not exist" >&2
        exit 1
    fi
    r1_name=$(echo $(basename -- "$r1_path") | awk -F 'R1' '{print $1}')
    r2_name=$(echo $(basename -- "$r2_path") | awk -F 'R2' '{print $1}')
    if [[ "$r1_name" != "$r2_name" ]] ; then
        echo "R1 and R2 FASTQ files don't match (\"$r1_name\" and \"$r2_name\")" >&2
        exit 1
    fi
done

# samtools tview requires the loci to be zero padded
padded_check_loci_dir=$(mktemp -d -t "check_lociXXXXXXX")
trap "rm -rf $padded_check_loci_dir" EXIT
for file in data/check_loci/* ; do
    awk '{printf "%0*d\n", 3, $1}' "$file" > "$padded_check_loci_dir/$(basename -- "$file")";
done

for ((i = 0; i < ${#sorted_args[@]}; i+=2)) ; do
    run_number=$((i / 2 + 1))
    total_runs=$((${#sorted_args[@]} / 2))

    # Setup directories for this run
    r1_path="${sorted_args[i]}"
    r2_path="${sorted_args[i+1]}"
    name=$(echo $(basename -- "$r1_path") | awk -F 'R1' '{print $1}')
    name="${name%%_}"
    name=$(echo "$name" | tr -d '\n')
    display "processing: ${name} (${run_number} / ${total_runs})"
    out="$output_dir/$name"
    mkdir -p "$out"
    mkdir -p "$out/analysis"

    # Trim adapters
    trimmed_paired_R1="${out}/trimmed_paired_R1.fastq.gz"
    trimmed_paired_R2="${out}/trimmed_paired_R2.fastq.gz"
    java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 "$r1_path" "$r2_path" "${trimmed_paired_R1}" "${out}/trimmed_unpaired_R1_UP.fastq.gz" "${trimmed_paired_R2}" "${out}/trimmed_unpaired_R2_UP.fastq.gz" ILLUMINACLIP:data/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:15 MINLEN:36;

    # Generate bowtie indexes only on the first run
    bt_converted="${output_dir}/bt_amplicons_converted"
    if [[ -z $gen_index ]] ; then
        bowtie2-build data/bowtie_amplicons.fa "${bt_converted}";
        gen_index=1
    fi

    # Merge
    # A different merging script has to be run if amp1 is included because there is gap between amp1 and the rest
    # amp1|...|amp2|amp3|amp4|amp5
    paired_merged="${out}/paired_merged.fastq"
    rm -f "${paired_merged}"
    if [[ -n $has_amp1 ]] ; then
        display "HAS AMP 1"
        fuse.sh in1="${trimmed_paired_R1}" in2="${trimmed_paired_R2}" pad=7 out="${paired_merged}" fusepairs;
    else
        bbmerge.sh -Xmx2048m in1="${trimmed_paired_R1}" in2="${trimmed_paired_R2}" out="${paired_merged}" outu1="${out}/paired_unmerged1.fastq.gz" outu2="${out}/paired_unmerged2.fastq.gz" rsem extend2=50 ecct forcetrimleft=4 mininsert0=180 mininsert=180 pfilter=1;
    fi

    bt2_bam="${out}/bowtie2.bam"
    bowtie2 -p 10 -x "${bt_converted}" -U "${paired_merged}" | samtools sort -T "${out}/sorted.temp" -O bam > "${bt2_bam}";
    samtools index "${bt2_bam}";

    for amp in "${amp_values[@]}"; do
        display "amp: ${amp}"

        bam_amplicon="${out}/bam_amplicon${amp}.bam"
        samtools view -bh "${bt2_bam}" "Amplicon_${amp}" > "${bam_amplicon}";
        samtools sort -n "${bam_amplicon}" -o namesorted.bam;

        bam_random_subset="${out}/random_subset_amplicon${amp}.bam"
        samtools view -h -q 8 namesorted.bam | head -1007 > "${bam_random_subset}";
        rm namesorted.bam;
        bam_sorted_random_subset="${out}/random_subset_amplicon${amp}.sorted.bam"
        samtools sort "${bam_random_subset}" -o "${bam_sorted_random_subset}";
        mv "${bam_sorted_random_subset}" "${bam_random_subset}";
        vcf_random_subset="${out}/random_subset_amplicon${amp}.vcf.gz"
        bcftools mpileup -f data/bowtie_amplicons.fa -x -A -d 10000 -p -t DP -t AD "${bam_random_subset}" | bcftools call -O z -m -T "data/targets/amplicon${amp}.tsv" --ploidy 2 > "${vcf_random_subset}";
        python3 chutake_vcf_counter5m.py "${vcf_random_subset}" "data/confirm_loci/amplicon${amp}.txt" "data/check_loci/amplicon${amp}.txt" "${bam_random_subset}";

        samtools bam2fq "${bam_random_subset}" | seqtk seq -A > "${out}/fasta_amplicon${amp}.fa";
        samtools index "${bam_amplicon}";
        for locus in $(cat "$padded_check_loci_dir/amplicon${amp}.txt") ; do
            samtools tview -d t -p "Amplicon_${amp}:${locus}" "${bam_amplicon}" data/bowtie_amplicons.fa | awk '{FS=""; print $1}' > "temp_${locus}"
        done
        paste -d '\0' temp_* > "${out}/amplicon${amp}.aligned"
        rm temp_*

        python3 visualize.py "${out}/amplicon${amp}.aligned" -n ${selection_count} -a ${amp} -o "$out/analysis"
    done
done
display "COMPLETE"
# ./run_pipeline.sh -a2 -a3 -a4 -a5 -o temp/out/ temp/in/amp2-5/*.fastq.gz