#__________________________________________________________
# PROCESSING OF DIRECT RNA SEQUENCING ANALYSIS

ml minimap2
ml samtools
ml dorado

#__________________________________________________________
# BASECALLING

function rna004_basecalling {
    dorado basecaller \
        sup,pseU $1 \
        -r \
        --device $3 >"$2"_basecall.unaligned.bam
}

#__________________________________________________________
# MAPPING 

function rna004_mapping {
    samtools fastq -T '*' "$1" \
	| minimap2 -y --MD -ax splice -uf -k14 -t 96 $2 - \
	| samtools view -bh \
	| samtools sort -@32 - -o ${1/.unaligned.bam/}.$3.bam
	samtools index ${1/.unaligned.bam/}.$3.bam
}

#__________________________________________________________
# MAIN FUNCTION 

function main {
    pod5s=$1
    outdir=$2
    cuda=$3
    ref=$4
    ref_name=$5
    ROI=$6
    rna004_basecalling $pod5s $outdir $cuda
    rna004_mapping "$outdir"_basecall.unaligned.bam $ref $ref_name
    samtools view -hb -F 256 "$outdir"_basecall."$ref_name".bam $ROI | samtools sort -@50 - -o "$outdir"_basecall.ROI.bam
    samtools index "$outdir"_basecall.ROI.bam
}


#___________________________________________________________
# REFERENCES 
oligo_ref="oligos_100.fa"
egfp_mcherry_ref="reporter_sequences.fa"

#___________________________________________________________
# INPUT POD5 FILES
gfppsi_pod5s="./Oligo1/pod5s"
mcherrypsi_pod5s="./Oligo2/pod5s"

sample_A_pods="./sample_A/pod5s"
sample_B_pods="./sample_B/pod5s"
sample_C_pods="./sample_C/pod5s"

#____________________________________________________________
# RUN PREPROCESSING
mkdir -p Results_RNA004_oligos

## Oligo 1
outname="RNA004_Oligo1"
main $gfppsi_pod5s $outdir ./Results_RNA004_oligos/"$outname" "cuda:1,2" $oligo_ref oligos Oligo1:73-115

## Oligo 2
outname="RNA004_Oligo2"
main $gfppsi_pod5s $outdir ./Results_RNA004_oligos/"$outname" "cuda:1,2" $oligo_ref oligos Oligo2:73-115

## Sample A
outname="RNA004_SampleA"
main $sample_A_pods $outdir ./Results_RNA004_oligos/"$outname" "cuda:1,2" $egfp_mcherry_ref vector "EGFP:115 mCherry:565"

## Sample B
outname="RNA004_SampleB"
main $sample_B_pods $outdir ./Results_RNA004_oligos/"$outname" "cuda:1,2" $egfp_mcherry_ref vector "EGFP:115 mCherry:565"

## Sample C
outname="RNA004_SampleC"
main $sample_C_pods $outdir ./Results_RNA004_oligos/"$outname" "cuda:1,2" $egfp_mcherry_ref vector "EGFP:115 mCherry:565"


#____________________________________________________________
# RUN QUALITY CHECK FOR UNALIGNED AND ALIGNED BAM FILES
NanoComp --ubam \
    RNA004_Oligo1_basecall.unaligned.bam RNA004_Oligo2_basecall.unaligned.bam RNA004_SampleA_basecall.unaligned.bam RNA004_SampleB_basecall.unaligned.bam RNA004_SampleC_basecall.unaligned.bam\
    -t 40 \
    -n RNA004_Oligo1 RNA004_Olgio2 RNA004_SampleA RNA004_SampleB RNA004_Control \
    --outdir ./NanoComp_unaligned --raw --store --tsv_stats

NanoComp --bam \
    RNA004_Oligo1_basecall.ROI.bam RNA004_Oligo2_basecall.ROI.bam RNA004_SampleA_basecall.ROI.bam RNA004_SampleB_basecall.ROI.bam RNA004_SampleC_basecall.ROI.bam\
    -t 40 \
    -n RNA004_Oligo1 RNA004_Olgio2 RNA004_SampleA RNA004_SampleB RNA004_Control \
    --outdir ./NanoComp_aligned --raw --store --tsv_stats
