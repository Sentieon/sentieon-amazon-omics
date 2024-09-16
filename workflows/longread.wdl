version 1.0

# Germline variant calling with long reads

workflow sentieon_longread {
  input {
    # Input fastq files
    Array[File] fastq
    Array[String] read_groups

    # Reference genome build
    String reference_name

    # Workflow model files
    File? dnascope_lr_model
    File longreadsv_model

    # Optional process arguments
    String minimap2_xargs = "-L -x map-hifi -Y"
    String sort_xargs = "--cram_write_options version=3.0,compressor=rans"
    String rw_xargs = "--cram_write_options version=3.0,compressor=rans"
    String longreadsv_xargs = ""

    # Sentieon license configuration
    String canonical_user_id
    String sentieon_license = "aws-omics.sentieon.com:9011"

    # Execution
    String n_threads = "32"
    String memory = "64 GiB"
    Int preemptible_tries = 2
    String sentieon_docker
  }
  # Perform a license check
  call SentieonLicense {
    input:
      canonical_user_id = canonical_user_id,
      sentieon_license = sentieon_license,
      sentieon_docker = sentieon_docker,
  }

  # Download the reference genome
  call DownloadReference {
    input:
      reference_name = reference_name,
      sentieon_docker = sentieon_docker,
  }

  call SentieonLongRead {
    input:
      fastq = fastq,
      read_groups = read_groups,

      ref_fasta = DownloadReference.ref_fasta,
      ref_fai = DownloadReference.ref_fai,

      dbsnp_vcf = DownloadReference.dbsnp_vcf,
      dbsnp_vcf_tbi = DownloadReference.dbsnp_vcf_tbi,

      dnascope_lr_model = dnascope_lr_model,
      longreadsv_model = longreadsv_model,

      minimap2_xargs = minimap2_xargs,
      sort_xargs = sort_xargs,
      rw_xargs = rw_xargs,
      longreadsv_xargs = longreadsv_xargs,

      license_ok = SentieonLicense.license_ok,
      canonical_user_id = canonical_user_id,
      sentieon_license = sentieon_license,

      n_threads = n_threads,
      memory = memory,
      preemptible_tries = preemptible_tries,
      sentieon_docker = sentieon_docker,
  }
  output {
    # Alignment files
    File aligned_reads = SentieonLongRead.aligned_reads
    File aligned_index = SentieonLongRead.aligned_index

    # DNAscope LR
    File? calls_vcf = SentieonLongRead.calls_vcf
    File? calls_vcf_tbi = SentieonLongRead.calls_vcf_tbi

    # LongRead SV
    File sv_vcf = SentieonLongRead.sv_vcf
    File sv_vcf_tbi = SentieonLongRead.sv_vcf_tbi
  }
}

task SentieonLicense {
  input {
    String canonical_user_id
    String sentieon_license = "aws-omics.sentieon.com:9011"
    String sentieon_docker
  }
  command <<<
    set -xv
    source /opt/sentieon/omics_credentials.sh "~{sentieon_license}" "~{canonical_user_id}"
    set -exvuo pipefail

    sentieon licclnt ping && echo "Ping is OK"
    sentieon licclnt query Haplotyper
    echo "License OK" >license_ok.txt
    unset http_proxy
  >>>
  runtime {
    preemptible: 3
    docker: sentieon_docker
    memory: "1 GB"
    cpu: 1
  }
  output {
    File license_ok = "license_ok.txt"
  }
}

task DownloadReference {
  input {
    String reference_name
    String sentieon_docker
  }
  command <<<
    set -exvuo pipefail

    if [ -z "${AWS_DEFAULT_REGION-}" ]; then
      set -e
      curl "${ECS_CONTAINER_METADATA_URI_V4}/task" > tmp_task.json
      AWS_DEFAULT_REGION=$(<tmp_task.json jq -rM '.AvailabilityZone' | sed 's/.$//')
      export AWS_DEFAULT_REGION
      rm tmp_task.json
      set +e
    fi

    ref_build=~{reference_name}
    s3_bucket_basename="s3://omics-$AWS_DEFAULT_REGION/reference"
    ref_idx_files=(
      ""
      ".fai"
      ".amb"
      ".ann"
      ".bwt"
      ".pac"
      ".sa"
    )

    case "$ref_build" in
      hg38_alt|hg38_gatk)
        ref_base="hg38/Homo_sapiens_assembly38.fasta"
        has_alt=true
        VCFS=(
          "hg38/dbsnp_138.hg38.vcf.gz"
          "hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"
          "hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        )
        ;;
      hg38|hg38_noalt|hs38)
        ref_base="hg38/hs38.fa"
        has_alt=false
        VCFS=(
          "hg38/dbsnp_138.hg38.vcf.gz"
          "hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"
          "hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        )
        ;;
      b37_gatk)
        ref_base="b37/Homo_sapiens_assembly19.fasta"
        has_alt=false
        VCFS=(
          "b37/dbsnp_138.b37.vcf.gz"
          "b37/1000G_phase1.indels.b37.vcf.gz"
          "b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
        )
        ;;
      b37|hs37d5)
        ref_base="b37/hs37d5.fa"
        has_alt=false
        VCFS=(
          "b37/dbsnp_138.b37.vcf.gz"
          "b37/1000G_phase1.indels.b37.vcf.gz"
          "b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
        )
        ;;
      hg19|ucsc_hg19)
        ref_base="hg19/ucsc.hg19.fasta"
        has_alt=false
        VCFS=(
          "hg19/dbsnp_138.hg19.vcf.gz"
          "hg19/1000G_phase1.indels.hg19.sites.vcf.gz"
          "hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
        )
        ;;
      quickstart)
        ref_base="quickstart/ucsc.hg19_chr22.fasta"
        has_alt=false
        VCFS=(
          "quickstart/dbsnp_135.hg19_chr22.vcf.gz"
          "quickstart/1000G_phase1.snps.high_confidence.hg19_chr22.sites.vcf.gz"
          "quickstart/Mills_and_1000G_gold_standard.indels.hg19_chr22.sites.vcf.gz"
        )
        ;;
      *)
        echo "ERROR: unknown genome build"
        exit 1
        ;;
    esac

    for idx in "${ref_idx_files[@]}"; do
      aws s3 cp "$s3_bucket_basename/$ref_base$idx" "./reference.fa$idx"
    done
    if [[ "$has_alt" == true ]]; then
      aws s3 cp "$s3_bucket_basename/${ref_base}.alt" "./reference.fa.alt"
    fi
    for i in $(seq 1 "${#VCFS[@]}"); do
      i=$((i - 1))
      vcf="${VCFS[$i]}"
      aws s3 cp "$s3_bucket_basename/$vcf" "./sites_${i}.vcf.gz"
      aws s3 cp "$s3_bucket_basename/${vcf}.tbi" "./sites_${i}.vcf.gz.tbi"
    done
  >>>
  runtime {
    preemptible: 3
    docker: sentieon_docker
    memory: "4 GB"
    cpu: 1
  }
  output {
    File ref_fasta = "reference.fa"
    File ref_fai = "reference.fa.fai"
    File? ref_alt = "reference.fa.alt"
    File ref_bwt = "reference.fa.bwt"
    File ref_sa = "reference.fa.sa"
    File ref_amb = "reference.fa.amb"
    File ref_ann = "reference.fa.ann"
    File ref_pac = "reference.fa.pac"

    File dbsnp_vcf = "sites_0.vcf.gz"
    File dbsnp_vcf_tbi = "sites_0.vcf.gz.tbi"
    Array[File] bqsr_vcfs = glob("sites_*.vcf.gz")
    Array[File] bqsr_vcf_tbis = glob("sites_*.vcf.gz.tbi")
  }
}

task SentieonLongRead {
  input {
    # Input fastq files
    Array[File] fastq
    Array[String] read_groups

    # Reference genome files
    File ref_fasta
    File ref_fai

    # Known sites VCFs
    File? dbsnp_vcf
    File? dbsnp_vcf_tbi

    # Workflow arguments
    File? junc_bed

    # Workflow model files
    File? dnascope_lr_model
    File longreadsv_model

    # Optional process arguments
    String minimap2_xargs = "-L -x map-hifi -Y"
    String sort_xargs = "--cram_write_options version=3.0,compressor=rans"
    String rw_xargs = "--cram_write_options version=3.0,compressor=rans"
    String longreadsv_xargs = ""

    # Sentieon license configuration
    File license_ok
    String canonical_user_id
    String sentieon_license = "aws-omics.sentieon.com:9011"

    # Execution
    String n_threads = "32"
    String memory = "64 GiB"
    Int preemptible_tries = 2
    String sentieon_docker
  }
  command <<<
    set -xv
    source /opt/sentieon/omics_credentials.sh "~{sentieon_license}" "~{canonical_user_id}"
    set -exvuo pipefail

    # Configuration
    nt=$(nproc)

    fastq=("~{sep='" "' fastq}")
    read_groups=("~{sep='" "' read_groups}")
    junc_bed="~{default='' junc_bed}"

    # Sanity check the input
    if [[ ${#fastq[@]} -ne ${#read_groups[@]} ]]; then
      echo "The number of readgroups does not match the number of input fastq files"
      exit 1
    fi

    # Set the BQSR and calling intervals
    first_chrom=$(head -n 1 ~{ref_fai} | cut -f 1)
    set +e +o pipefail
    case "$first_chrom" in
      chrM)
        cat "~{ref_fai}" | head -n 23 | tail -n +2 | awk -v 'OFS=\t' '{print $1,0,$2}' > diploid_intervals.bed
        CALLING_INTERVALS="chrM,chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
        ;;
      1)
        cat "~{ref_fai}" | head -n 22 | awk -v 'OFS=\t' '{print $1,0,$2}' > diploid_intervals.bed
        CALLING_INTERVALS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT"
        ;;
      chr1)
        cat "~{ref_fai}" | head -n 22 | awk -v 'OFS=\t' '{print $1,0,$2}' > diploid_intervals.bed
        CALLING_INTERVALS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
        ;;
      chr22)
        echo "chr22	0	51304566" > diploid_intervals.bed
        CALLING_INTERVALS="chr22"
        ;;
      *)
        echo "ERROR: unknown first reference chromosome"
        exit 1
        ;;
    esac
    set -exvuo pipefail

    # Alignment
    alignment_output=()
    for i in $(seq 1 ${#fastq[@]}); do
      i=$((i - 1))
      fq="${fastq[$i]}"
      rg="${read_groups[$i]}"
      sentieon minimap2 -R "$rg" ${junc_bed:+--junc-bed "$junc_bed"} \
        ~{minimap2_xargs} -t $nt -a "~{ref_fasta}" "$fq" | \
        sentieon util sort -t $nt --sam2bam -o "sample_${i}_sorted.cram" \
        -i - ~{sort_xargs} --reference "~{ref_fasta}"
      alignment_output+=("-i" "sample_${i}_sorted.cram")

      (rm "$fq" || (exit 0)) &  # Try removing the input files to save disk space
    done

    # Merge multiple input fastq
    sentieon driver "${alignment_output[@]}" -r "~{ref_fasta}" --algo ReadWriter \
      ~{rw_xargs} "sample_aligned.cram"
    (rm "${alignment_output[@]}" || (exit 0))

    # DNAscope LR
    dnascope_lr_model="~{default='' dnascope_lr_model}"
    if [[ -n "$dnascope_lr_model" ]]; then
      # Unpack the tar file
      mkdir pipeline
      tar -C pipeline -zxf "$dnascope_lr_model"
      dnascope_lr=$(ls ./pipeline/*/dnascope_HiFi.sh)
      dnascope_lr_model=$(ls ./pipeline/*/*.model)

      # Variant calling
      dbsnp_vcf="~{default='' dbsnp_vcf}"
      bash "$dnascope_lr" ${dbsnp_vcf:+-d "$dbsnp_vcf"} \
        -m "$dnascope_lr_model" -b "diploid_intervals.bed" -t "$nt" \
        -r "~{ref_fasta}" -i "sample_aligned.cram" "sample_calls.vcf.gz"
      rm -r pipeline
    fi

    # LongReadSV
    sentieon driver -r "~{ref_fasta}" -i "sample_aligned.cram" \
      --interval "$CALLING_INTERVALS" \
      -t $nt --algo LongReadSV --model "~{longreadsv_model}" \
      ~{longreadsv_xargs} "sample_svs.vcf.gz"
    unset http_proxy
    exit 0
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: sentieon_docker
    memory: memory
    cpu: n_threads
  }
  output {
    # Alignment files
    File aligned_reads = "sample_aligned.cram"
    File aligned_index = "sample_aligned.cram.crai"

    # DNAscope LR
    File? calls_vcf = "sample_calls.vcf.gz"
    File? calls_vcf_tbi = "sample_calls.vcf.gz.tbi"

    # LongRead SV
    File sv_vcf = "sample_svs.vcf.gz"
    File sv_vcf_tbi = "sample_svs.vcf.gz.tbi"
  }
}
