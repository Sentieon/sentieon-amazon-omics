version 1.0

# Somatic variant calling with Sentieon TNhaplotyper2

workflow sentieon_somatic {
  input {
    # Input tumor fastq files
    Array[File] r1_fastq = []
    Array[File] r2_fastq = []
    Array[String] read_groups = []

    # Input normal fastq files
    Array[File] normal_r1_fastq = []
    Array[File] normal_r2_fastq = []
    Array[String] normal_read_groups = []
    
    # Reference genome files
    String reference_name

    # Known sites VCFs
    File? pon_vcf
    File? pon_vcf_tbi
    File? germline_vcf
    File? germline_vcf_tbi
    File contamination_vcf
    File contamination_vcf_tbi

    # Optional process arguments
    String bwa_xargs = ""
    String bwa_karg = "10000000"
    String sort_xargs = "--bam_compression 1"
    String lc_xargs = ""
    String alnstat_xargs = "--adapter_seq ''"
    String dedup_xargs = "--cram_write_options version=3.0,compressor=rans"
    String qcal_xargs = ""
    String calling_driver_xargs = ""
    String calling_algo_xargs = ""

    # Sentieon license configuration
    String canonical_user_id
    String sentieon_license = "aws-omics.sentieon.com:9011"

    # Execution
    String n_threads = "32"
    String memory = "64 GiB"
    Int preemptible_tries = 3
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

  call SentieonSomatic {
    input:
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      read_groups = read_groups,

      normal_r1_fastq = normal_r1_fastq,
      normal_r2_fastq = normal_r2_fastq,
      normal_read_groups = normal_read_groups,

      ref_fasta = DownloadReference.ref_fasta,
      ref_fai = DownloadReference.ref_fai,
      ref_alt = DownloadReference.ref_alt,
      ref_bwt = DownloadReference.ref_bwt,
      ref_sa = DownloadReference.ref_sa,
      ref_amb = DownloadReference.ref_amb,
      ref_ann = DownloadReference.ref_ann,
      ref_pac = DownloadReference.ref_pac,

      bqsr_vcfs = DownloadReference.bqsr_vcfs,
      bqsr_vcf_tbis = DownloadReference.bqsr_vcf_tbis,
      pon_vcf = pon_vcf,
      pon_vcf_tbi = pon_vcf_tbi,
      germline_vcf = germline_vcf,
      germline_vcf_tbi = germline_vcf_tbi,
      contamination_vcf = contamination_vcf,
      contamination_vcf_tbi = contamination_vcf_tbi,

      bwa_xargs = bwa_xargs,
      bwa_karg = bwa_karg,
      sort_xargs = sort_xargs,
      lc_xargs = lc_xargs,
      alnstat_xargs = alnstat_xargs,
      dedup_xargs = dedup_xargs,
      qcal_xargs = qcal_xargs,
      calling_driver_xargs = calling_driver_xargs,
      calling_algo_xargs = calling_algo_xargs,

      license_ok = SentieonLicense.license_ok,
      canonical_user_id = canonical_user_id,
      sentieon_license = sentieon_license,

      n_threads = n_threads,
      memory = memory,
      preemptible_tries = preemptible_tries,
      sentieon_docker = sentieon_docker,
  }
  output {
    # Core alignment files
    File aligned_reads = SentieonSomatic.aligned_reads
    File aligned_index = SentieonSomatic.aligned_index
    File? normal_reads = SentieonSomatic.normal_reads
    File? normal_index = SentieonSomatic.normal_index

    # QualCal output
    File? bqsr_table = SentieonSomatic.bqsr_table
    File? normal_bqsr_table = SentieonSomatic.normal_bqsr_table

    # VCF outputs
    File calls_vcf = SentieonSomatic.calls_vcf
    File calls_vcf_tbi = SentieonSomatic.calls_vcf_tbi

    # QC output metrics
    File dedup_metrics = SentieonSomatic.dedup_metrics
    File mq_metrics = SentieonSomatic.mq_metrics
    File qd_metrics = SentieonSomatic.qd_metrics
    File gc_summary = SentieonSomatic.gc_summary
    File gc_metrics = SentieonSomatic.gc_metrics
    File as_metrics = SentieonSomatic.as_metrics
    File is_metrics = SentieonSomatic.is_metrics

    # QC output plots
    File mq_plot = SentieonSomatic.mq_plot
    File qd_plot = SentieonSomatic.qd_plot
    File gc_plot = SentieonSomatic.gc_plot
    File is_plot = SentieonSomatic.is_plot

    # Normal QC output metrics
    File? normal_dedup_metrics = SentieonSomatic.normal_dedup_metrics
    File? normal_mq_metrics = SentieonSomatic.normal_mq_metrics
    File? normal_qd_metrics = SentieonSomatic.normal_qd_metrics
    File? normal_gc_summary = SentieonSomatic.normal_gc_summary
    File? normal_gc_metrics = SentieonSomatic.normal_gc_metrics
    File? normal_as_metrics = SentieonSomatic.normal_as_metrics
    File? normal_is_metrics = SentieonSomatic.normal_is_metrics

    # Normal QC output plots
    File? normal_mq_plot = SentieonSomatic.normal_mq_plot
    File? normal_qd_plot = SentieonSomatic.normal_qd_plot
    File? normal_gc_plot = SentieonSomatic.normal_gc_plot
    File? normal_is_plot = SentieonSomatic.normal_is_plot
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

task SentieonSomatic {
  input {
    # Input tumor fastq files
    Array[File] r1_fastq = []
    Array[File] r2_fastq = []
    Array[String] read_groups = []

    # Input normal fastq files
    Array[File] normal_r1_fastq = []
    Array[File] normal_r2_fastq = []
    Array[String] normal_read_groups = []
    
    # Reference genome files
    File ref_fasta
    File ref_fai
    File? ref_alt
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac

    # Known sites VCFs
    Array[File] bqsr_vcfs
    Array[File] bqsr_vcf_tbis
    File? pon_vcf
    File? pon_vcf_tbi
    File? germline_vcf
    File? germline_vcf_tbi
    File contamination_vcf
    File contamination_vcf_tbi

    # Optional process arguments
    String bwa_xargs = ""
    String bwa_karg = "10000000"
    String sort_xargs = "--bam_compression 1"
    String lc_xargs = ""
    String alnstat_xargs = "--adapter_seq ''"
    String dedup_xargs = "--cram_write_options version=3.0,compressor=rans"
    String qcal_xargs = ""
    String calling_driver_xargs = ""
    String calling_algo_xargs = ""

    # Sentieon license configuration
    File license_ok
    String canonical_user_id
    String sentieon_license = "aws-omics.sentieon.com:9011"

    # Execution
    String n_threads = "32"
    String memory = "64 GiB"
    Int preemptible_tries = 3
    String sentieon_docker

    # Expressions
    Boolean has_normal = length(normal_r1_fastq) > 0
  }
  command <<<
    set -xv
    source /opt/sentieon/omics_credentials.sh "~{sentieon_license}" "~{canonical_user_id}"
    set -exvuo pipefail

    # Get the NUMA configuration
    numa_nodes=$(lscpu | grep "NUMA node(s):" | sed 's/^NUMA node.* //')
    numa_cpulist=()
    for i in $(seq 1 "$numa_nodes"); do
        i=$((i - 1))
        numa_cpulist+=($(lscpu | grep "NUMA node$i CPU" | sed 's/^NUMA.* //'))
    done

    # Configuration
    nt=$(nproc)
    n_threads=$((nt / numa_nodes))

    mem_kb=$(cat /proc/meminfo | grep "MemTotal" | awk '{print $2}')
    bwt_mem=$((mem_kb / 1024 / 1024 / numa_nodes - 10))
    export bwt_max_mem="$bwt_mem"G

    r1_fastq=("~{sep='" "' r1_fastq}")
    r2_fastq=("~{sep='" "' r2_fastq}")
    read_groups=("~{sep='" "' read_groups}")

    normal_r1_fastq=()
    normal_r2_fastq=()
    normal_read_groups=()
    has_normal=~{true="true" false="" has_normal}
    if [[ "$has_normal" == "true" ]]; then
      normal_r1_fastq=("~{sep='" "' normal_r1_fastq}")
      normal_r2_fastq=("~{sep='" "' normal_r2_fastq}")
      normal_read_groups=("~{sep='" "' normal_read_groups}")
    fi

    # Sanity check the input
    if [[ ${#r1_fastq[@]} -ne ${#read_groups[@]} ]]; then
      echo "The number of readgroups does not match the number of input fastq files"
      exit 1
    fi
    if [[ ${#r1_fastq[@]} -ne ${#r2_fastq[@]} ]]; then
      echo "The number of r1 fastq does not equal the number of r2 fastq"
      exit 1
    fi
    if [[ ${#normal_r1_fastq[@]} -ne ${#normal_read_groups[@]} ]]; then
      echo "The number of normal readgroups does not match the number of input fastq files"
      exit 1
    fi
    if [[ ${#normal_r1_fastq[@]} -ne ${#normal_r2_fastq[@]} ]]; then
      echo "The number of normal r1 fastq does not equal the number of r2 fastq"
      exit 1
    fi

    # Set the BQSR and calling intervals
    first_chrom=$(head -n 1 ~{ref_fai} | cut -f 1)
    case "$first_chrom" in
      chrM)
        BQSR_INTERVALS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"
        CALLING_INTERVALS="chrM,chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
        ;;
      1)
        BQSR_INTERVALS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"
        CALLING_INTERVALS="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT"
        ;;
      chr1)
        BQSR_INTERVALS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"
        CALLING_INTERVALS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"
        ;;
      chr22)
        BQSR_INTERVALS="chr22"
        CALLING_INTERVALS="chr22"
        ;;
      *)
        echo "ERROR: unknown first reference chromosome"
        exit 1
        ;;
    esac

    # Ensure VCFs are adjecent to their index
    vcf_files=("~{sep='" "' bqsr_vcfs}")
    tbi_files=("~{sep='" "' bqsr_vcf_tbis}")
    bqsr_args=()
    for i in $(seq 1 ${#vcf_files[@]}); do
      i=$((i - 1))
      ln -s "${vcf_files[$i]}" ./input_vcf_"$i".vcf.gz
      ln -s "${tbi_files[$i]}" ./input_vcf_"$i".vcf.gz.tbi
      bqsr_args+=("-k" "./input_vcf_${i}.vcf.gz")
    done

    pon_vcf="~{default='' pon_vcf}"
    if [[ -n "$pon_vcf" ]]; then
      ln -s "$pon_vcf" ./input_pon.vcf.gz
      ln -s "~{default='' pon_vcf_tbi}" ./input_pon.vcf.gz.tbi
    fi

    germline_vcf="~{default='' germline_vcf}"
    if [[ -n "$germline_vcf" ]]; then
      ln -s "$germline_vcf" ./input_germine.vcf.gz
      ln -s "~{default='' germline_vcf_tbi}" ./input_germine.vcf.gz.tbi
    fi

    contamination_vcf="~{default='' contamination_vcf}"
    if [[ -n "$contamination_vcf" ]]; then
      ln -s "$contamination_vcf" ./input_contamination.vcf.gz
      ln -s "~{default='' contamination_vcf_tbi}" ./input_contamination.vcf.gz.tbi
    fi

    # Alignment with BWA
    alignment_output=()
    for i in $(seq 1 ${#r1_fastq[@]}); do
      i=$((i - 1))
      r1="${r1_fastq[$i]}"
      r2="${r2_fastq[$i]}"
      rg="${read_groups[$i]}"

      for j in $(seq 1 "$numa_nodes"); do
        j=$((j - 1))
        cpulist="${numa_cpulist[$j]}"

        # Alignment command
        perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)';
        taskset -c "$cpulist" sentieon bwa mem -R "$rg" \
          ~{bwa_xargs} -K ~{bwa_karg} -t $n_threads -p "~{ref_fasta}" \
          <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
            sentieon fqidx extract -F "$j"/"$numa_nodes" -K ~{bwa_karg} \
            <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
              igzip -dc "$r1") \
            <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
              igzip -dc "$r2")) | \
          taskset -c "$cpulist" sentieon util sort -t $n_threads --sam2bam \
          -o "sample_sorted_${i}_${j}.bam" -i - ~{sort_xargs} &
        alignment_output+=("sample_sorted_${i}_${j}.bam")
      done
      wait

      (rm "$r1" "$r2" || (exit 0)) &  # Try removing the input files to save disk space
    done

    # Alignment with BWA - normal
    normal_alignment_output=()
    for i in $(seq 1 ${#normal_r1_fastq[@]}); do
      i=$((i - 1))
      r1="${normal_r1_fastq[$i]}"
      r2="${normal_r2_fastq[$i]}"
      rg="${normal_read_groups[$i]}"

      for j in $(seq 1 "$numa_nodes"); do
        j=$((j - 1))
        cpulist="${numa_cpulist[$j]}"

        # Alignment command
        perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)';
        taskset -c "$cpulist" sentieon bwa mem -R "$rg" \
          ~{bwa_xargs} -K ~{bwa_karg} -t $n_threads -p "~{ref_fasta}" \
          <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
            sentieon fqidx extract -F "$j"/"$numa_nodes" -K ~{bwa_karg} \
            <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
              igzip -dc "$r1") \
            <(perl -MFcntl -e 'fcntl(STDOUT, 1031, 268435456)'; \
              igzip -dc "$r2")) | \
          taskset -c "$cpulist" sentieon util sort -t $n_threads --sam2bam \
          -o "normal_sorted_${i}_${j}.bam" -i - ~{sort_xargs} &
        normal_alignment_output+=("normal_sorted_${i}_${j}.bam")
      done
      wait

      (rm "$r1" "$r2" || (exit 0)) &  # Try removing the input files to save disk space
    done

    # Add the bam input files to the bam_str
    bam_str=()
    for i in $(seq 1 ${#alignment_output[@]}); do
      i=$((i - 1))
      bam_str+=("-i")
      bam_str+=("${alignment_output[$i]}")
    done
    normal_bam=()
    for i in $(seq 1 ${#normal_alignment_output[@]}); do
      i=$((i - 1))
      normal_bam+=("-i")
      normal_bam+=("${normal_alignment_output[$i]}")
    done

    # Dedup - tumor
    tumor_deduped="sample_aligned.cram"
    sentieon driver "${bam_str[@]}" -r ~{ref_fasta} \
      --algo LocusCollector ~{lc_xargs} "sample_score.txt.gz" \
      --algo MeanQualityByCycle "sample_mq_metrics.txt" \
      --algo QualDistribution "sample_qd_metrics.txt" \
      --algo GCBias --summary "sample_gc_summary.txt" "sample_gc_metrics.txt" \
      --algo AlignmentStat ~{alnstat_xargs} "sample_aln_metrics.txt" \
      --algo InsertSizeMetricAlgo "sample_is_metrics.txt"

    sentieon driver "${bam_str[@]}" -r ~{ref_fasta} --algo Dedup \
      ~{dedup_xargs} --score_info "sample_score.txt.gz" \
      --metrics "sample_dedup_metrics.txt" \
      "$tumor_deduped"

    # Plot the metrics output
    sentieon plot GCBias -o "sample_gc-report.pdf" "sample_gc_metrics.txt" &
    sentieon plot QualDistribution -o "sample_qd-report.pdf" "sample_qd_metrics.txt" &
    sentieon plot MeanQualityByCycle -o "sample_mq-report.pdf" "sample_mq_metrics.txt" &
    sentieon plot InsertSizeMetricAlgo -o "sample_is-report.pdf" "sample_is_metrics.txt" &

    (rm "${alignment_output[@]}" || (exit 0)) &   # Remove intermediate files to save space

    # Dedup - normal
    normal_deduped=""
    if [[ "$has_normal" == "true" ]]; then
      normal_deduped="normal_aligned.cram"
      sentieon driver "${normal_bam[@]}" -r ~{ref_fasta} \
        --algo LocusCollector ~{lc_xargs} "normal_score.txt.gz" \
        --algo MeanQualityByCycle "normal_mq_metrics.txt" \
        --algo QualDistribution "normal_qd_metrics.txt" \
        --algo GCBias --summary "normal_gc_summary.txt" "normal_gc_metrics.txt" \
        --algo AlignmentStat ~{alnstat_xargs} "normal_aln_metrics.txt" \
        --algo InsertSizeMetricAlgo "normal_is_metrics.txt"

      sentieon driver "${normal_bam[@]}" -r ~{ref_fasta} --algo Dedup \
        ~{dedup_xargs} --score_info "normal_score.txt.gz" \
        --metrics "normal_dedup_metrics.txt" \
        "$normal_deduped"

      # Plot the metrics output
      sentieon plot GCBias -o "normal_gc-report.pdf" "normal_gc_metrics.txt" &
      sentieon plot QualDistribution -o "normal_qd-report.pdf" "normal_qd_metrics.txt" &
      sentieon plot MeanQualityByCycle -o "normal_mq-report.pdf" "normal_mq_metrics.txt" &
      sentieon plot InsertSizeMetricAlgo -o "normal_is-report.pdf" "normal_is_metrics.txt" &

      (rm "${normal_alignment_output[@]}" || (exit 0)) &   # Remove intermediate files to save space
    fi

    # BQSR
    bqsr_out="sample_recal.table"
    normal_bqsr_out=""
    sentieon driver -r ~{ref_fasta} -i "$tumor_deduped" \
      --interval "$BQSR_INTERVALS" --algo QualCal \
      "${bqsr_args[@]}" \
      ~{qcal_xargs} "$bqsr_out"

    if [[ "$has_normal" == "true" ]]; then
      # BQSR - Normal
      normal_bqsr_out="normal_recal.table"
      sentieon driver -r ~{ref_fasta} -i "$normal_deduped" \
        --interval "$BQSR_INTERVALS" --algo QualCal \
        "${bqsr_args[@]}" \
        ~{qcal_xargs} "$normal_bqsr_out"
    fi

    # Extract the samples
    tumor_sm=$(samtools samples "$tumor_deduped" | cut -f 1 | head -n 1)
    normal_sm=""
    if [[ "$has_normal" == "true" ]]; then
      normal_sm=$(samtools samples "$normal_deduped" | cut -f 1 | head -n 1)
    fi

    # Variant calling
    sentieon driver -r "~{ref_fasta}" \
      -i "$tumor_deduped" \
      ${normal_deduped:+-i "$normal_deduped"} \
      ${bqsr_out:+-q "$bqsr_out"} \
      ${normal_bqsr_out:+-q "$normal_bqsr_out"} \
      --interval "$CALLING_INTERVALS" \
      --algo TNhaplotyper2 \
        --tumor_sample "$tumor_sm" \
        ${normal_sm:+--normal_sample "$normal_sm"} \
        ${pon_vcf:+--pon "./input_pon.vcf.gz"} \
        ${germline_vcf:+--germline_vcf "./input_germine.vcf.gz"} \
        "sample_tnhap2_tmp.vcf.gz" \
      --algo OrientationBias \
        --tumor_sample "$tumor_sm" \
        "sample_orientation" \
      --algo ContaminationModel \
        --tumor_sample "$tumor_sm" \
        ${normal_sm:+--normal_sample "$normal_sm"} \
        ${contamination_vcf:+--vcf "./input_contamination.vcf.gz"} \
        --tumor_segments "sample_segments" \
        "sample_contamination"

    # TNfilter
    sentieon driver -r "~{ref_fasta}" \
      --algo TNfilter \
      -v "sample_tnhap2_tmp.vcf.gz" \
      --tumor_sample "$tumor_sm" \
      ${normal_sm:+--normal_sample "$normal_sm"} \
      --contamination "sample_contamination" \
      --tumor_segments "sample_segments" \
      --orientation_priors "sample_orientation" \
      "sample_tnhap2.vcf.gz"
    wait
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
    # Core alignment files
    File aligned_reads = "sample_aligned.cram"
    File aligned_index = "sample_aligned.cram.crai"
    File? normal_reads = "normal_aligned.cram"
    File? normal_index = "normal_aligned.cram.crai"

    # QualCal output
    File? bqsr_table = "sample_recal.table"
    File? normal_bqsr_table = "normal_recal.table"

    # VCF outputs
    File calls_vcf = "sample_tnhap2.vcf.gz"
    File calls_vcf_tbi = "sample_tnhap2.vcf.gz.tbi"

    # QC output metrics
    File dedup_metrics = "sample_dedup_metrics.txt"
    File mq_metrics = "sample_mq_metrics.txt"
    File qd_metrics = "sample_qd_metrics.txt"
    File gc_summary = "sample_gc_summary.txt"
    File gc_metrics = "sample_gc_metrics.txt"
    File as_metrics = "sample_aln_metrics.txt"
    File is_metrics = "sample_is_metrics.txt"

    # QC output plots
    File mq_plot = "sample_mq-report.pdf"
    File qd_plot = "sample_qd-report.pdf"
    File gc_plot = "sample_gc-report.pdf"
    File is_plot = "sample_is-report.pdf"

    # Normal QC output metrics
    File? normal_dedup_metrics = "normal_dedup_metrics.txt"
    File? normal_mq_metrics = "normal_mq_metrics.txt"
    File? normal_qd_metrics = "normal_qd_metrics.txt"
    File? normal_gc_summary = "normal_gc_summary.txt"
    File? normal_gc_metrics = "normal_gc_metrics.txt"
    File? normal_as_metrics = "normal_aln_metrics.txt"
    File? normal_is_metrics = "normal_is_metrics.txt"

    # Normal QC output plots
    File? normal_mq_plot = "normal_mq-report.pdf"
    File? normal_qd_plot = "normal_qd-report.pdf"
    File? normal_gc_plot = "normal_gc-report.pdf"
    File? normal_is_plot = "normal_is-report.pdf"
  }
}
