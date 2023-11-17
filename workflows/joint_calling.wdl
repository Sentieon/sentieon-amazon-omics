version 1.0

# Sentieon's pipeline for distributed joint calling

workflow distributed_gvcftyper {
  input {
    # Input gVCFs
    #Array[File] gvcfs
    File gvcf_list  # Input gVCF files; one file per line

    # Reference genome files
    File ref_fasta
    File ref_fai

    # Settings for distribution
    String? region  # Can be used to limit the joint call to a specific region. Defaults to the whole genome
    Int compute_workers = 4  # Acceptable for small to medium joint calls. May need to be increased for larger joint calls
    Int shard_size = 100000000  # The shard size, decrease for large joint calls
    String gvcftyper_memory = "32 GiB"
    String gvcftyper_threads = "16"
    Int concurrent_downloads = 5

    # Settings for the merged output
    Int output_splits = 5  # Number of VCFs to split the output into
    String merge_memory = "32 GiB"
    String merge_threads = "16"

    # GVCFtyper arguments
    File? dbsnp_vcf
    File? dbsnp_vcf_tbi
    String driver_xargs = "--traverse_param 10000/200"
    String gvcftyper_xargs = "--genotype_model multinomial --max_alt_alleles 12"

    # Sentieon license configuration
    String canonical_user_id
    String sentieon_license = "aws-omics.sentieon.com:9011"

    # Execution
    String sentieon_docker
  }

  # GenerateShards output:
  #  GenerateShards.gvcftyper_shards - [[shard1, shard2, shard3], [shard4, shard5, ...], ...]
  #  GenerateShards.merge_start - [0, 5, 10, 15, ...] - start index of shards for the merge
  #  GenerateShards.merge_size - number of shards for the merge
  #  GenerateShards.merge_shards = [merge_shard1, merge_shard2, ... ]
  call GenerateShards {
    input:
      ref_fai = ref_fai,
      region = region,
      compute_workers = compute_workers,
      shard_size = shard_size,
      output_splits = output_splits,
      sentieon_docker = sentieon_docker,
  }

  scatter(i in range(length(GenerateShards.gvcftyper_shards))) {
    # DistributedGvcfyper output:
    #  DistributedGvcfyper.shard_vcfs - [vcf1, vcf2, ...]
    #  DistributedGvcfyper.shard_tbis - [tbi1, tbi2, ...]
    call DistributedGvcfyper {
      input:
        gvcf_list = gvcf_list,
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        gvcftyper_shards = GenerateShards.gvcftyper_shards[i],
        gvcftyper_memory = gvcftyper_memory,
        gvcftyper_threads = gvcftyper_threads,
        concurrent_downloads = concurrent_downloads,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_tbi = dbsnp_vcf_tbi,
        driver_xargs = driver_xargs,
        gvcftyper_xargs = gvcftyper_xargs,
        sentieon_docker = sentieon_docker,

        canonical_user_id = canonical_user_id,
        sentieon_license = sentieon_license,
    }
  }
  Array[File] flat_shard_vcfs = flatten(DistributedGvcfyper.shard_vcfs)
  Array[File] flat_shard_tbis = flatten(DistributedGvcfyper.shard_tbis)

  # MergeGvcftyper output:
  #  MergeGvcftyper.vcf
  #  MergeGvcftyper.tbi
  scatter(i in range(length(GenerateShards.merge_start))) {
    Int this_merge_start = GenerateShards.merge_start[i]
    String this_merge_shard = GenerateShards.merge_shards[i]
    call MergeGvcftyper {
      input:
        shard_vcfs = flat_shard_vcfs,
        shard_tbis = flat_shard_tbis,
        merge_start = this_merge_start,
        merge_size = GenerateShards.merge_size,
        merge_shard = this_merge_shard,
        merge_memory = merge_memory,
        merge_threads = merge_threads,
        sentieon_docker = sentieon_docker,

        canonical_user_id = canonical_user_id,
        sentieon_license = sentieon_license,
    }
  }
  Array[File] merged_vcfs = MergeGvcftyper.vcf
  Array[File] merged_tbis = MergeGvcftyper.tbi

  output {
    Array[Array[String]] gvcftyper_shards = GenerateShards.gvcftyper_shards
    Array[Int] merge_start = GenerateShards.merge_start
    Int merge_size = GenerateShards.merge_size
    Array[String] merge_shards = GenerateShards.merge_shards
    Array[File] out_shard_vcfs = flat_shard_vcfs
    Array[File] out_shard_tbis = flat_shard_tbis
    Array[File] out_merge_vcfs = merged_vcfs
    Array[File] out_merge_tbis = merged_tbis
  }
}

task GenerateShards {
  input {
    # Reference genome files
    File ref_fai

    # Settings for distribution
    String? region  # Can be used to limit the joint call to a specific region. Defaults to the whole genome
    Int compute_workers
    Int shard_size
    Int output_splits
    String sentieon_docker
  }

  command <<<
    set -exvuo pipefail

    region_arg=()
    region=~{default="" region}
    if [[ -n "$region" ]]; then
      region_arg=("--region" "$region")
    fi

    set +u
    /opt/sentieon/generate_shards.py \
      --fai "~{ref_fai}" \
      "${region_arg[@]}" \
      --n_workers "~{compute_workers}" \
      --shard_size "~{shard_size}" \
      --n_merge_splits "~{output_splits}" \
      --shards_json "gvcftyper_shards.json" \
      --merge_info "merge_info.json"
    set -u

    cat "merge_info.json" | jq '.size' > "merge_size.txt"
    cat "merge_info.json" | jq '.shards | [.[][0]]' > "merge_starts.json"
    cat "merge_info.json" | jq '.shards | [.[][1]]' > "merge_shards.json"
  >>>
  runtime {
    docker: sentieon_docker
    memory: "2 GiB"
    cpu: 1
  }
  output {
    Array[Array[String]] gvcftyper_shards = read_json("gvcftyper_shards.json")
    Int merge_size = read_int("merge_size.txt")
    Array[Int] merge_start = read_json("merge_starts.json")
    Array[String] merge_shards = read_json("merge_shards.json")
  }
}

task DistributedGvcfyper {
  input {
    # Input gVCFs
    #Array[File] gvcfs
    File gvcf_list  # Input gVCF files; one file per line

    # Reference genome files
    File ref_fasta
    File ref_fai

    Array[String] gvcftyper_shards

    # Settings for distribution
    String gvcftyper_memory = "64 GiB"
    String gvcftyper_threads = "16"
    Int concurrent_downloads = 5
    String sentieon_docker

    # GVCFtyper arguments
    File? dbsnp_vcf
    File? dbsnp_vcf_tbi
    String driver_xargs = "--traverse_param 10000/200"
    String gvcftyper_xargs = "--genotype_mode multinomial --max_alt_alleles 12"

    # Sentieon license configuration
    String canonical_user_id
    String sentieon_license = "aws-omics.sentieon.com:9011"
  }

  command <<<
    set -xv
    source /opt/sentieon/omics_credentials.sh "~{sentieon_license}" "~{canonical_user_id}"
    set -exvuo pipefail

    dbsnp_arg=()
    dbsnp_vcf=~{default="" dbsnp_vcf}
    if [[ -n "$dbsnp_vcf" ]]; then
      ln -s "$dbsnp_vcf" ./dbsnp.vcf.gz
      ln -s '~{default="" dbsnp_vcf_tbi}' ./dbsnp.vcf.gz.tbi
      dbsnp_arg=("--dbsnp" "./dbsnp.vcf.gz")
    fi

    ln -s "~{ref_fasta}" ./ref.fa
    ln -s "~{ref_fai}" ./ref.fa.fai

    set +u
    /opt/sentieon/sharded_joint_call.py \
      --ref "./ref.fa" \
      --gvcf_list "~{gvcf_list}" \
      --shards "~{sep='" "' gvcftyper_shards}" \
      --debug \
      "${dbsnp_arg[@]}" \
      --driver_xargs "~{driver_xargs}" \
      --gvcftyper_xargs "~{gvcftyper_xargs}" \
      --output_basename "shard_vcfs" \
      --concurrent_downloads "~{concurrent_downloads}"
    set -u
  >>>
  runtime {
    docker: sentieon_docker
    memory: gvcftyper_memory
    cpu: gvcftyper_threads
  }
  output {
    Array[File] shard_vcfs = glob("shard_vcfs*.vcf.gz")
    Array[File] shard_tbis = glob("shard_vcfs*.vcf.gz.tbi")
  }
}

task MergeGvcftyper {
  input {
    Array[File] shard_vcfs
    Array[File] shard_tbis

    # Merge info
    Int merge_start
    Int merge_size
    String merge_shard

    String merge_memory
    String merge_threads
    String sentieon_docker

    # Sentieon license configuration
    String canonical_user_id
    String sentieon_license = "aws-omics.sentieon.com:9011"
  }

  command <<<
    set -xv
    source /opt/sentieon/omics_credentials.sh "~{sentieon_license}" "~{canonical_user_id}"
    set -exvuo pipefail

    sharded_vcfs=("~{sep='" "' shard_vcfs}")
    sharded_tbis=("~{sep='" "' shard_tbis}")

    subset_vcfs=(${sharded_vcfs[@]:~{merge_start}:~{merge_size}})
    subset_tbis=(${sharded_tbis[@]:~{merge_start}:~{merge_size}})

    all_vcfs=()
    for i in $(seq 1 ${#subset_vcfs[@]}); do
      i=$((i - 1))
      ln -s $(realpath "${subset_vcfs[$i]}") ./called_${i}.vcf.gz
      ln -s $(realpath "${subset_tbis[$i]}") ./called_${i}.vcf.gz.tbi
      all_vcfs+=("./called_${i}.vcf.gz")
    done

    /opt/sentieon/merge_sharded_vcfs.py \
      --vcfs "${all_vcfs[@]}" \
      --region "~{merge_shard}" \
      --output merged_region.vcf.gz
  >>>
  runtime {
    docker: sentieon_docker
    memory: merge_memory
    cpu: merge_threads
  }
  output {
    File vcf = "merged_region.vcf.gz"
    File tbi = "merged_region.vcf.gz.tbi"
  }
}

