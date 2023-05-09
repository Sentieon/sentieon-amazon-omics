version 1.1

workflow sentieon_germline {
  input {
    # Sentieon license configuration
    String canonical_user_id
    String sentieon_docker
    String sentieon_license = "aws-omics.sentieon.com:9011"
  }
  # Perform a license check
  call SentieonLicense {
    input:
      canonical_user_id = canonical_user_id,
      sentieon_docker = sentieon_docker,
      sentieon_license = sentieon_license,
  }
  output {
    File license_ok = SentieonLicense.license_ok
  }
}

task SentieonLicense {
  input {
    String canonical_user_id
    String sentieon_docker
    String sentieon_license = "aws-omics.sentieon.com:9011"
  }
  command <<<
    set -xv

    set +e
    source /opt/sentieon/omics_credentials.sh "~{sentieon_license}" "~{canonical_user_id}"
    set -e

    # Test Sentieon commands
    sentieon licclnt ping && echo "Ping is OK"
    sentieon licclnt query Haplotyper

    # Can run custom Sentieon commands here:
    #  sentieon bwa mem ...
    #  sentieon driver ... --algo Haplotyper ...

    echo "License OK" >license_ok.txt
    unset http_proxy https_proxy
    sleep 10
  >>>
  runtime {
    container: sentieon_docker
    memory: "1 GB"
    cpu: 1
  }
  output {
    File license_ok = "license_ok.txt"
  }
}
