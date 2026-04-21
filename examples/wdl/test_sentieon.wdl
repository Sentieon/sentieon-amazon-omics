version 1.1

workflow sentieon_germline {
  input {
    # Sentieon license configuration
    String sentieon_docker
    String sentieon_license
  }
  # Perform a license check
  call SentieonLicense {
    input:
      sentieon_docker = sentieon_docker,
      sentieon_license = sentieon_license,
  }
  output {
    File license_ok = SentieonLicense.license_ok
  }
}

task SentieonLicense {
  input {
    String sentieon_docker
    String sentieon_license
  }
  command <<<
    set -exvuo pipefail
    export SENTIEON_LICENSE="~{sentieon_license}"

    # Test Sentieon commands
    sentieon licclnt ping && echo "Ping is OK"
    sentieon licclnt query Haplotyper

    # Can run custom Sentieon commands here:
    #  sentieon bwa mem ...
    #  sentieon driver ... --algo Haplotyper ...

    echo "License OK" >license_ok.txt
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
