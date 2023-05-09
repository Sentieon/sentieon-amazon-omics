nextflow.enable.dsl=2

process SentieonLicence {
    label 'sentieon'

    output:
        path "license_ok.txt", emit: license_ok
    
    script:
        """
        set -xv

        set +e
        source /opt/sentieon/omics_credentials.sh "${params.sentieon_license}" "${params.canonical_user_id}"
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
        """
}

// implicit workflow
workflow {
    SentieonLicence()
}
