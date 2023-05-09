#!/usr/bin/env bash

set -xvu
set +e

SENTIEON_LICENSE="${1-}"
CANONICAL_USER_ID="${2-}"
IS_DAEMON="${3-}"

mkdir -p ~/.sentieon

if [ -z "$SENTIEON_LICENSE" ] || [ -z "$CANONICAL_USER_ID" ]; then
    echo "Usage: source ${BASH_SOURCE[0]} <SENTIEON_LICENSE> <AWS_canonical_user_ID>"
    exit 1
fi
if [ "${#CANONICAL_USER_ID}" -ne 64 ]; then
    echo -e "Error: The CANONICAL_USER_ID does not have the expected length.\nPlease confirm your canonical user id following the instructions at, https://docs.aws.amazon.com/AmazonS3/latest/userguide/finding-canonical-user-id.html.\nA correct canonical user ID is required for license validation."
    exit 1
fi


LICENSE_TMP_PATH=~/.sentieon/sentieon_license
AUTH_DATA_PATH=~/.sentieon/sentieon_license_encode
SENTIEON_BUCKET_BASENAME="sentieon-omics-license"
SLEEP_TIME=900

# Find the AWS Region
if [ -z "${AWS_DEFAULT_REGION-}" ]; then
    set -e
    curl "${ECS_CONTAINER_METADATA_URI_V4}/task" > tmp_task.json
    AWS_DEFAULT_REGION=$(<tmp_task.json jq -rM '.AvailabilityZone' | sed 's/.$//')
    export AWS_DEFAULT_REGION
    rm tmp_task.json
    set +e
fi

# The license file URI
LICENSE_URI="s3://$SENTIEON_BUCKET_BASENAME-$AWS_DEFAULT_REGION/${CANONICAL_USER_ID}.lic"

if [ -n "$IS_DAEMON" ]; then
    # Refresh the license
    while true; do
        sleep "$SLEEP_TIME"
        aws s3 cp "$LICENSE_URI" "$LICENSE_TMP_PATH"
        <"$LICENSE_TMP_PATH" base64 > "$AUTH_DATA_PATH"
    done
else
    # Check for a license token
    aws s3 cp "$LICENSE_URI" "$LICENSE_TMP_PATH"
    token_ok="$?"
    if [ "$token_ok" -ne 0 ]; then
        echo -e "ERROR: Unable to access a license token at $LICENSE_URI.\nPlease email support@sentieon.com to request access to the Sentieon license for Amazon Omics"
        exit 1
    fi
    <"$LICENSE_TMP_PATH" base64 > "$AUTH_DATA_PATH"

    # Configure the environment for the Sentieon tools
    export SENTIEON_LICENSE
    export SENTIEON_AUTH_MECH=aws_omics_service
    export SENTIEON_AUTH_DATA=~/.sentieon/sentieon_license_encode

    # Set `http_proxy` to the proxy server's IP
    ## Parse the proxy
    proxy_host_port=${OMICS_NETWORK_PROXY}
    proxy_port="${proxy_host_port##*:}"
    proxy_prefix_host="${proxy_host_port%:*}"
    proxy_host="${proxy_prefix_host##http://}"

    # Get the proxy ip address(es)
    getent ahosts "$proxy_host" | cut -f 1 -d ' ' | uniq > test_dns.txt
    sleep 1

    # Test the proxy ip with curl
    proxy_ip=""
    while IFS= read -r line; do
      echo "Testing ip address: $line"
      https_proxy=http://"$line":"$proxy_port" curl -v -m 5 -k https://aws-omics.sentieon.com:9011
      retcode="$?"  # curl sets a returncode of 28 for a timeout while 35 is expected for a handshake failure
      if [ "$retcode" -ne 28 ]; then
        proxy_ip="$line"
      fi
    done < test_dns.txt
    rm test_dns.txt
    if [ -z "$proxy_ip" ]; then
      echo "Error: unable to determine proxy IP"
      exit 1
    fi

    export OMICS_NETWORK_PROXY=http://"$proxy_ip":"$proxy_port"
    export http_proxy=${OMICS_NETWORK_PROXY}

    # Test the connection to the license server
    if sentieon licclnt ping; then
        echo "Successfully connected to the license server"
    else
        echo "Error: Unable to connect to the license server"
        exit 1
    fi

    # Kickoff the token refresh script
    bash -c "cd /; setsid bash \"${BASH_SOURCE[0]}\" \"$1\" \"$2\" 1 </dev/null &>/dev/null &"
fi
