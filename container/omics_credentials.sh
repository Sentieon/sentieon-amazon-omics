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

INITIAL_LICENSE_ATTEMPTS=5
INITIAL_LICENSE_SLEEP=15

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
        LICENSE_SLEEP=$(awk -v seed=$RANDOM -v s_time=$INITIAL_LICENSE_SLEEP 'BEGIN{srand(seed); print rand() * s_time}')
        sleep "$LICENSE_SLEEP"
        if aws s3 cp "$LICENSE_URI" "$LICENSE_TMP_PATH"; then
            <"$LICENSE_TMP_PATH" base64 > "$AUTH_DATA_PATH"
        fi
    done
else
    # Check for a license token
    i=0
    while [ $i -lt $INITIAL_LICENSE_ATTEMPTS ]; do
        if aws s3 cp "$LICENSE_URI" "$LICENSE_TMP_PATH"; then
            break
        fi
        LICENSE_SLEEP=$(awk -v seed=$RANDOM -v s_time=$INITIAL_LICENSE_SLEEP 'BEGIN{srand(seed); print rand() * s_time}')
        sleep "$LICENSE_SLEEP"
        i=$((i+1))
    done
    if [ $i -ge $INITIAL_LICENSE_ATTEMPTS ]; then
        echo -e "ERROR: Unable to access a license token at $LICENSE_URI.\nPlease email support@sentieon.com to request access to the Sentieon license for Amazon Omics"
        exit 1
    fi
    <"$LICENSE_TMP_PATH" base64 > "$AUTH_DATA_PATH"

    # Configure the environment for the Sentieon tools
    export SENTIEON_LICENSE
    export SENTIEON_AUTH_MECH=aws_omics_service
    export SENTIEON_AUTH_DATA=~/.sentieon/sentieon_license_encode
    SENTIEON_JOB_TAG=$(<"$LICENSE_TMP_PATH" jq '.job_tag')
    export SENTIEON_JOB_TAG

    # Check that OMICS_NETWORK_PROXY is set
    OMICS_NETWORK_PROXY="${OMICS_NETWORK_PROXY-}"
    if [ -z "$OMICS_NETWORK_PROXY" ]; then
        echo "Error: the variable OMICS_NETWORK_PROXY is required but is unset. Please reach out to the AWS support team so that your account can be allow-listed with the proxy server."
        exit 1
    fi

    # Configure the software to use the proxy server
    export http_proxy="${OMICS_NETWORK_PROXY}"

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
