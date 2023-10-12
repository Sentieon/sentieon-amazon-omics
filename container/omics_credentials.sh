#!/usr/bin/env bash

# Do nothing if we are not in Omics
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
if [ -z "${OMICS_NETWORK_PROXY-}" ]; then
    . $SCRIPT_DIR/$(basename ${BASH_SOURCE[0]}).orig "$@"
    exit $!
fi

#set -xvu
set +e

mkdir -p ~/.sentieon

if [ -z "${SENTIEON_LICENSE-}" ]; then
    echo "Please set the environment variable SENTIEON_LICENSE."
    exit 1
fi
if [ -z "${CANONICAL_USER_ID-}" ]; then
    echo "Please set the environment variable CANONICAL_USER_ID."
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

if [ -n "${IS_DAEMON-}" ]; then
    # The license file URI
    LICENSE_URI="s3://$SENTIEON_BUCKET_BASENAME-$AWS_DEFAULT_REGION/${CANONICAL_USER_ID}.lic"

    # Refresh the license
    while true; do
        sleep "$SLEEP_TIME"
        aws s3 cp "$LICENSE_URI" "$LICENSE_TMP_PATH"
        <"$LICENSE_TMP_PATH" base64 > "$AUTH_DATA_PATH"
    done
    exit 0
fi

# License setup if the encoded token is not found
if [ ! -f "$AUTH_DATA_PATH" ]; then
    # Find the AWS Region
    if [ -z "${AWS_DEFAULT_REGION-}" ]; then
        set -e
        curl "${ECS_CONTAINER_METADATA_URI_V4}/task" > ~/.sentieon/tmp_task.json
        AWS_DEFAULT_REGION=$(<~/.sentieon/tmp_task.json jq -rM '.AvailabilityZone' | sed 's/.$//')
        export AWS_DEFAULT_REGION
        rm ~/.sentieon/tmp_task.json
        set +e
    fi

    # The license file URI
    LICENSE_URI="s3://$SENTIEON_BUCKET_BASENAME-$AWS_DEFAULT_REGION/${CANONICAL_USER_ID}.lic"

    # Check for a license token
    aws s3 cp "$LICENSE_URI" "$LICENSE_TMP_PATH"
    token_ok="$?"
    if [ "$token_ok" -ne 0 ]; then
        echo -e "ERROR: Unable to access a license token at $LICENSE_URI.\nPlease email support@sentieon.com to request access to the Sentieon license for Amazon Omics"
        exit 1
    fi
    <"$LICENSE_TMP_PATH" base64 > "$AUTH_DATA_PATH"

    # Kickoff the token refresh script
    IS_DAEMON=true bash -c "cd /; setsid bash \"${BASH_SOURCE[0]}\" </dev/null &>/dev/null &"
fi

# Configure the environment
export SENTIEON_AUTH_MECH=aws_omics_service
export SENTIEON_AUTH_DATA="$AUTH_DATA_PATH"
SENTIEON_JOB_TAG=$(<"$LICENSE_TMP_PATH" jq '.job_tag')
export SENTIEON_JOB_TAG
export http_proxy=${OMICS_NETWORK_PROXY}

# Run the original command
. $SCRIPT_DIR/$(basename ${BASH_SOURCE[0]}).orig "$@"
