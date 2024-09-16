# Sentieon-Amazon-Omics
Sentieon pipelines for AWS HealthOmics

## Introduction

Sentieon supports bioinformatic workflows running on [AWS HealthOmics](https://aws.amazon.com/healthomics/). The files in this repository can be used to run Sentieon pipelines as private workflows on AWS HealthOmics or you can use this repository as a starting point for developing customized pipelines that utilize the Sentieon software.

The Sentieon software is a commercial software package and a license is required to run the software. To support workflows running on AWS HealthOmics, Sentieon operates a dedicated license server for AWS HealthOmics workflows.

To use the license server for AWS HealthOmics, you will need to provide Sentieon (support@sentieon.com) your AWS Canonical User ID. You can find your canonical ID by following the instructions at the following link, https://docs.aws.amazon.com/accounts/latest/reference/manage-acct-identifiers.html#FindCanonicalId.

## Running Sentieon pipelines as private workflows

### Requirements
* Docker cli or another container implementation (Podman, etc.)
* AWS CLI v2

### Step 0: Add your account to the allowlist for the AWS HealthOmics proxy server

Generate an AWS support case to get access to the Sentieon license server proxy. To create a support case, navigate to https://support.console.aws.amazon.com/. Provide your AWS account and Region in the support case. Your account will be added to the allowlist for the license server proxy.

### Step 1: build the Sentieon container image

The following files are in the [`container`](/container) directory:
* `sentieon_omics.dockerfile`: A dockerfile that can be used to create a Sentieon container image for AWS HealthOmics
* `omics_credentials.sh`: a shell script to perform license authentication on AWS HealthOmics

To build the container image for the latest version of Sentieon, run:
```bash
cd ./container
docker build --platform linux/amd64 --build-arg SENTIEON_VERSION=202112.07 -t sentieon:omics-1 -f sentieon_omics.dockerfile .
```

### Step 2: push the container image to an Amazon ECR private repository

Create a private repository in AWS ECR

```bash
aws ecr create-repository --repository-name sentieon
```

Login to the registry

```bash
aws ecr get-login-password --region <region-name> | docker login --username AWS --password-stdin <account-id>.dkr.ecr.<region-name>.amazonaws.com
```

Tag the custom Sentieon container and push the container image to the repository

```bash
docker tag sentieon:omics-1 <account-id>.dkr.ecr.<region-name>.amazonaws.com/sentieon:omics-1
docker push <account-id>.dkr.ecr.<region-name>.amazonaws.com/sentieon:omics-1
```

Grant the HealthOmics service permission to interact with the repository using the policy in the `assets` directory

```bash
aws ecr set-repository-policy --repository-name sentieon --policy-text file://assets/omics-ecr-repository-policy.json
```

### Step 3: grant the HealthOmics service role read access to the Sentieon license bucket in AWS s3

As part of the license validation, the `omics_credentials.sh` script will obtain a license token from AWS s3 for your workflow. Adding the following policy to your AWS HealthOmics service role will grant the workflow read access to files in the license bucket for your region:
```json
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "s3:GetObjectAcl",
                "s3:GetObject"
            ],
            "Resource": [
                "arn:aws:s3:::sentieon-omics-license-<region>/*"
            ]
        }
    ]
}
```

### Step 4: create an example workflow on AWS HealthOmics

We are now ready to create Sentieon workflows on AWS HealthOmics. Running the following command at the start of the workflow will configure the environment for the Sentieon software:

```bash
source /opt/sentieon/omics_credentials.sh <SENTIEON_LICENSE> <CANONICAL_USER_ID>
```
Where `<SENTIEON_LICENSE>` is the FQDN and port of the Sentieon license server and  `<CANONICAL_USER_ID>` is the AWS canonical user ID of the account running the workflow.

Example workflows can be found in the [`examples`](/examples) directory and complete workflow implementations can be found in the [`workflows`](/workflows) directory.

#### WDL

```bash
(cd examples/wdl && zip test_sentieon.wdl.zip test_sentieon.wdl)

aws omics create-workflow \
    --name test-sentieon-wdl \
    --engine WDL \
    --definition-zip fileb://examples/wdl/test_sentieon.wdl.zip \
    --parameter-template file://examples/parameter.template.json
```

#### Nextflow
```bash
(cd examples/nextflow && zip -r ${OLDPWD}/test_sentieon.nextflow.zip .)

aws omics create-workflow \
    --name test-sentieon-nextflow \
    --engine NEXTFLOW \
    --main test_sentieon.nf \
    --definition-zip fileb://test_sentieon.nextflow.zip \
    --parameter-template file://examples/parameter.template.json
```

The `create-workflow` command will output some information including the workflow-id.

### Step 5: run the example workflow

To run the example workflow, modify the `examples/test.parameters.json` file replacing `<canonical-id>`, `<account-id>`, and `<region-name>` to match your environment. Then run the following, using the `workflow-id` from the `create-workflow` command and the `role-name` for your AWS HealthOmics service role:

```bash
aws omics start-run \
    --role-arn "arn:aws:iam::<account-id>:role/<role-name>" \
    --workflow-id <workflow_id> \
    --name "test $(date +%Y%m%d-%H%M%S)" \
    --output-uri <s3-uri> \
    --parameters file://examples/test.parameters.json
```

After ~20min, verify that the test workflow completes successfully:

```bash
aws omics get-run --id <run-id>
```

You should see a response like:
```json
{
    "arn": "arn:aws:omics:<region>:<account-id>:run/<run-id>",
    "creationTime": "2023-04-24T17:14:38.880864+00:00",
    "digest": "sha256:<sha256>",
    "id": "<run-id>",
    "name": "test 20230424-101437",
    "outputUri": "<s3-output-uri>",
    "parameters": {
        "canonical_user_id": "<canonical_id>",
        "sentieon_docker": "<account-id>.dkr.ecr.<region>.amazonaws.com/sentieon:omics"
    },
    "resourceDigests": {
        "<account-id>.dkr.ecr.<region>.amazonaws.com/sentieon:omics": "sha256:<sha256>"
    },
    "roleArn": "arn:aws:iam::<account-id>:role/<role-name>",
    "startTime": "2023-04-24T17:25:06.021000+00:00",
    "startedBy": "arn:aws:iam::<account-id>:<user>",
    "status": "COMPLETED",
    "stopTime": "2023-04-24T17:39:32.138095+00:00",
    "tags": {},
    "workflowId": "<workflow-id>",
    "workflowType": "PRIVATE"
}
```

The example workflow has only one task called `SentieonLicense`. Locate this task:

```bash
aws omics list-run-tasks --id <run-id>
```

You should see a response like:
```json
{
    "items": [
        {
            "cpus": 1,
            "creationTime": "2023-04-24T17:25:40.323030+00:00",
            "memory": 4,
            "name": "SentieonLicence",
            "startTime": "2023-04-24T17:29:42.881000+00:00",
            "status": "COMPLETED",
            "stopTime": "2023-04-24T17:30:24.442000+00:00",
            "taskId": "<task-id>"
        }
    ]
}
```

Get the log-stream for the task:
```bash
aws logs get-log-events --log-group-name /aws/omics/WorkflowLog --log-stream-name run/<run-id>/task/<task-id> --output text
```
*Note that `get-log-events` is paginated, and may not return the full log stream for workflows with verbose logs*


If license verification is successful, you should see event lines like:
```text
EVENTS  1682357406252   sentieon licclnt ping && echo "Ping is OK"      1682357399707
EVENTS  1682357406252   + sentieon licclnt ping 1682357399707
EVENTS  1682357406252   + echo 'Ping is OK'     1682357400013
EVENTS  1682357406252   sentieon licclnt query Haplotyper       1682357400013
EVENTS  1682357406252   + sentieon licclnt query Haplotyper     1682357400015
EVENTS  1682357406252   Ping is OK      1682357400015
EVENTS  1682357406252   499968  1682357400539
```

## Next steps

Congratulations! You've successfully run a test workflow with the Sentieon software on AWS HealthOmics. Feel free to update/extend the example workflow to implement your own custom pipelines with the Sentieon software.

Alternatively, you can find full pipeline implementations in the [`workflows`](/workflows) directory that you can modify or implement as private workflows.
