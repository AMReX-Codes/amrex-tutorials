#!/usr/bin/env bash

GITHUB_EVENT_API_URL="https://api.github.com/repos/AMReX-Codes/amrex-tutorials/events"
GITHUB_PR_API="https://api.github.com/repos/AMReX-Codes/amrex-tutorials/pulls"
MIRROR_TRIGGER_URL="https://software.nersc.gov/api/v4/projects/307/trigger/pipeline"


set -eu -o pipefail

sudo apt-get update

sudo apt-get install -y --no-install-recommends\
  jq   \
  curl


# indetify the PR associated with second most recent comment

PR_ID=$(curl -H "Accept: application/vnd.github.v3+json" ${GITHUB_EVENT_API_URL} | jq .[1].payload.issue.number)


# identify forked branch from PR

CONTENT=$(curl -H "Accept: application/vnd.github.v3+json" ${GITHUB_PR_API}/${PR_ID})


PR_BRANCH=$(jq -r '.head.ref' <<< "${CONTENT}")

echo PR ID Number: $PR_ID
echo PR Branchname: $PR_BRANCH

#Trigger GitLab CI, pass variables
# succeeds at command line

curl -s -X POST\
     -F "token=$1" \
     -F "ref=main" \
     -F "variables[PR_ID]=${PR_ID}" \
     -F "variables[PR_BRANCH]=${PR_BRANCH}" \
     ${MIRROR_TRIGGER_URL}






