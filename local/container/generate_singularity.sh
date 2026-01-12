#!/usr/bin/env bash

set -e  # Exit on error
set -u  # Exit on undefined variable

# Use provided name or default to transcript_analysis.sif
CONTAINER_NAME="${1:-transcript_analysis.sif}"

if [ -f "${CONTAINER_NAME}" ]; then
    echo "Container already exists: ${CONTAINER_NAME}. Skipping build."
    exit 0
fi

echo "Building Singularity container: ${CONTAINER_NAME}..."
singularity build --remote "${CONTAINER_NAME}" transcript.def

echo "Container built successfully: ${CONTAINER_NAME}"