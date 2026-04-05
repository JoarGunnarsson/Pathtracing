#!/bin/bash
set -e

echo "Setting up the python environment..."

PROJECT_DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
cd $PROJECT_DIR

python -m venv venv
source venv/bin/activate

pip install --upgrade pip
pip install -r requirements.txt
