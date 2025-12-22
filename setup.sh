#!/bin/bash

echo "Setting up the project environment..."

SOURCE_DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
cd $SOURCE_DIR

python -m venv venv
source venv/bin/activate

pip install --upgrade pip
pip install -r requirements.txt

echo "Building project..."

./build.sh
