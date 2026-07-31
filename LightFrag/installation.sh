#!/bin/bash

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Creating virtual environment..."
python3.11 -m venv "$SCRIPT_DIR/venv"

# Need pandas for the .csv files we make and to nicely send results to the output
echo "Installing dependencies..."
"$SCRIPT_DIR/venv/bin/python3.11" -m pip install --upgrade pip
"$SCRIPT_DIR/venv/bin/python3.11" -m pip install -r "$SCRIPT_DIR/requirements.txt"

echo "Installation complete"