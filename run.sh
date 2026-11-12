#!/bin/sh
# Runs the engine offline. Order: portable Python in ../_runtime/python, then python3, then python.
cd "$(dirname "$0")" || exit 1
for p in ../_runtime/python/bin/python3 python3 python; do
  if command -v "$p" >/dev/null 2>&1; then exec "$p" hvac_engine.py "$@"; fi
done
echo "Python 3.6+ not found. Install Python, or place a portable Python in ../_runtime/python/"; exit 1
