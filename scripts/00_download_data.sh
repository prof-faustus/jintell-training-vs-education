#!/usr/bin/env bash
# Downloads the Open LLM Leaderboard contents table at the exact revision analysed and verifies it,
# then unpacks the Epoch AI benchmark data (CC BY 4.0) shipped in data/epoch/.
set -euo pipefail
REV=9c09a7cae43334062a82cb164f2ef255013dafa2
SHA=f13a7061e63133bc8745d89f8c7a56148584940ede22d3b27e21cad1e078dca5
mkdir -p data/hf
if [ ! -f data/hf/contents.parquet ]; then
  curl -L -sS -o data/hf/contents.parquet \
    "https://huggingface.co/datasets/open-llm-leaderboard/contents/resolve/${REV}/data/train-00000-of-00001.parquet"
fi
echo "${SHA}  data/hf/contents.parquet" | sha256sum -c -
echo "7fd923b837532422a0a35aabec72c753b6f554fd425221fa2bb3aeb9948ffe34  data/epoch/benchmark_data.zip" | sha256sum -c -
rm -rf data/epoch/extracted && mkdir -p data/epoch/extracted && unzip -q data/epoch/benchmark_data.zip -d data/epoch/extracted
