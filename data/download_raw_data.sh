#!/bin/bash

while read -r sample url; do
  echo "Downloading $sample from ENA..."
  curl $url -o $sample"_raw.fq.gz"
done < samples.txt
echo "Done"

