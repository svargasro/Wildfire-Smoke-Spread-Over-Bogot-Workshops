#!/bin/bash
for dir in */; do
  for file in "$dir"*_*; do
    mv "$file" "$(echo "$file" | sed 's/_[0-9]\{2\}-[0-9]\{2\}-[0-9]\{4\}//')"
  done
done
