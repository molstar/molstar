#!/bin/bash

# Find all __tests__ directories and create 3 subdirectories in each
find . -type d -name "__tests__" | while read tests_dir; do
  mkdir -p "$tests_dir/jest"
  mkdir -p "$tests_dir/vitest"
  echo "Created jest/ and vitest/ inside $tests_dir"
done

