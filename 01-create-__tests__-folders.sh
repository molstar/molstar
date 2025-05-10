#!/bin/bash

# Add __tests__ directories to appropriate locations
find src -type f -name '*.ts' ! -name '*.d.ts' -exec dirname {} \; | sort -u | while read dir; do 
  mkdir -p "$dir/__tests__"; 
done