/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import {
  GeneratedBedFile,
  generatedGenesBed,
  MisreadBedFile,
  MisreadGffFile,
  misreadBed,
  misreadGff,
} from '../genomes/files.js';
import { PackageFolder } from '../node-assets.js';

// The genome formats' generated files are their generators' output, byte for byte (the provenance rule's "a file a
// generator writes and a test compares").

describe('the genome formats’ generated files', () => {
  it('generated-genes.bed, misread.bed and misread.gff3 are what their generators write', () => {
    expect(readFileSync(`${PackageFolder}${GeneratedBedFile}`, 'utf8')).toBe(generatedGenesBed());
    expect(readFileSync(`${PackageFolder}${MisreadBedFile}`, 'utf8')).toBe(misreadBed());
    expect(readFileSync(`${PackageFolder}${MisreadGffFile}`, 'utf8')).toBe(misreadGff());
  });
});
