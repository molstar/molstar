/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync, writeFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import {
  GeneratedBedFile,
  GeneratedVcfFile,
  generatedGenesBed,
  GeneratedGenotypesFile,
  generatedGenotypesVcf,
  generatedVariantsVcf,
  MisreadBedFile,
  MisreadGffFile,
  MisreadVcfFile,
  misreadBed,
  misreadGff,
  misreadVcf,
  MisreadGenbankFile,
  misreadGenbank,
} from '../genomes/files.js';
import { PackageFolder } from '../node-assets.js';

// The genome formats' generated files are their generators' output, byte for byte (the provenance rule's "a file a
// generator writes and a test compares").

describe('the genome formats’ generated files', () => {
  it('generated-genes.bed, generated-variants.vcf, generated-genotypes.vcf and the misread files are what their generators write', () => {
    expect(readFileSync(`${PackageFolder}${GeneratedBedFile}`, 'utf8')).toBe(generatedGenesBed());
    expect(readFileSync(`${PackageFolder}${MisreadBedFile}`, 'utf8')).toBe(misreadBed());
    expect(readFileSync(`${PackageFolder}${MisreadGffFile}`, 'utf8')).toBe(misreadGff());
    expect(readFileSync(`${PackageFolder}${GeneratedVcfFile}`, 'utf8')).toBe(generatedVariantsVcf());
    if (process.env.SEQ_RECORDS === 'update')
      writeFileSync(`${PackageFolder}${GeneratedGenotypesFile}`, generatedGenotypesVcf());
    expect(readFileSync(`${PackageFolder}${GeneratedGenotypesFile}`, 'utf8')).toBe(generatedGenotypesVcf());
    expect(readFileSync(`${PackageFolder}${MisreadVcfFile}`, 'utf8')).toBe(misreadVcf());
    if (process.env.SEQ_RECORDS === 'update') writeFileSync(`${PackageFolder}${MisreadGenbankFile}`, misreadGenbank());
    expect(readFileSync(`${PackageFolder}${MisreadGenbankFile}`, 'utf8')).toBe(misreadGenbank());
  });
});
