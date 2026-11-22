/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createHash } from 'node:crypto';
import { existsSync, readdirSync, readFileSync, statSync } from 'node:fs';
import { join } from 'node:path';
import { describe, expect, it } from 'vitest';
import { AllowedLicences } from '../files/record-trim.js';
import { GeneFile } from '../p53/data.js';
import { ConservationFile } from '../p53/conservation.js';
import { RegionFile } from '../p53/region.js';
import {
  GeneratedBedFile,
  GeneratedVcfFile,
  MisreadBedFile,
  MisreadGffFile,
  MisreadVcfFile,
} from '../genomes/files.js';
import { PackageFolder } from '../node-assets.js';

// Provenance is data (examples.md, "One rule for every fixture"): each dataset's `sources.json` lists every file of its
// folder (a file a generator writes and a test compares instead is named here), each row with its entry, URL, date,
// licence (one allowed), attribution, what was trimmed and derived, and the size and SHA-256 of the file as committed;
// a file under 1 MB, a dataset under 2 MB.

interface Row {
  readonly file: string;
  readonly entry: string;
  readonly url: string;
  readonly recorded: string;
  readonly licence: readonly string[];
  readonly attribution: string;
  readonly trimmed: string;
  readonly derived: string;
  readonly bytesDownloaded: number;
  readonly bytesFixture: number;
  readonly sha256Download: string;
  readonly sha256Fixture: string;
}

const Fixtures = join(PackageFolder, 'fixtures');
/** The datasets recorded by the script (the `files` folder is the Files feature's own tables). */
const Datasets = ['p53', 'globins', 'complex', 'antibody', 'services', 'genomes'];
/** Files a generator writes, compared with it by a test (not recorded): by their path under `fixtures/`. */
const Generated = new Set(
  [
    GeneFile,
    RegionFile,
    ConservationFile,
    GeneratedBedFile,
    MisreadBedFile,
    MisreadGffFile,
    GeneratedVcfFile,
    MisreadVcfFile,
  ].map((f) => f.slice('fixtures/'.length)),
);

const MaxFile = 1024 * 1024;
const MaxDataset = 2 * 1024 * 1024;

describe.each(Datasets)('fixtures/%s/sources.json', (dataset) => {
  const dir = join(Fixtures, dataset);
  const rows = JSON.parse(readFileSync(join(dir, 'sources.json'), 'utf8')) as Row[];

  it('lists every file of the folder once, and only those', () => {
    const files = readdirSync(dir)
      .filter((f) => f !== 'sources.json')
      .map((f) => `${dataset}/${f}`)
      .filter((f) => !Generated.has(f))
      .sort();
    // a row whose file is `../src/…` is held in a source module (not a file of the folder), its checksum the module's
    expect(
      rows
        .map((r) => r.file)
        .filter((f) => !f.startsWith('../'))
        .sort(),
    ).toEqual(files);
  });

  it('each row: its provenance said, a licence allowed, its size and checksum the file’s', () => {
    for (const r of rows) {
      const text = readFileSync(join(Fixtures, r.file));
      expect(r.entry, r.file).not.toBe('');
      expect(r.url, r.file).toMatch(/^https:\/\//);
      expect(r.recorded, r.file).toMatch(/^\d{4}-\d{2}-\d{2}$/);
      expect(r.attribution, r.file).not.toBe('');
      expect(r.trimmed, r.file).not.toBe('');
      expect(r.derived, r.file).not.toBe('');
      expect(r.licence.length, r.file).toBeGreaterThan(0);
      for (const l of r.licence) expect(AllowedLicences, r.file).toContain(l);
      expect(r.bytesFixture, r.file).toBe(text.byteLength);
      expect(r.sha256Fixture, r.file).toBe(createHash('sha256').update(text).digest('hex'));
      expect(r.sha256Download, r.file).toMatch(/^[0-9a-f]{64}$/);
      expect(r.bytesDownloaded, r.file).toBeGreaterThan(0);
    }
  });

  it('a file under 1 MB, the dataset under 2 MB', () => {
    let total = 0;
    for (const f of readdirSync(dir)) {
      const size = statSync(join(dir, f)).size;
      expect(size, f).toBeLessThan(MaxFile);
      total += size;
    }
    expect(total).toBeLessThan(MaxDataset);
  });
});

it('every generated file named is there', () => {
  for (const f of Generated) expect(existsSync(join(Fixtures, f)), f).toBe(true);
});
