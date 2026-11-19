/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createHash } from 'node:crypto';
import { readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { decodeSeq } from '@molstar/seqx-core/data';
import { readStockholm } from '@molstar/seqx-core/data/alignment-readers';
import { describe, expect, it } from 'vitest';
import { PackageFolder } from '../node-assets.js';
import { DerivedData } from '../trna/derived-data.js';
import { derive, derivedModule } from '../trna/derive.js';
import {
  covariationGrid,
  MoleculeFile,
  PairsFile,
  Pairs,
  SeedFile,
  scaleGrid,
  SystemName,
  TrnaLength,
  TrnaLetters,
} from '../trna/data.js';
import { cifRows, kindOfSaenger } from '../trna/record-trim.js';

// The tRNA fixtures (examples.md, tRNA): what the recorder wrote is what the dataset's code holds, its checksums are
// the sources file's, the Rfam seed is as large as it was recorded, and what is computed from them is what is committed
// (`SEQ_RECORDS=update` writes `derived-data.ts` again).

const text = (f: string) => readFileSync(`${PackageFolder}${f}`, 'utf8');
const update = process.env.SEQ_RECORDS === 'update';

/** A column of a ranges table by name. */
function field(table: any, name: string): any[] {
  const f = table.fields.find((x: { name: string }) => x.name === name).column;
  return f.dict ? Array.from(f.codes as ArrayLike<number>, (c: number) => f.dict[c]) : Array.from(f);
}

function ranges(file: string) {
  const { data, issues } = decodeSeq(text(file));
  expect(issues).toEqual([]);
  return data!.pairs[0].answer as any;
}

describe('the tRNA fixtures', () => {
  it('1EHZ: 76 letters, its 14 modified residues as ranges with their codes and parent letters', () => {
    const a = ranges(MoleculeFile);
    const letters: string = a.tables.letters.letters;
    expect(letters).toBe(TrnaLetters);
    expect(letters).toHaveLength(TrnaLength);
    const r = a.tables.ranges;
    expect(Array.from(r.start as Int32Array)).toEqual([9, 15, 16, 25, 31, 33, 36, 38, 39, 45, 48, 53, 54, 57]);
    expect(field(r, 'label').join(' ')).toBe('2MG H2U H2U M2G OMC OMG YYG PSU 5MC 7MG 5MC 5MU PSU 1MA');
    // a modified residue is written as its parent letter; the anticodon is Gm34-A35-A36
    expect(field(r, 'category').every((p: string, k: number) => p === letters[r.start[k]])).toBe(true);
    expect(letters.slice(33, 36)).toBe('GAA');
    expect(decodeSeq(text(MoleculeFile)).data!.coverage.systems[0].name).toBe(SystemName);
  });

  it('the pairs the code holds are the file’s: 27, six tertiary, the kinds from the entry’s own classes', () => {
    const r = ranges(PairsFile).tables.ranges;
    const rows = Array.from(r.start as Int32Array, (s, k) => [
      s,
      field(r, 'partner')[k],
      field(r, 'type')[k],
      field(r, 'tertiary')[k],
    ]);
    expect(rows).toEqual(
      Pairs.map(([i, j, k, t]) => [i, j, { C: 'canonical', W: 'wobble', N: 'non-canonical' }[k], t ? 'yes' : 'no']),
    );
    expect(Pairs.filter((p) => p[3]).map((p) => `${p[0] + 1}-${p[1] + 1}`)).toEqual([
      '8-14',
      '15-48',
      '18-55',
      '19-56',
      '26-44',
      '54-58',
    ]);
    expect(Pairs.filter((p) => p[2] === 'W').map((p) => `${p[0] + 1}-${p[1] + 1}`)).toEqual(['4-69']);
  });

  it('the seed is as published: 954 rows of 118 columns with SS_cons’ nested pairs', () => {
    const read = readStockholm(text(SeedFile));
    expect(read.rows).toHaveLength(954);
    expect(read.columns).toBe(118);
    expect(read.pairs!.length).toBeGreaterThanOrEqual(20);
    expect(read.issues).toEqual([]);
  });

  it('sources.json: every file’s bytes and checksum are the file’s', () => {
    const rows = JSON.parse(text('fixtures/trna/sources.json')) as {
      file: string;
      bytesFixture: number;
      sha256Fixture: string;
      licence: string[];
    }[];
    expect(rows.map((r) => r.file)).toEqual(['trna/1EHZ.seq', 'trna/1EHZ-pairs.seq', 'trna/RF00005.sto']);
    for (const r of rows) {
      const t = text(`fixtures/${r.file}`);
      expect(new TextEncoder().encode(t).byteLength).toBe(r.bytesFixture);
      expect(createHash('sha256').update(t).digest('hex')).toBe(r.sha256Fixture);
      expect(r.licence).toEqual(['CC0-1.0']);
    }
  });
});

describe('what is computed from them', () => {
  const letters = ranges(MoleculeFile).tables.letters.letters as string;
  const derived = derive(text(SeedFile), letters);

  it('is what is committed', () => {
    if (update)
      writeFileSync(fileURLToPath(new URL('../trna/derived-data.ts', import.meta.url)), derivedModule(derived));
    expect(JSON.parse(JSON.stringify(DerivedData))).toEqual(JSON.parse(JSON.stringify(derived)));
  });

  it('places 1EHZ on a seed row of high identity, the seed’s pairs landing on its stems', () => {
    expect(derived.identity).toBeGreaterThan(0.7);
    expect(derived.columns.filter((c) => c >= 0).length).toBeGreaterThan(60);
    // SS_cons projected onto 1EHZ: its pairs land on the structure’s pairs
    const read = readStockholm(text(SeedFile));
    const at = new Map(derived.columns.map((c, i) => [c, i] as const));
    const projected = read
      .pairs!.filter(([a, b]) => at.has(a) && at.has(b))
      .map(([a, b]) => `${at.get(a)}-${at.get(b)}`);
    const entry = new Set(Pairs.map(([i, j]) => `${i}-${j}`));
    expect(projected.length).toBeGreaterThanOrEqual(18);
    expect(projected.filter((p) => entry.has(p)).length / projected.length).toBeGreaterThan(0.85);
  });

  it('covaries across the stems more than elsewhere', () => {
    const stems = new Set(Pairs.filter((p) => !p[3]).map(([i, j]) => `${i}-${j}`));
    const on = derived.covariation.filter(([i, j]) => stems.has(`${i}-${j}`));
    const off = derived.covariation.filter(([i, j]) => !stems.has(`${i}-${j}`));
    const mean = (v: typeof on) => v.reduce((s, e) => s + e[2], 0) / v.length;
    expect(on.length).toBeGreaterThan(10);
    expect(mean(on)).toBeGreaterThan(mean(off) * 2);
    const g = covariationGrid();
    expect(g[0 * TrnaLength + 71]).toBe(g[71 * TrnaLength + 0]);
  });

  it('the scale grid is symmetric and from a seed', () => {
    const g = scaleGrid(64);
    expect(g[3 * 64 + 40]).toBe(g[40 * 64 + 3]);
    expect(Array.from(scaleGrid(64))).toEqual(Array.from(g));
  });
});

describe('the mmCIF reader', () => {
  const cif = [
    'data_x',
    '_entity_poly.entity_id 1',
    '_entity_poly.pdbx_seq_one_letter_code',
    ';GC(2MG)A',
    ';',
    '_entity_poly.pdbx_seq_one_letter_code_can  GCGA',
    '#',
    'loop_',
    '_t.a',
    '_t.b',
    '_t.c',
    "1 'x y' ?",
    '2 "p q"',
    ';a block',
    ';',
    '#',
  ].join('\n');

  it('reads a category of items and a loop, quoted values and text blocks', () => {
    const e = cifRows(cif, '_entity_poly')[0];
    expect(e.pdbx_seq_one_letter_code_can).toBe('GCGA');
    expect(e.pdbx_seq_one_letter_code).toBe('GC(2MG)A');
    expect(cifRows(cif, '_t')).toEqual([
      { a: '1', b: 'x y', c: '?' },
      { a: '2', b: 'p q', c: 'a block' },
    ]);
    expect(cifRows(cif, '_none')).toEqual([]);
  });

  it('classes a pair by its Saenger type', () => {
    expect([19, 20, 28, 22, '?'].map((c) => kindOfSaenger(String(c)))).toEqual([
      'canonical',
      'canonical',
      'wobble',
      'non-canonical',
      'non-canonical',
    ]);
  });
});
