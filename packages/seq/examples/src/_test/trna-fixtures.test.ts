/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { decodeSeq } from '@molstar/seqx-core/data';
import { readStockholm } from '@molstar/seqx-core/data/alignment-readers';
import { describe, expect, it } from 'vitest';
import { PackageFolder } from '../node-assets.js';
import { DerivedData } from '../trna/derived-data.js';
import { derive, mutualInformation } from '../trna/derive.js';
import { writeDerived } from '../trna/generated.js';
import {
  covariationGrid,
  ModifiedCompounds,
  ModifiedResidues,
  MoleculeFile,
  PairsFile,
  Pairs,
  SeedFile,
  scaleGrid,
  SystemName,
  TrnaLength,
  TrnaLetters,
} from '../trna/data.js';
import { readCif } from '../../scripts/cif.js';
import { kindOfSaenger } from '../trna/record-trim.js';

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
    // the sequence's `modified` table the views write is the file's: each component at its position, its parent the letter
    expect(ModifiedResidues.map((m) => m[0])).toEqual(Array.from(r.start as Int32Array));
    expect(ModifiedResidues.map((m) => m[1])).toEqual(field(r, 'label'));
    for (const [p, code] of ModifiedResidues) expect(ModifiedCompounds[code].parent).toBe(letters[p]);
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

  it('sources.json: every file’s bytes and checksum are the file’s, the licences as the section says', () => {
    const rows = JSON.parse(text('fixtures/trna/sources.json')) as {
      file: string;
      entry: string;
      bytesFixture: number;
      sha256Fixture: string;
      licence: string[];
    }[];
    expect(rows.map((r) => r.file).sort()).toEqual([
      'src/trna/data.ts#scaleGrid',
      'src/trna/derived-data.ts',
      'trna/1EHZ-pairs.seq',
      'trna/1EHZ.seq',
      'trna/RF00005.sto',
    ]);
    for (const r of rows) {
      expect(r.licence).toEqual([r.file.startsWith('src/') ? 'LicenseRef-seqx-generated' : 'CC0-1.0']);
      if (!r.bytesFixture) continue;
      const t = text(r.file.startsWith('src/') ? r.file : `fixtures/${r.file}`);
      expect(new TextEncoder().encode(t).byteLength).toBe(r.bytesFixture);
      expect(createHash('sha256').update(t).digest('hex')).toBe(r.sha256Fixture);
    }
    // the release is read from Rfam's README, not written in the script
    expect(rows.find((r) => r.file === 'trna/RF00005.sto')!.entry).toMatch(/release \d+\.\d+/);
  });
});

describe('what is computed from them', () => {
  const letters = ranges(MoleculeFile).tables.letters.letters as string;
  // computed once for the file, not once a test: 954 alignments and 2,850 pairs of columns, measured at 150 ms
  const derived = derive(text(SeedFile), letters);

  it('is what is committed', () => {
    // written as the formatter does: a regeneration leaves `format:check` clean
    if (update)
      writeDerived(fileURLToPath(new URL('../trna/derived-data.ts', import.meta.url)), text(SeedFile), letters);
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

  it('covaries across the stems more than elsewhere: the grid’s exact cells', () => {
    const g = covariationGrid();
    const at = (i: number, j: number) => g[i * TrnaLength + j];
    const stems = Pairs.filter((p) => !p[3] && p[2] !== 'N');
    // a stem pair's cell holds its value, symmetric, and the entries are the upper triangle's (i < j)
    for (const [i, j] of stems) expect(at(i, j)).toBe(at(j, i));
    expect(derived.covariation.every(([i, j]) => i < j)).toBe(true);
    // the cells of the stems above the rest: the lowest of those with evidence beats the median of the others
    const onStem = new Set(stems.map(([i, j]) => `${i}-${j}`));
    const withValue = stems.map(([i, j]) => at(i, j)).filter((v) => !Number.isNaN(v));
    expect(withValue.length).toBeGreaterThanOrEqual(15);
    const others = derived.covariation
      .filter(([i, j]) => !onStem.has(`${i}-${j}`))
      .map((e) => e[2])
      .sort((a, b) => a - b);
    expect(withValue.sort((a, b) => a - b)[Math.floor(withValue.length / 4)]).toBeGreaterThan(
      others[Math.floor(others.length / 2)],
    );
    // the best cell of the grid is a stem’s
    let best = [0, 0];
    for (const [i, j, v] of derived.covariation)
      if (v > at(best[0], best[1]) || Number.isNaN(at(best[0], best[1]))) best = [i, j];
    expect(onStem.has(`${best[0]}-${best[1]}`)).toBe(true);
  });

  it('the scale grid is symmetric and from a seed', () => {
    const g = scaleGrid(64);
    expect(g[3 * 64 + 40]).toBe(g[40 * 64 + 3]);
    expect(Array.from(scaleGrid(64))).toEqual(Array.from(g));
  });
});

describe('mutual information', () => {
  const n = 100;
  const column = (f: (k: number) => number) => Uint8Array.from({ length: n }, (_, k) => f(k));

  it('a perfectly covarying pair is 2 bits (never past the columns’ entropy); independent columns 0', () => {
    const a = column((k) => k % 4);
    expect(
      mutualInformation(
        a,
        column((k) => ((k % 4) + 1) % 4),
      ),
    ).toBeCloseTo(2, 10);
    // each letter of a meets each letter of b equally: no information, and the correction takes it below zero
    expect(
      mutualInformation(
        a,
        column((k) => Math.floor(k / 4) % 4),
      ),
    ).toBe(0);
  });

  it('a pair that agrees for three rows in four: the plug-in value less the Miller–Madow bias', () => {
    const a = column((k) => k % 4);
    const b = column((k) => (k % 4 === 0 && Math.floor(k / 4) % 4 === 0 ? 1 : k % 4));
    // joint: 100 rows, the diagonal 25 each but one row in four of letter 0 moved: counts by hand
    const joint = new Map<number, number>();
    for (let k = 0; k < n; k++) joint.set(a[k] * 4 + b[k], (joint.get(a[k] * 4 + b[k]) ?? 0) + 1);
    const pa = [0, 1, 2, 3].map((x) =>
      [...joint].filter(([c]) => Math.floor(c / 4) === x).reduce((s, [, v]) => s + v, 0),
    );
    const pb = [0, 1, 2, 3].map((y) => [...joint].filter(([c]) => c % 4 === y).reduce((s, [, v]) => s + v, 0));
    let mi = 0;
    for (const [c, v] of joint) mi += (v / n) * Math.log2((v * n) / (pa[Math.floor(c / 4)] * pb[c % 4]));
    const used = (v: number[]) => v.filter((x) => x > 0).length;
    const bias = (joint.size - used(pa) - used(pb) + 1) / (2 * n * Math.LN2);
    expect(mutualInformation(a, b)).toBeCloseTo(Math.max(0, mi - bias), 10);
    expect(mutualInformation(a, b)).toBeLessThan(2);
  });

  it('under 30 rows with both letters is NaN; a gap or another letter leaves the row out', () => {
    expect(mutualInformation(new Uint8Array(29), new Uint8Array(29))).toBeNaN();
    const a = column((k) => (k < 40 ? 4 : k % 4));
    expect(
      mutualInformation(
        a,
        column((k) => ((k % 4) + 1) % 4),
      ),
    ).toBeCloseTo(2, 10);
  });
});

describe('the record scripts’ mmCIF reader (`@molstar/io`’s parser)', () => {
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

  it('reads a category of items and a loop, quoted values and text blocks, `.` and `?` as written', async () => {
    const parsed = await readCif(cif);
    const cifRows = (_: string, c: string) => parsed.rows(c);
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
