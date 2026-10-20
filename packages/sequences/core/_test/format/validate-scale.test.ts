/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Validation at genome scale and its cost (V5): a table's `target` column is a ref per row, so a dense values
// table over 3M positions names a chromosome 3M times. The check reads each cell once, records one ref per distinct
// value, and a snapshot re-checks only the nodes it touched (each node's own checks memoized by its params object).
// Times are compared between two sizes, the best of five runs each, never against a fixed bound (the heavy lane: they
// read the clock); where the work can be counted (cell reads), it is counted instead, in the quick lane.

import { SeqSets } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import {
  int,
  RequiredField,
  register,
  type SeqViewSpecDoc,
  SimpleParamsSchema,
  validate,
} from '../../src/format/index.js';
import { kernel } from './fixtures.js';
import { CHROMOSOMES, codedDoc, countingCodes, freshParams, gcNode } from './table-fixtures.js';

type J = any;
const reg = kernel();
const GC = '/root/children/0/children/24';
const best = (run: () => void) => {
  let t = Number.POSITIVE_INFINITY;
  for (let k = 0; k < 5; k++) {
    const t0 = performance.now();
    run();
    t = Math.min(t, performance.now() - t0);
  }
  return t;
};
const brief = (d: SeqViewSpecDoc, o = {}) => validate(d, reg, o).map((i) => `${i.code} ${i.path}`);

describe('validation scales with table rows', () => {
  it('200,000 target refs from one annotation validate, coded or as strings, in linear time', {
    tags: ['heavy'],
  }, () => {
    for (const plain of [false, true]) {
      const doc = (n: number) => {
        const d = codedDoc(n);
        // The uncoded form: a string per row.
        if (plain) gcNode(d).params.value.target = Array.from({ length: n }, (_, i) => CHROMOSOMES[((i * 24) / n) | 0]);
        return d;
      };
      const small = doc(50_000);
      const large = doc(200_000);
      expect(validate(large, reg)).toEqual([]);
      const ratio = best(() => validate(freshParams(large), reg)) / best(() => validate(freshParams(small), reg));
      // Quadratic in refs (the old ref graph), 4× the rows took 16×; linear, about 4×.
      expect(ratio).toBeLessThan(8);
    }
  });

  it('the 3M-row coded genome table validates in linear time (V5)', { tags: ['heavy'] }, () => {
    const small = codedDoc(750_000);
    const large = codedDoc(3_000_000);
    expect(validate(large, reg)).toEqual([]);
    const ratio = best(() => validate(freshParams(large), reg)) / best(() => validate(freshParams(small), reg));
    expect(ratio).toBeLessThan(8);
  });

  it('10 snapshots cost one validation, not eleven: counted in cell reads, not timed', () => {
    // VW2 timed this (a ratio under 2), which failed on a loaded machine; the work is the cells read.
    const n = 100_000;
    const reads = (snapshots: number) => {
      const { codes, reads } = countingCodes(Array.from({ length: n }, (_, i) => Math.floor((i * 24) / n)));
      expect(validate(codedDoc(n, { snapshots, codes }), reg)).toEqual([]);
      return reads();
    };
    // Each cell read once by the checks and once by the envelope's depth bound (CM-03), which remembers the list.
    expect([reads(0), reads(10)]).toEqual([2 * n, 2 * n]);
  });
});

describe('the per-node memo', () => {
  it('reads each cell once across 10 snapshots, and not at all on revalidation', () => {
    const n = 1000;
    const { codes, reads } = countingCodes(Array.from({ length: n }, (_, i) => i % 24));
    const d = codedDoc(n, { snapshots: 10, codes });
    expect(validate(d, reg)).toEqual([]);
    // Once by the checks, once by the depth bound (CM-03); both remember it.
    expect(reads()).toBe(2 * n);
    expect(validate(d, reg)).toEqual([]);
    expect(reads()).toBe(2 * n);
  });

  it('a registry version bump invalidates it', () => {
    const r = kernel();
    const d = codedDoc(100);
    gcNode(d).custom = { 'org.example': { score: 'high' } };
    expect(validate(d, r)).toEqual([]);
    const schema = SimpleParamsSchema({ score: RequiredField(int, 'a score') });
    register(r, { custom: [{ namespace: 'org.example', kind: 'annotation', schema }] });
    expect(validate(d, r).map((i) => `${i.code} ${i.path}`)).toEqual([
      `svs.custom.invalid ${GC}/custom/org.example/score`,
    ]);
  });

  it('an issue moves with its node: a sibling inserted before it shifts the path', () => {
    const d: J = codedDoc(100);
    gcNode(d).params.label = 5;
    expect(brief(d)).toEqual([`svs.params.type ${GC}/params/label`]);
    const data = d.root.children[0];
    const extra = { kind: 'sequence', ref: 'chrM', params: { alphabet: 'dna', value: 'ACGT' } };
    const shifted = {
      ...d,
      root: { ...d.root, children: [{ ...data, children: [extra, ...data.children] }, d.root.children[1]] },
    };
    expect(brief(shifted)).toEqual(['svs.params.type /root/children/0/children/25/params/label']);
    expect(brief(d)).toEqual([`svs.params.type ${GC}/params/label`]);
  });

  it('one glyph key is an error under 1.0 and a warning under 1.1, from the same params object', () => {
    const d: J = codedDoc(100);
    const layer = d.root.children[1].children[0].children[24].children[0];
    layer.params.glyph = { theme: 'uniform', params: { glyph: { char: 'x', zoom: 2 } } };
    const at = '/root/children/1/children/0/children/24/children/0/params/glyph/params/glyph/zoom';
    const run = (version: string) => validate({ ...d, version }, reg).map((i) => `${i.code} ${i.severity} ${i.path}`);
    expect(run('1.0.0')).toEqual([`svs.glyph.invalid error ${at}`]);
    expect(run('1.1.0')).toEqual(['svs.version.newer_minor warning /version', `svs.version.newer_minor warning ${at}`]);
    expect(run('1.0.0')).toEqual([`svs.glyph.invalid error ${at}`]);
  });

  it('refs to a set are resolved on every call: a new set lacking a molecule reports it, the old one clears it', () => {
    const d: J = codedDoc(240);
    d.root.children[0].children = [gcNode(d)];
    const full = SeqSets.build({ molecules: CHROMOSOMES.map((key) => ({ key, residues: 'ACGT' })) }).set;
    const short = SeqSets.build({ molecules: CHROMOSOMES.slice(0, 23).map((key) => ({ key, residues: 'ACGT' })) }).set;
    expect(validate(d, reg, { set: full })).toEqual([]);
    const issues = validate(d, reg, { set: short });
    expect(issues.map((i) => i.code)).toEqual(['svs.ref.unresolved', 'svs.ref.unresolved']);
    expect(issues.map((i) => i.path)).toContain('/root/children/0/children/0/params/value/target/dict/23');
    expect(validate(d, reg, { set: full })).toEqual([]);
  });
});

describe('repeated issues', () => {
  it('one unknown chrX in a million rows is one issue with its count', () => {
    const n = 1_000_000;
    const codes = new Uint8Array(n);
    codes.fill(24, n / 2);
    const d = codedDoc(n, { codes, dict: [...CHROMOSOMES, 'chrX'] });
    const issues = validate(d, reg);
    expect(issues.map((i) => `${i.code} ${i.path}`)).toEqual([`svs.ref.unresolved ${GC}/params/value/target/dict/24`]);
    expect(issues[0].message).toContain('(500000 rows)');
  });

  it('as strings too, and a bad cell repeated down a column is one issue', () => {
    const d: J = codedDoc(1000);
    const value = gcNode(d).params.value;
    value.target = Array.from({ length: 1000 }, (_, i) => (i < 10 ? 'chr1' : 'chrX'));
    value.value = Array.from({ length: 1000 }, (_, i) => (i < 3 ? 0.5 : 'high'));
    const issues = validate(d, reg);
    expect(issues.map((i) => `${i.code} ${i.path}`)).toEqual([
      `svs.columns.type ${GC}/params/value/value/3`,
      `svs.ref.unresolved ${GC}/params/value/target/10`,
    ]);
    expect(issues[0].message).toContain('(997 rows)');
    expect(issues[1].message).toContain('(990 rows)');
  });
});
