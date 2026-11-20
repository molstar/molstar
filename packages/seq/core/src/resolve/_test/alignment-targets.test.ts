/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Pairs, ranges and maps on an alignment's columns (_spec/README.md "Language", Systems; "Resolve", Systems, sources,
// annotations): `pairs.a`/`b`, `ranges.target` and `map.to` name an alignment, its columns system its address; a
// family's consensus pairs held on the columns, a chain mapped onto the seed row it is closest to (a segments map
// from the gapped row), what is on the chain projected onto the columns and back.

import { applyPatch } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { checkSeqSpec, createSeqSpec, SeqLanguage, SeqValidateOptions } from '@molstar/seqx-core/lang';
import { checkSeqChanges } from '@molstar/seqx-core/lang/changes';
import { partOf, rangesLocation } from '@molstar/seqx-core/query';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import { resolveSpec, rowsSystemOf, sourceLengths, type StaticModel } from '../model.js';

/** A small seed: three rows over 12 columns; `b` is the row the chain is closest to. */
const Seed = { label: ['a', 'b', 'c'], letters: ['GCG-AUUUAGCU', 'GCGGA-UU.GCU', 'GC-GAUUUAGC-'] };
/** The consensus pairs on the columns (an `SS_cons`'s): [0, 11], [1, 10], [2, 9]. */
const Consensus = { a_start: [0, 1, 2], a_end: [1, 2, 3], b_start: [11, 10, 9], b_end: [12, 11, 10] };
/** Row `b`'s residues, ungapped: the chain. */
const Chain = Seed.letters[1].replace(/[-.]/g, '');

/** The segments of a gapped row onto its columns: residue k of the row at the column holding its k-th letter. */
function rowSegments(gapped: string) {
  const from: number[] = [];
  const to: number[] = [];
  const length: number[] = [];
  let k = 0;
  for (let c = 0; c < gapped.length; c++) {
    if (gapped[c] === '-' || gapped[c] === '.') continue;
    const last = from.length - 1;
    if (last >= 0 && from[last] + length[last] === k && to[last] + length[last] === c) length[last]++;
    else {
      from.push(k);
      to.push(c);
      length.push(1);
    }
    k++;
  }
  return { from, to, length };
}

function spec(o: { map?: boolean; past?: boolean; values?: boolean } = {}): Node {
  const b = createSeqSpec({ title: 'alignment targets' });
  const seed = b.data.alignment('seed', { alphabet: 'rna', rows: Seed });
  const chain = b.data.sequence('chain', { alphabet: 'rna', letters: Chain });
  b.data.pairs('ss', {
    a: seed,
    rows: o.past ? { ...Consensus, b_end: [13, 11, 10] } : Consensus,
  });
  b.data.ranges('stems', { target: seed, rows: { start: [0, 9], end: [3, o.past ? 13 : 12] } });
  if (o.map !== false)
    b.data.map('onto', { from: chain, to: seed, method: 'segments', segments: rowSegments(Seed.letters[1]) });
  if (o.values) b.data.values('v', { target: seed, values: new Float32Array(12) } as never);
  b.views.view('v', (v) => {
    v.section('family', { of: seed }).track('arcs').layer('arcs', { type: 'pairs', data: '/data/ss' });
    v.section('chain', { of: chain }).track('stems').layer('stems', { type: 'ranges', data: '/data/stems' });
  });
  return b.build();
}

const ranges = (m: StaticModel, from: string, r: number[], to: string) =>
  Array.from(partOf(m.project(rangesLocation(from, r), to).location, to)?.ranges ?? []);

describe('an alignment named by pairs, ranges and maps (EC-48, in part: an alignment its columns)', () => {
  it('the checks take it, and bound its ends by the columns; `values` still names a system only', () => {
    expect(checkSeqSpec(spec())).toEqual([]);
    const past = checkSeqSpec(spec({ past: true })).map((i) => [i.code, i.address]);
    expect(past).toEqual([
      ['seq.lang.pairs.bounds', '/data/ss'],
      ['seq.lang.ranges.bounds', '/data/stems'],
    ]);
    const values = checkSeqSpec(spec({ values: true })).map((i) => [i.code, i.address]);
    expect(values).toEqual([['seq.lang.ref.target', '/data/v']]);
  });

  it('a change narrowing the alignment checks what names its columns again', () => {
    const root = spec();
    const narrow = Seed.letters.map((l) => l.slice(0, 10));
    const r = applyPatch(
      root,
      { ops: [{ op: 'set', at: ['data', 'seed'], key: 'rows', value: { label: Seed.label, letters: narrow } }] },
      { language: SeqLanguage, validate: SeqValidateOptions },
    );
    if (!r.ok) throw new Error(r.issues.map((i) => i.code).join());
    const codes = checkSeqChanges(root, r.root, r.changes).map((i) => [i.code, i.address]);
    expect(codes).toEqual([
      ['seq.lang.pairs.bounds', '/data/ss'],
      ['seq.lang.ranges.bounds', '/data/stems'],
    ]);
  });

  it('pairs and ranges held on the columns system, the alignment’s address', () => {
    const m = resolveSpec(spec());
    expect(m.issues).toEqual([]);
    const ss = m.annotations.get('/data/ss')!;
    expect([ss.target, ss.ySystem]).toEqual(['/data/seed', '/data/seed']);
    expect(Array.from(ss.pairs!.bStart)).toEqual([11, 10, 9]);
    expect(m.annotations.get('/data/stems')!.target).toBe('/data/seed');
    expect(m.systems.get('/data/seed')!.length).toBe(12);
    expect(m.systems.get(rowsSystemOf('/data/seed'))!.length).toBe(3);
  });

  it('a chain’s map onto its closest seed row: residues to their columns, a column the row gaps to nothing', () => {
    const m = resolveSpec(spec());
    // b: G C G G A - U U . G C U: residue 5 (U) at column 6, residue 7 (G) at column 9
    expect(ranges(m, '/data/chain', [5, 6], '/data/seed')).toEqual([6, 7]);
    expect(ranges(m, '/data/chain', [7, 8], '/data/seed')).toEqual([9, 10]);
    expect(ranges(m, '/data/chain', [0, Chain.length], '/data/seed')).toEqual([0, 5, 6, 8, 9, 12]);
    // back: the gapped columns 5 and 8 map nowhere; the stems on the columns reach the chain's residues
    expect(ranges(m, '/data/seed', [5, 6], '/data/chain')).toEqual([]);
    expect(ranges(m, '/data/seed', [9, 12], '/data/chain')).toEqual([7, 10]);
    // the path: one step, the map read forward from the chain onto the columns, its runs the row's
    const path = m.path('/data/chain', '/data/seed').path!;
    expect([path.from, path.to, path.step]).toEqual(['/data/chain', '/data/seed', 0]);
    expect(path.edges.map((e) => [e.map, e.from, e.to, e.backward])).toEqual([
      ['/data/onto', '/data/chain', '/data/seed', false],
    ]);
    expect(path.linear.runs()?.map((r) => [r.from, r.to, r.length])).toEqual([
      [0, 0, 5],
      [5, 6, 2],
      [7, 9, 3],
    ]);
  });

  it('with no map, the ranges on the columns shown on the chain’s section: unplaced, said once', () => {
    const issues = resolveViews(resolveSpec(spec({ map: false }))).flatMap((v) => v.issues);
    expect(issues.map((i) => i.code)).toEqual(['seq.view.layer.unplaced']);
    expect(issues.map((i) => i.address)).toEqual(['/views/v/chain/stems/stems']);
  });

  it('a ranges source shared by two alignments reads each one’s columns by its name, `<alignment>/columns`', () => {
    const b = createSeqSpec({ title: 'one source over two alignments’ columns' });
    const seed = b.data.alignment('seed', { alphabet: 'rna', rows: Seed });
    const full = b.data.alignment('full', { alphabet: 'rna', rows: Seed });
    const file = b.data.source('gc-file', { format: 'bed', uri: 'families.bed' });
    b.data.ranges('gc-seed', { target: seed, source: file });
    b.data.ranges('gc-full', { target: full, source: file });
    const m = resolveSpec(b.build());
    expect([m.systems.get('/data/seed')!.name, m.systems.get('/data/full')!.name]).toEqual([
      'seed/columns',
      'full/columns',
    ]);
    const access = (a: string) => m.annotations.get(a)!.access;
    expect([access('/data/gc-seed'), access('/data/gc-full')]).toEqual([
      { kind: 'source', source: '/data/gc-file', system: 'seed/columns', shape: 'ranges' },
      { kind: 'source', source: '/data/gc-file', system: 'full/columns', shape: 'ranges' },
    ]);
  });

  it('an alignment’s own source names its columns `columns`: their length read from it', () => {
    const b = createSeqSpec({ title: 'an alignment from a source' });
    const file = b.data.source('sto', { format: 'stockholm', uri: 'family.sto' });
    b.data.alignment('family', { alphabet: 'rna', source: file });
    const m = resolveSpec(b.build());
    const coverage = {
      systems: [
        { name: 'rows', length: 7 },
        { name: 'columns', length: 42 },
      ],
    };
    const lengths = sourceLengths(m, () => coverage as never);
    const again = resolveSpec(b.build(), { lengths });
    expect([
      again.systems.get('/data/family')!.length,
      again.systems.get(rowsSystemOf('/data/family'))?.length,
    ]).toEqual([42, 7]);
  });
});
