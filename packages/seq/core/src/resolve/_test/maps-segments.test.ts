/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `segments` and `structure` maps and a numbering's `offset` and `segments` methods (_spec/README.md "Resolve", Maps
// and Numberings): the same answers on every encoding, what lies past either end mapping
// nowhere (EC-21), a run wrapping a cyclic origin (EC-11), a run read downward (EC-20's mechanism), the structure's
// residue status with the unit winning on unobserved, and author numbering with gaps and insertion codes (EC-3 to
// EC-5).

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { partOf, rangesLocation, resolveAddressIn } from '@molstar/seqx-core/query';
import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { describe, expect, it, vi } from 'vitest';
import { linearsEqual } from '../linear.js';
import { structureProvenance, structureStatus } from '../maps-segments.js';
import { type ResolveConfig, resolveSpec, type StaticModel } from '../model.js';

const Encodings: ResolveConfig['encoding'][] = ['cheapest', 'runs', 'expanded'];

/** What an address names, as `[start, end, …]` on its first part. */
function ranges(root: Node, model: StaticModel, address: string): number[] | undefined {
  const r = resolveAddressIn(root, model, address);
  return r.result?.type === 'location' ? Array.from(r.result.location.parts[0]?.ranges ?? []) : undefined;
}

function segmentsSpec(segments: Record<string, unknown>, cyclic = false): Node {
  const b = createSeqSpec();
  const a = b.data.sequence('a', { alphabet: 'dna', letters: 'ACGTACGTACGTACGTACGT' });
  const t = b.data.sequence('t', { alphabet: 'dna', letters: 'ACGTACGTAC', ...(cyclic ? { cyclic } : {}) });
  b.data.map('m', { from: a, to: t, method: 'segments', segments: segments as never });
  return b.build();
}

describe('segments maps (every encoding)', () => {
  for (const encoding of Encodings) {
    const resolve = (n: Node) => resolveSpec(n, { config: { encoding } });

    it(`runs forward and back, one read downward; one run within its target is arithmetic (EC-20; ${encoding})`, () => {
      const m = resolve(segmentsSpec({ from: [0, 6], to: [2, 9], length: [3, 4], strand: ['+', '-'] }));
      expect(m.issues).toEqual([]);
      const map = m.maps.get('/data/m')!;
      expect(map.kind).toBe('segments');
      expect(map.priority).toBe(0);
      expect(map.linear.runs()).toEqual([
        { from: 0, to: 2, length: 3, dir: 1 },
        { from: 6, to: 9, length: 4, dir: -1 },
      ]);
      expect(map.linear.backward(6, 8)).toEqual([{ from: 6, to: 9, length: 2, dir: -1 }]);
      const one = resolve(segmentsSpec({ from: [4], to: [1], length: [5] })).maps.get('/data/m')!;
      if (encoding === 'cheapest') expect(one.encoding).toBe('implicit');
      expect(one.linear.runs()).toEqual([{ from: 4, to: 1, length: 5, dir: 1 }]);
    });

    it(`what lies past either sequence maps nowhere, and is no issue (EC-21; ${encoding})`, () => {
      const m = resolve(segmentsSpec({ from: [18], to: [7], length: [6] }));
      expect(m.issues).toEqual([]);
      expect(m.maps.get('/data/m')!.linear.runs()).toEqual([{ from: 18, to: 7, length: 2, dir: 1 }]);
      const none = resolve(segmentsSpec({ from: [0], to: [12], length: [3] }));
      expect(none.issues.map((i) => i.code)).toEqual(['seq.resolve.map.empty']);
    });

    it(`onto a cyclic target a run wraps, split at the origin (EC-11; ${encoding})`, () => {
      const m = resolve(segmentsSpec({ from: [0], to: [8], length: [4] }, true));
      expect(m.maps.get('/data/m')!.linear.runs()).toEqual([
        { from: 0, to: 8, length: 2, dir: 1 },
        { from: 2, to: 0, length: 2, dir: 1 },
      ]);
    });
  }

  it('every encoding gives the same map', () => {
    const spec = segmentsSpec({ from: [0, 5, 12], to: [9, 0, 6], length: [3, 4, 2], strand: ['-', '+', '+'] });
    const [a, ...rest] = Encodings.map((encoding) => resolveSpec(spec, { config: { encoding } }).maps.get('/data/m')!);
    for (const b of rest) expect(linearsEqual(a.linear, b.linear)).toBe(true);
  });
});

/** UniProt-like `up` onto unit A of an entity: a tag on the chain, a conflict, an unobserved loop. */
function structureSpec(unobserved: { start: number[]; end: number[] }): Node {
  const b = createSeqSpec();
  const up = b.data.sequence('up', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQ' });
  const ent = b.data.sequence('ent', { alphabet: 'protein', letters: 'GSHMEEPQSDPSVEPP' });
  let a = '';
  b.data.complex('c', {}, (c) => {
    a = c.unit('A', { of: ent, chain: 'A', observed: { start: [2, 12], end: [10, 16] } });
  });
  b.data.map('sifts', {
    from: up,
    to: a,
    method: 'structure',
    // the tag GSH is unpaired; M from the tag is paired by an alignment
    segments: { from: [0, 4], to: [3, 7], length: [4, 9], provenance: ['alignment', 'sifts'] },
    status: {
      start: [6, ...unobserved.start],
      end: [7, ...unobserved.end],
      status: ['conflict', ...unobserved.start.map(() => 'unobserved')],
    },
    priority: 2,
  });
  return b.build();
}

describe('structure maps: provenance and residue status', () => {
  it('pairs, conflicts and unpaired residues; each run’s provenance; the unit’s observed taken', () => {
    const m = resolveSpec(structureSpec({ start: [10], end: [12] }));
    expect(m.issues).toEqual([]);
    const map = m.maps.get('/data/sifts')!;
    expect([map.kind, map.priority]).toEqual(['structure', 2]);
    const unit = m.molecular.units.get('/data/c/A')!;
    const status = Array.from({ length: 16 }, (_, k) => structureStatus(map, k, unit.observed));
    expect(status).toEqual([
      'unobserved',
      'unobserved',
      'unpaired',
      'paired',
      'paired',
      'paired',
      'conflict',
      'paired',
      'paired',
      'paired',
      'unobserved',
      'unobserved',
      'paired',
      'paired',
      'paired',
      'paired',
    ]);
    // the map's own word where no unit is given
    expect(structureStatus(map, 10)).toBe('unobserved');
    expect(structureStatus(map, 0)).toBe('unpaired');
    expect([0, 5, 14].map((k) => structureProvenance(map.structure!, k))).toEqual(['alignment', 'sifts', undefined]);
    // a sequence position on the chain through the map: up 5 (S) is the chain's 8
    const r = m.project(rangesLocation('/data/up', [5, 6]), '/data/c/A').location;
    expect(Array.from(partOf(r, '/data/c/A')?.ranges ?? [])).toEqual([8, 9]);
  });

  it('a run’s provenance is its own row’s, past a row that holds nothing (an unchecked spec)', () => {
    const b = createSeqSpec();
    const up = b.data.sequence('up', { alphabet: 'protein', letters: 'MEEPQSDPSVEPPLSQ' });
    const ent = b.data.sequence('ent', { alphabet: 'protein', letters: 'GSHMEEPQSDPSVEPP' });
    b.data.map('m', {
      from: up,
      to: ent,
      method: 'structure',
      segments: { from: [0, 0, 4], to: [0, 3, 7], length: [0, 4, 9], provenance: ['none', 'alignment', 'sifts'] },
    });
    const map = resolveSpec(b.build()).maps.get('/data/m')!;
    expect([0, 5].map((k) => structureProvenance(map.structure!, k))).toEqual(['alignment', 'sifts']);
  });

  it('reads a unit’s observed runs once, not per residue asked', () => {
    const m = resolveSpec(structureSpec({ start: [10], end: [12] }));
    const map = m.maps.get('/data/sifts')!;
    const observed = Int32Array.of(2, 10, 12, 16);
    const normalize = vi.spyOn(HalfOpenRanges, 'normalize');
    for (let k = 0; k < 16; k++) structureStatus(map, k, observed);
    expect(normalize).toHaveBeenCalledTimes(1);
  });

  it('a map’s unobserved residues that are not the unit’s are reported once; the unit wins', () => {
    const m = resolveSpec(structureSpec({ start: [9], end: [12] }));
    expect(m.issues.map((i) => i.code)).toEqual(['seq.resolve.map.observed']);
    const map = m.maps.get('/data/sifts')!;
    expect(structureStatus(map, 9, m.molecular.units.get('/data/c/A')!.observed)).toBe('paired');
  });
});

describe('numberings: offset and segments', () => {
  function numbered(): Node {
    const b = createSeqSpec();
    b.data.sequence('s', { alphabet: 'protein', letters: 'MEEPQSDPSVEP' }, (s) => {
      s.numbering('uniprot', { method: 'offset', by: 94, scheme: 'uniprot' });
      s.numbering('auth', {
        method: 'segments',
        scheme: 'auth',
        // 1–4, a gap, 52 and 52A, then 60–64
        segments: {
          start: [0, 5, 6, 7],
          end: [4, 6, 7, 12],
          first: [1, 52, 52, 60],
          insertion: [null, null, 'A', null],
        },
      });
    });
    return b.build();
  }

  it('an offset numbers a position plus `by`, answering a label without a scan (EC-3)', () => {
    const m = resolveSpec(numbered());
    expect(m.issues).toEqual([]);
    const n = m.numberingOf('/data/s', 'uniprot')!.numbering!;
    expect([n.labelAt(0), n.labelAt(11), n.labelAt(12)]).toEqual(['94', '105', undefined]);
    expect([n.indexOf('100'), n.indexOf('93')]).toEqual([6, -1]);
    expect(n.snap?.('99', 'start')).toBe(5);
    expect(ranges(numbered(), m, '/data/s/(at :label 100 :numbering uniprot)')).toEqual([6, 7]);
  });

  it('segments number runs, leave a gap unlabelled and carry insertion codes (EC-4, EC-5)', () => {
    const m = resolveSpec(numbered());
    const n = m.numberingOf('/data/s', 'auth')!.numbering!;
    expect(Array.from({ length: 12 }, (_, i) => n.labelAt(i) ?? null)).toEqual([
      '1',
      '2',
      '3',
      '4',
      null,
      '52',
      '52A',
      '60',
      '61',
      '62',
      '63',
      '64',
    ]);
    expect([n.indexOf('52'), n.indexOf('52A'), n.indexOf('52a'), n.indexOf('5'), n.indexOf('x')]).toEqual([
      5, 6, 6, -1, -1,
    ]);
    expect(ranges(numbered(), m, '/data/s/(range :from 52A :to 61 :numbering auth)')).toEqual([6, 9]);
    expect(ranges(numbered(), m, '/data/s/(at :label 60 :numbering auth)')).toEqual([7, 8]);
  });
});
