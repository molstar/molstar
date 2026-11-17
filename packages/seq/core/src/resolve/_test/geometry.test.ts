/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The static model of p53 on chr17 (its real coding geometry), every encoding through the model's configuration,
// closure along the map itself, path issues and their places, sources named and bound, lengths not known yet.

import type { Node } from '@molstar/core/repr/tree';
import type { FormatRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { evaluateQuery, type Location, rangesLocation } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import { GeneticCodes } from '../codon.js';
import { createResolver, resolveSpec } from '../model.js';
import { p53GeometrySpec } from './p53-geometry.js';

const Fasta: FormatRuntime = { name: 'fasta', formats: ['fasta'], open: () => Promise.reject(new Error('no')) };
const Spec = p53GeometrySpec();
const ranges = (l: Location) => l.parts.map((p) => [p.system, [...p.ranges]]);
const codes = (issues: readonly { code: string }[]) => issues.map((i) => i.code);

describe('static model: p53 on chr17', () => {
  const model = resolveSpec(Spec, { runtimes: [Fasta] });

  it('reads the gene by its record’s name in its source', () => {
    expect(model.issues).toEqual([]);
    const chr17 = model.systems.get('/data/chr17')!;
    expect([chr17.name, chr17.length, chr17.access]).toEqual([
      '17',
      undefined,
      { kind: 'source', source: '/data/grch38', system: '17', shape: 'letters' },
    ]);
    expect(model.maps.get('/data/p53-on-chr17')!.code).toBe(GeneticCodes[1]);
  });

  it('carries the design’s query onto the gene: the DNA-binding domain over exons 4 to 8', () => {
    const q = evaluateQuery(Spec, model, '/data/p53/(range :from 94 :to 292)', { limit: 10 });
    const dbd = q.result?.type === 'location' ? q.result.location : undefined;
    expect(dbd && ranges(dbd)).toEqual([['/data/p53', [93, 292]]]);
    const onGene = model.project(dbd!, '/data/chr17').location.parts[0].ranges;
    expect([onGene.length / 2, onGene[0], onGene[onGene.length - 1]]).toEqual([5, 7673743, 7676089]);
  });

  it('answers alike whatever encoding the configuration holds its maps in', () => {
    const loc = rangesLocation('/data/p53', [20, 30, 380, 393]);
    const want = ranges(model.project(loc, '/data/chr17').location);
    for (const encoding of ['runs', 'expanded'] as const) {
      const m = resolveSpec(Spec, { runtimes: [Fasta], config: { encoding } });
      expect(m.maps.get('/data/p53-on-chr17')!.encoding).toBe(encoding);
      expect(ranges(m.project(loc, '/data/chr17').location)).toEqual(want);
      const back = m.project(rangesLocation('/data/chr17', [7676520, 7676522]), '/data/p53').location;
      expect(ranges(back)).toEqual([['/data/p53', [24, 25]]]);
    }
  });
});

describe('static model: composition through an expanded map', () => {
  it('a codon path through an expanded offset composes as through its runs', () => {
    const b = createSeqSpec();
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MKW' });
    const g = b.data.sequence('g', { alphabet: 'dna', letters: 'A'.repeat(20) });
    const chr = b.data.sequence('chr', { alphabet: 'dna', letters: 'A'.repeat(100) });
    b.data.map('codon', { from: p, to: g, method: 'codon', exons: { start: [2, 12], end: [6, 17] }, strand: '-' });
    b.data.map('place', { from: g, to: chr, method: 'offset', by: 40 });
    const spec = b.build();
    const answers = (['cheapest', 'runs', 'expanded'] as const).map((encoding) => {
      const m = resolveSpec(spec, { config: { encoding } });
      const path = m.path('/data/p', '/data/chr').path!;
      return [ranges(m.project(rangesLocation('/data/p', [1, 2]), '/data/chr').location), path.linear.runs()];
    });
    expect(answers[1]).toEqual(answers[0]);
    expect(answers[2]).toEqual(answers[0]);
    expect(answers[0][0]).toEqual([['/data/chr', [45, 46, 52, 54]]]);
  });
});

describe('static model: closure along the map itself (EC-30)', () => {
  it('closes through each codon map, though no single path joins the gene to either protein', () => {
    const b = createSeqSpec();
    const g = b.data.sequence('g', { alphabet: 'dna', letters: 'A'.repeat(30) });
    const p1 = b.data.sequence('p1', { alphabet: 'protein', letters: 'MKW' });
    const p2 = b.data.sequence('p2', { alphabet: 'protein', letters: 'MKW' });
    b.data.map('m1', { from: p1, to: g, method: 'codon', exons: { start: [0], end: [9] }, strand: '+' });
    // two maps of p2 onto the gene, differently: the path between them is ambiguous
    b.data.map('m2', { from: p2, to: g, method: 'codon', exons: { start: [10], end: [19] }, strand: '+' });
    b.data.map('m3', { from: p2, to: g, method: 'codon', exons: { start: [20], end: [29] }, strand: '+' });
    const m = resolveSpec(b.build());
    expect(codes([m.path('/data/g', '/data/p2').issue!])).toEqual(['seq.resolve.path.ambiguous']);
    const closed = m.close(rangesLocation('/data/g', [4, 5, 25, 26]));
    expect(ranges(closed)).toEqual([
      ['/data/g', [3, 6, 23, 26]],
      ['/data/p1', [1, 2]],
      ['/data/p2', [1, 2]],
    ]);
  });
});

describe('static model: path issues', () => {
  /** Systems a and b joined by `n` equal two-step paths through middles. */
  function fan(n: number): Node {
    const b = createSeqSpec();
    const a = b.data.sequence('a', { alphabet: 'protein', letters: 'MKMK' });
    const z = b.data.sequence('z', { alphabet: 'protein', letters: 'MKMK' });
    for (let i = 0; i < n; i++) {
      const mid = b.data.sequence(`m${i}`, { alphabet: 'protein', letters: 'MKMK' });
      b.data.map(`a${i}`, { from: a, to: mid, method: 'offset', by: 0 });
      b.data.map(`z${i}`, { from: mid, to: z, method: 'offset', by: 0 });
    }
    return b.build();
  }

  it('names both systems, the first where the path starts, read either way', () => {
    const b = createSeqSpec();
    b.data.sequence('a', { alphabet: 'protein', letters: 'MK' });
    b.data.sequence('b', { alphabet: 'protein', letters: 'MK' });
    const m = resolveSpec(b.build());
    const ab = m.path('/data/a', '/data/b').issue!;
    const ba = m.path('/data/b', '/data/a').issue!;
    expect([ab.code, ab.address, ab.places]).toEqual([
      'seq.resolve.path.none',
      '/data/a',
      [{ address: '/data/a' }, { address: '/data/b' }],
    ]);
    expect([ba.address, ba.places]).toEqual(['/data/b', [{ address: '/data/b' }, { address: '/data/a' }]]);
  });

  it('reports more candidates than are compared as its own code, a bound of the configuration', () => {
    expect(resolveSpec(fan(3)).path('/data/a', '/data/z').path).toBeDefined();
    const limited = resolveSpec(fan(3), { config: { maxPathCandidates: 2 } });
    expect(codes([limited.path('/data/a', '/data/z').issue!])).toEqual(['seq.resolve.path.limit']);
    // the reverse, found from the forward search, names its own start
    expect(limited.path('/data/z', '/data/a').issue!.address).toBe('/data/z');
  });

  it('carries paths across resolves while the maps their search read are the same', () => {
    const resolver = createResolver();
    const a = resolver.resolve(fan(2));
    const pa = a.path('/data/a', '/data/z').path!;
    // another sequence with no map: the search read nothing of it, so the path is kept
    const spec = fan(2);
    const b = resolver.resolve({
      ...spec,
      children: spec.children!.map((c) =>
        c.kind === 'data'
          ? {
              ...c,
              children: [...c.children!, { kind: 'sequence', id: 'lone', params: { alphabet: 'dna', letters: 'A' } }],
            }
          : c,
      ),
    });
    expect(b).not.toBe(a);
    expect(b.path('/data/a', '/data/z').path).toBe(pa);
    // a third way from a: searched again
    const c = resolver.resolve(fan(3));
    expect(c.path('/data/a', '/data/z').path).not.toBe(pa);
  });
});

describe('static model: sources and lengths not known yet', () => {
  it('a source with no format is reported and unavailable', () => {
    const spec: Node = {
      kind: 'spec',
      params: { version: '0.0.0' },
      children: [{ kind: 'data', id: 'data', children: [{ kind: 'source', id: 's', params: { uri: 'x' } }] }],
    };
    const m = resolveSpec(spec, { runtimes: [Fasta] });
    expect(codes(m.issues)).toEqual(['seq.resolve.source.format']);
    expect(m.sources.get('/data/s')).toEqual({
      address: '/data/s',
      entry: { address: '/data/s', format: '', uri: 'x' },
    });
  });

  it('a genetic code not built is reported', () => {
    const b = createSeqSpec();
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'M' });
    const g = b.data.sequence('g', { alphabet: 'dna', letters: 'ATG' });
    b.data.map('m', { from: p, to: g, method: 'codon', exons: { start: [0], end: [3] }, strand: '+', table: 4 });
    b.data.map('mito', { from: p, to: g, method: 'codon', exons: { start: [0], end: [3] }, strand: '+', table: 2 });
    const m = resolveSpec(b.build());
    expect(codes(m.issues)).toEqual(['seq.resolve.map.table']);
    expect(m.maps.get('/data/m')!.code).toBeUndefined();
    expect(m.maps.get('/data/mito')!.code).toBe(GeneticCodes[2]);
  });

  it('a row or an exon over the origin of a cyclic sequence of unknown length waits for it', () => {
    const b = createSeqSpec();
    const src = b.data.source('src', { format: 'fasta', uri: 'p.fa' });
    const pl = b.data.sequence('plasmid', { alphabet: 'dna', source: src, cyclic: true });
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MKW' });
    b.data.ranges('f', { target: pl, rows: { start: [8, 1], end: [2, 4], id: ['ori', 'x'] } });
    b.data.map('m', { from: p, to: pl, method: 'codon', exons: { start: [3, 17], end: [6, 3] }, strand: '+' });
    const m = resolveSpec(b.build(), { runtimes: [Fasta] });
    expect(codes(m.issues)).toEqual(['seq.resolve.length.unknown', 'seq.resolve.length.unknown']);
    expect(m.issues.map((i) => i.address)).toEqual(['/data/m', '/data/f']);
    expect(m.annotations.get('/data/f')!.items).toEqual([{ key: 'x', start: 1, end: 4 }]);
  });
});
