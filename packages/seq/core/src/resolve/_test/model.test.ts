/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The static model of the p53 example and of small specs: systems, maps, sources bound to runtimes, access plans, the
// query forms over it, projection both ways through the codon map on either strand, closure, caching by hashes.

import type { Node } from '@molstar/core/repr/tree';
import { type FormatRuntime, readSeq, SeqRuntime } from '@molstar/seqx-core/data';
import { checkSeqSpec, createSeqSpec } from '@molstar/seqx-core/lang';
import { evaluateQuery, type Location, rangesLocation } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import { EXONS, P53_LETTERS, p53Spec } from '../../lang/_test/p53.js';
import { createResolver, OpenLength, resolveSpec, type StaticModel } from '../model.js';

const Fasta: FormatRuntime = {
  name: 'fasta',
  formats: ['fasta'],
  open: () => Promise.reject(new Error('not in tests')),
};

const P53 = p53Spec();
const model = resolveSpec(P53, { runtimes: [Fasta] });

const ranges = (l: Location) => l.parts.map((p) => [p.system, [...p.ranges]]);
const codes = (issues: readonly { code: string }[]) => issues.map((i) => i.code);

/** The p53 example with the gene's strand, or its letters inline. */
function variant(strand: '+' | '-', geneLetters?: string): Node {
  const b = createSeqSpec({ title: 'p53' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53_LETTERS });
  const gene = geneLetters
    ? b.data.sequence('tp53', { alphabet: 'dna', letters: geneLetters })
    : b.data.sequence('tp53', { alphabet: 'dna', source: b.data.source('chr17', { format: 'fasta', uri: 'x.fa' }) });
  b.data.map('p53-on-gene', { from: p53, to: gene, method: 'codon', exons: EXONS, strand });
  return b.build();
}

describe('static model: the p53 example', () => {
  it('declares its systems: the protein inline, the gene open until its source is read', () => {
    expect(checkSeqSpec(P53).filter((i) => i.severity === 'error')).toEqual([]);
    const p53 = model.systems.get('/data/p53')!;
    expect([p53.alphabet, p53.length, p53.cyclic, p53.numberings]).toEqual([
      'protein',
      30,
      false,
      ['index', 'ordinal'],
    ]);
    const gene = model.systems.get('/data/tp53')!;
    expect([gene.alphabet, gene.length, gene.source]).toEqual(['dna', undefined, '/data/chr17']);
    expect(model.issues).toEqual([]);
  });

  it('resolves the codon map to its runs, minus strand', () => {
    const m = model.maps.get('/data/p53-on-gene')!;
    expect([m.kind, m.from, m.to, m.encoding, m.step]).toEqual(['codon', '/data/p53', '/data/tp53', 'runs', 1]);
    expect(m.linear.runs()).toEqual([
      { from: 0, to: 199, length: 40, dir: -1 },
      { from: 40, to: 149, length: 50, dir: -1 },
    ]);
  });

  it('binds each source to a runtime by its format, and plans each node’s reads', () => {
    const chr17 = model.sources.get('/data/chr17')!;
    expect(chr17.runtime).toBe(Fasta);
    expect(chr17.entry).toEqual({
      address: '/data/chr17',
      format: 'fasta',
      uri: 'data/chr17-region.fa',
      index: 'data/chr17-region.fa.fai',
    });
    expect(model.systems.get('/data/tp53')!.access).toEqual({
      kind: 'source',
      source: '/data/chr17',
      system: 'tp53',
      shape: 'letters',
    });
    // inline data is a seq source of one tile over the whole system
    const domains = model.annotations.get('/data/domains')!;
    expect([domains.shape, domains.target, domains.fields]).toEqual(['ranges', '/data/p53', ['label']]);
    const plan = domains.access!;
    expect(plan.kind === 'inline' && [plan.runtime, plan.system, plan.shape]).toEqual([
      SeqRuntime,
      '/data/p53',
      'ranges',
    ]);
    const read = plan.kind === 'inline' ? readSeq(plan.source.value!) : undefined;
    expect(read?.data?.pairs.map((p) => [p.request.start, p.request.end])).toEqual([[0, 30]]);
    const values = model.annotations.get('/data/conservation')!;
    expect([values.shape, values.access?.kind, values.access?.shape]).toEqual(['values', 'inline', 'values']);
  });

  it('an inline alignment’s rows by key for `(row :key k)`: its labels, positions of its rows system', () => {
    const b = createSeqSpec();
    const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: ['h', 'm'], letters: ['MK-', 'MKL'] } });
    const m = resolveSpec(b.build());
    expect(m.row?.(`${aln}/(rows)`, 'm')).toBe(1);
    expect(m.row?.(`${aln}/(rows)`, 'x')).toBeUndefined();
    const r = evaluateQuery(b.build(), m, `${aln}/(rect :a ${aln}/(row :key m) :b (range :from 2 :to 3))`, {
      limit: 100,
    });
    expect(r.issues).toEqual([]);
    expect(r.result?.type === 'location' && [...r.result.location.rects[0].a.ranges]).toEqual([1, 2]);
  });

  it('reports a source no runtime reads, and keeps it unavailable', () => {
    const bare = resolveSpec(P53);
    expect(bare.sources.get('/data/chr17')!.runtime).toBeUndefined();
    expect(codes(bare.issues)).toEqual(['seq.resolve.source.runtime']);
    expect(bare.issues[0].address).toBe('/data/chr17');
  });

  it('answers the query forms over the spec', () => {
    const q = (a: string) => evaluateQuery(P53, model, a, { limit: 100 });
    const r = q('/data/p53/(range :from 2 :to 14)');
    expect(r.result?.type === 'location' && ranges(r.result.location)).toEqual([['/data/p53', [1, 14]]]);
    // the design's example query on a protein of 30 residues: nothing of it is left
    const far = q('/data/p53/(range :from 94 :to 292)');
    expect(codes(far.issues)).toContain('seq.query.range.empty');
    expect(far.result?.type === 'location' && far.result.location.parts).toEqual([]);
    // the open gene answers labels anywhere
    const g = q('/data/tp53/(range :from 94 :to 292)');
    expect(g.result?.type === 'location' && ranges(g.result.location)).toEqual([['/data/tp53', [93, 292]]]);
    expect(model.system('/data/tp53')!.length).toBe(OpenLength);
    const items = q('/data/domains/(items :overlapping /data/p53/(range :from 14 :to 20))');
    expect(items.result?.type === 'keys' && items.result.keys.keys.map((k) => k.key)).toEqual([0, 1]);
    const label = q('/data/domains/(item :label TAD)');
    expect(label.result?.type === 'keys' && label.result.keys.keys.map((k) => k.key)).toEqual([1]);
  });

  it('answers no chain and no dynamic node: no structures and no overlay in tier 1', () => {
    const r = evaluateQuery(P53, model, '/data/p53/(chain :id A)', { limit: 10 });
    expect(r.result?.type === 'nodes' && r.result.nodes).toEqual([]);
    expect(codes(r.issues)).toEqual(['seq.query.node.missing']);
    expect(model.annotation('/data/conservation')).toBeUndefined();
  });
});

describe('static model: projection through the codon map', () => {
  const at = (m: StaticModel, l: Location, to: string, close = false) => ranges(m.project(l, to, { close }).location);

  it('minus strand: a protein range onto the gene and back', () => {
    expect(at(model, rangesLocation('/data/p53', [0, 1]), '/data/tp53')).toEqual([['/data/tp53', [197, 200]]]);
    expect(at(model, rangesLocation('/data/p53', [13, 14]), '/data/tp53')).toEqual([
      ['/data/tp53', [148, 150, 160, 161]],
    ]);
    expect(at(model, rangesLocation('/data/tp53', [160, 161]), '/data/p53')).toEqual([['/data/p53', [13, 14]]]);
    expect(at(model, rangesLocation('/data/tp53', [150, 160]), '/data/p53')).toEqual([]);
  });

  it('plus strand: in genomic order', () => {
    const plus = resolveSpec(variant('+'), { runtimes: [Fasta] });
    expect(at(plus, rangesLocation('/data/p53', [0, 1]), '/data/tp53')).toEqual([['/data/tp53', [100, 103]]]);
    expect(at(plus, rangesLocation('/data/p53', [0, 30]), '/data/tp53')).toEqual([
      ['/data/tp53', [100, 150, 160, 200]],
    ]);
    expect(at(plus, rangesLocation('/data/tp53', [101, 102]), '/data/p53')).toEqual([['/data/p53', [0, 1]]]);
  });

  it('closure: a nucleotide gains its residue and lifts to its whole codon (EC-30)', () => {
    const nt = rangesLocation('/data/tp53', [160, 161]);
    expect(ranges(model.close(nt))).toEqual([
      ['/data/p53', [13, 14]],
      ['/data/tp53', [148, 150, 160, 161]],
    ]);
    expect(at(model, nt, '/data/tp53', true)).toEqual([['/data/tp53', [148, 150, 160, 161]]]);
    // closure is one step, among the systems named
    expect(model.close(nt, ['/data/tp53'])).toBe(nt);
  });

  it('a part with no path is left out with the path’s issue', () => {
    const b = createSeqSpec();
    b.data.sequence('a', { alphabet: 'protein', letters: 'MK' });
    b.data.sequence('b', { alphabet: 'protein', letters: 'MK' });
    const m = resolveSpec(b.build());
    const r = m.project(rangesLocation('/data/a', [0, 1]), '/data/b');
    expect(r.location.parts).toEqual([]);
    expect(codes(r.issues)).toEqual(['seq.resolve.path.none']);
    expect(m.path('/data/a', '/data/b')).toBe(m.path('/data/a', '/data/b'));
  });
});

describe('static model: maps and their issues', () => {
  it('an offset map onto a cyclic sequence wraps, split at the origin (EC-11)', () => {
    const b = createSeqSpec();
    const a = b.data.sequence('a', { alphabet: 'dna', letters: 'ACGTAC' });
    const pl = b.data.sequence('plasmid', { alphabet: 'dna', letters: 'ACGTACGTAC', cyclic: true });
    b.data.map('m', { from: a, to: pl, method: 'offset', by: 7 });
    const m = resolveSpec(b.build());
    expect(m.maps.get('/data/m')!.linear.runs()).toEqual([
      { from: 0, to: 7, length: 3, dir: 1 },
      { from: 3, to: 0, length: 3, dir: 1 },
    ]);
    const r = m.project(rangesLocation('/data/a', [1, 5]), '/data/plasmid');
    expect(ranges(r.location)).toEqual([['/data/plasmid', [0, 2, 8, 10]]]);
  });

  it('a map that pairs nothing is reported; a partial one is not (EC-21)', () => {
    const b = createSeqSpec();
    const a = b.data.sequence('a', { alphabet: 'protein', letters: 'MKMK' });
    const c = b.data.sequence('c', { alphabet: 'protein', letters: 'MKMK' });
    b.data.map('partial', { from: a, to: c, method: 'offset', by: 2 });
    b.data.map('none', { from: a, to: c, method: 'offset', by: 9 });
    const m = resolveSpec(b.build());
    expect(codes(m.issues)).toEqual(['seq.resolve.map.empty']);
    expect(m.issues[0].address).toBe('/data/none');
  });

  it('a coding sequence that is not whole codons is reported (EC-17)', () => {
    const b = createSeqSpec();
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MKW' });
    const g = b.data.sequence('g', { alphabet: 'dna', letters: 'ATGAAATGGTAA' });
    b.data.map('short', { from: p, to: g, method: 'codon', exons: { start: [0], end: [8] }, strand: '+' });
    b.data.map('stop', { from: p, to: g, method: 'codon', exons: { start: [0], end: [12] }, strand: '+' });
    const m = resolveSpec(b.build());
    expect(codes(m.issues)).toEqual(['seq.resolve.map.codon_length']);
    expect(m.issues[0].address).toBe('/data/short');
  });

  it('a reference to no node is reported on an unchecked spec', () => {
    const spec: Node = {
      kind: 'spec',
      params: { version: '0.0.0' },
      children: [
        {
          kind: 'data',
          id: 'data',
          children: [{ kind: 'ranges', id: 'r', params: { target: '/data/nope' } }],
        },
      ],
    };
    expect(codes(resolveSpec(spec).issues)).toEqual(['seq.resolve.ref.missing']);
  });

  it('a wrapped row on a cyclic sequence is two rows of one item (EC-9)', () => {
    const b = createSeqSpec();
    const pl = b.data.sequence('plasmid', { alphabet: 'dna', letters: 'ACGTACGTAC', cyclic: true });
    b.data.ranges('f', { target: pl, rows: { start: [8], end: [2], id: ['ori'] } });
    const m = resolveSpec(b.build());
    expect(m.annotations.get('/data/f')!.items).toEqual([
      { key: 'ori', start: 8, end: 10 },
      { key: 'ori', start: 0, end: 2 },
    ]);
    const r = evaluateQuery(b.build(), m, '/data/f/(extent :of (item :id ori))', { limit: 10 });
    expect(r.result?.type === 'location' && ranges(r.result.location)).toEqual([['/data/plasmid', [0, 2, 8, 10]]]);
  });
});

describe('static model: inline ranges’ strands', () => {
  for (const strand of [null, [null, null], [null, '-']] as const)
    it(`a strand of ${JSON.stringify(strand)} reads as no strand where null`, () => {
      const b = createSeqSpec();
      const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTACGTAC' });
      b.data.ranges('f', { target: s, rows: { start: [0, 4], end: [2, 6], strand } as never });
      const m = resolveSpec(b.build());
      const a = m.annotations.get('/data/f')!;
      expect(a.items).toBeDefined();
      expect(a.items!.map((i) => i.key)).toEqual([0, 1]);
      expect(a.items!.map((i) => i.fields?.strand ?? null)).toEqual(Array.isArray(strand) ? strand : [null, null]);
    });
});

describe('static model: caching by hashes', () => {
  it('an equal spec gives the same model; an edit rebuilds only what it touched', () => {
    const resolver = createResolver({ runtimes: [Fasta] });
    const a = resolver.resolve(P53);
    expect(resolver.resolve(p53Spec())).toBe(a);
    // another source file: the source is new, the map over the same sequences is kept
    const b = resolver.resolve(variant('-'));
    expect(b).not.toBe(a);
    expect(b.sources.get('/data/chr17')).not.toBe(a.sources.get('/data/chr17'));
    expect(b.systems.get('/data/p53')).toBe(a.systems.get('/data/p53'));
    expect(b.maps.get('/data/p53-on-gene')).toBe(a.maps.get('/data/p53-on-gene'));
    const c = resolver.resolve(variant('+'));
    expect(c.maps.get('/data/p53-on-gene')).not.toBe(b.maps.get('/data/p53-on-gene'));
    expect(c.systems.get('/data/p53')).toBe(a.systems.get('/data/p53'));
  });
});
