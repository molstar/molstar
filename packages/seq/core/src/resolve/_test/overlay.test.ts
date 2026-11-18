/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { type NodeEntry, nodeKey, SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { resolveSpec } from '../model.js';
import { buildOverlay, itemAddress, nodeAddress, OverlayBuilder } from '../overlay.js';
import { projectRanges } from '../paths.js';

// The overlay over a static model: dynamic nodes checked by the language, refs to spec or held nodes, a transcript's
// item system and codon map made on demand, paths through them, and nested sources listed on need within a depth.

const S = '/data/genes';

function model() {
  const b = createSeqSpec();
  const genes = b.data.source('genes', { format: 'seq', uri: 'genes.seq' });
  b.data.sequence('gene', { alphabet: 'dna', letters: 'ACGT'.repeat(375) });
  b.data.sequence('other', { alphabet: 'dna', source: genes });
  return resolveSpec(b.build(), { runtimes: [SeqRuntime] });
}

/** A transcript on the gene: minus strand, two exons of 300 and 462 (254 codons). */
const transcript = (target = '/data/gene'): Node => ({
  kind: 'ranges',
  params: {
    target,
    map: 'codon',
    rows: { start: [100, 960], end: [400, 1422], id: ['T1', 'T1'], strand: ['-', '-'] },
  },
});

const source = (format = 'seq'): Node => ({ kind: 'source', params: { format, value: { coverage: {}, pairs: [] } } });

function reader(nodes: [source: string, key: string, node: Node][], version = 1) {
  const map = new Map<string, NodeEntry>();
  for (const [s, key, node] of nodes) map.set(nodeKey(s, key), { source: s, key, node, version });
  return { nodes: map, nodeVersion: { value: version } };
}

describe('the overlay', () => {
  it('addresses a node by its source and key, an item by its node and id', () => {
    expect(nodeAddress(S, 'T1')).toBe('/data/genes/(node :key T1)');
    expect(itemAddress(nodeAddress(S, 'T1'), 'T1')).toBe('/data/genes/(node :key T1)/(item :id T1)');
  });

  it('checks each node by the language: an invalid one, or one naming nothing held, is left out with an issue', () => {
    const m = model();
    const o = buildOverlay(
      m,
      reader([
        [S, 'T1', transcript()],
        [S, 'bad', { kind: 'ranges', params: { target: 'not an address' } }],
        [S, 'view', { kind: 'view', params: {} }],
        [S, 'lost', transcript('/data/none')],
        [S, 'nested', source()],
        ['/data/unknown', 'X', transcript()],
      ]),
    );
    expect(o.version).toBe(1);
    expect(o.nodes.map((n) => [n.key, n.kind, n.depth])).toEqual([
      ['T1', 'ranges', 0],
      ['nested', 'source', 0],
    ]);
    expect(o.issues.map((i) => [i.code, i.address])).toEqual([
      ['seq.resolve.node.invalid', '/data/genes/(node :key bad)'],
      ['seq.resolve.node.invalid', '/data/genes/(node :key view)'],
      ['seq.resolve.node.ref', '/data/genes/(node :key lost)'],
    ]);
    // the static model never holds a dynamic node
    expect(m.annotations.has(nodeAddress(S, 'T1'))).toBe(false);
    expect(o.annotation(nodeAddress(S, 'T1'))?.items?.length).toBe(2);
  });

  it('gives a coding transcript an item system and a codon map onto its target, on demand', () => {
    const o = buildOverlay(model(), reader([[S, 'T1', transcript()]]));
    const item = itemAddress(nodeAddress(S, 'T1'), 'T1');
    const sys = o.system(item);
    expect(sys).toMatchObject({ address: item, alphabet: 'protein', length: 254, name: 'T1' });
    expect(o.system(itemAddress(nodeAddress(S, 'T1'), 'T9'))).toBeUndefined();
    const down = o.path(item, '/data/gene').path!;
    expect(down.step).toBe(1);
    // residue 0 is the last three nucleotides of the top exon (minus strand); residues 153 and 154 are the first
    // codon of the top exon and the last of the lower one, in the protein's order
    expect(projectRanges(down, [0, 1])).toEqual([1419, 1422]);
    expect(projectRanges(down, [153, 155])).toEqual([960, 963, 397, 400]);
    const up = o.path('/data/gene', item).path!;
    expect(projectRanges(up, [960, 963])).toEqual([153, 154]);
    // between spec systems, the static model's own path
    const m = model();
    const o2 = buildOverlay(m, reader([[S, 'T1', transcript()]]));
    expect(o2.path('/data/gene', '/data/other')).toBe(m.path('/data/gene', '/data/other'));
  });

  it('checks each node by the language’s own rules too: a range ending before it starts is refused, no item made', () => {
    const backwards: Node = {
      kind: 'ranges',
      params: { target: '/data/gene', map: 'codon', rows: { start: [400, 960], end: [100, 1422], id: ['T1', 'T1'] } },
    };
    const o = buildOverlay(model(), reader([[S, 'T1', backwards]]));
    expect(o.nodes).toEqual([]);
    expect(o.issues.map((i) => [i.code, i.address, i.path])).toEqual([
      ['seq.lang.ranges.interval', nodeAddress(S, 'T1'), ['params', 'rows', 'start', 0]],
    ]);
    expect(o.system(itemAddress(nodeAddress(S, 'T1'), 'T1'))).toBeUndefined();
    expect(o.path(itemAddress(nodeAddress(S, 'T1'), 'T1'), '/data/gene').path).toBeUndefined();
    // past the end of its target (the bound read from the target's length)
    const past: Node = { kind: 'ranges', params: { target: '/data/gene', rows: { start: [10], end: [1600] } } };
    expect(buildOverlay(model(), reader([[S, 'X', past]])).issues.map((i) => i.code)).toEqual([
      'seq.lang.ranges.bounds',
    ]);
  });

  it('an item on both strands takes its first row’s, with an issue; its letters are its source’s', () => {
    const mixed = transcript();
    (mixed.params!.rows as Record<string, unknown>).strand = ['-', '+'];
    const o = buildOverlay(model(), reader([[S, 'T1', mixed]]));
    const sys = o.system(itemAddress(nodeAddress(S, 'T1'), 'T1'))!;
    expect(sys.access).toEqual({ kind: 'source', source: S, system: 'T1', shape: 'letters' });
    expect(o.issues.map((i) => i.code)).toEqual(['seq.resolve.item.strand']);
    expect(o.path(sys.address, '/data/gene').path?.step).toBe(1);
  });

  it('resolves a nested node on a shallower node’s item', () => {
    const item = itemAddress(nodeAddress(S, 'T1'), 'T1');
    const nested = nodeAddress(S, 'T1-domains');
    const domains: Node = { kind: 'ranges', params: { target: item, source: nested } };
    const o = buildOverlay(
      model(),
      reader([
        [S, 'T1', transcript()],
        [S, 'T1-domains', source()],
        [nested, 'domains', domains],
      ]),
    );
    const n = o.node(nodeAddress(nested, 'domains'))!;
    expect([n.depth, n.root]).toEqual([1, S]);
    expect(o.annotation(n.address)).toMatchObject({
      target: item,
      access: { kind: 'source', source: nested, system: 'T1' },
    });
    expect(o.issues).toEqual([]);
  });

  it('lists a nested source only when a template needs its depth, never past the limit (EC-179)', () => {
    const a = nodeAddress(S, 'a');
    const b = nodeAddress(a, 'b');
    const nodes = reader([
      [S, 'a', source()],
      [a, 'b', source()],
      [b, 'c', source('none')],
    ]);
    const o = buildOverlay(model(), nodes, { runtimes: [SeqRuntime], config: { maxNestedDepth: 2 } });
    expect(o.nested(new Map())).toEqual([]);
    expect(o.nested(new Map([[S, 1]])).map((n) => [n.address, n.depth, n.runtime?.name])).toEqual([[a, 1, 'seq']]);
    expect(o.nested(new Map([[S, 5]])).map((n) => n.address)).toEqual([a, b]);
    expect(o.issues.map((i) => [i.code, i.address])).toEqual([['seq.resolve.node.depth', nodeAddress(b, 'c')]]);
  });

  it('lists four levels of nested sources by default, and refuses the fifth (EC-179)', () => {
    const addresses = [nodeAddress(S, 'l1')];
    for (let i = 2; i <= 5; i++) addresses.push(nodeAddress(addresses[i - 2], `l${i}`));
    const nodes = reader([
      [S, 'l1', source()],
      ...addresses.slice(0, 4).map((a, i): [string, string, Node] => [a, `l${i + 2}`, source()]),
    ]);
    const o = buildOverlay(model(), nodes, { runtimes: [SeqRuntime] });
    expect(o.nested(new Map([[S, 99]])).map((n) => n.depth)).toEqual([1, 2, 3, 4]);
    expect(o.issues.map((i) => [i.code, i.address])).toEqual([['seq.resolve.node.depth', addresses[4]]]);
  });

  it('is one overlay a node version: the builder keeps it while the version and the model hold', () => {
    const m = model();
    const builder = new OverlayBuilder();
    const r = reader([[S, 'T1', transcript()]], 3);
    const o = builder.of(m, r);
    expect(builder.of(m, r)).toBe(o);
    const r2 = { ...r, nodeVersion: { value: 4 } };
    expect(builder.of(m, r2)).not.toBe(o);
    expect(builder.of(m, reader([], 5))).toMatchObject({ version: 5, nodes: [] });
  });
});
