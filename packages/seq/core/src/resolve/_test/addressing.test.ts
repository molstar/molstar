/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { applyPatch } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import { type NodeEntry, nodeKey, SeqRuntime } from '@molstar/seqx-core/data';
import { addRefPatch, createSeqSpec, SeqLanguage } from '@molstar/seqx-core/lang';
import { p53Spec } from '../../lang/_test/p53.js';
import { addressFor, locationOfResult, rangesLocation, resolveAddressIn } from '@molstar/seqx-core/query';
import { resolveView } from '@molstar/seqx-core/view';
import { resolveSpec } from '../model.js';
import { buildOverlay, EmptyOverlay } from '../overlay.js';
import { overlayQueryModel } from '../query-model.js';

// Addressing over a real spec: the static model and an overlay as the query model, so a ref to a dynamic node waits
// while its node is not held and resolves once it is; a layer drawing through a ref; and the performance
// representative, 1,000 addresses of ids, forms and refs on the p53 spec.

const S = '/data/genes';

function spec(): Node {
  const b = createSeqSpec();
  b.data.source('genes', { format: 'seq', uri: 'genes.seq' });
  b.data.sequence('gene', { alphabet: 'dna', letters: 'ACGT'.repeat(375) });
  b.refs.ref('t1', { to: `${S}/(node :key T1)`, label: 'T1' });
  return b.build();
}

/** A transcript on the gene: minus strand, two exons. */
const transcript: Node = {
  kind: 'ranges',
  params: {
    target: '/data/gene',
    map: 'codon',
    rows: { start: [100, 960], end: [400, 1422], id: ['T1', 'T1'], strand: ['-', '-'] },
  },
};

function reader(nodes: [source: string, key: string, node: Node][], version = 1) {
  const map = new Map<string, NodeEntry>();
  for (const [s, key, node] of nodes) map.set(nodeKey(s, key), { source: s, key, node, version });
  return { nodes: map, nodeVersion: { value: version } };
}

describe('a ref to a dynamic node', () => {
  it('waits while its node is not held, and resolves when it arrives', () => {
    const root = spec();
    const model = resolveSpec(root, { runtimes: [SeqRuntime] });
    const before = resolveAddressIn(root, overlayQueryModel(model, EmptyOverlay), '/refs/t1');
    expect(before.status).toBe('waiting');
    expect(before.issues.map((i) => i.code)).toEqual(['seq.address.waiting']);

    const overlay = buildOverlay(model, reader([[S, 'T1', transcript]]));
    const q = overlayQueryModel(model, overlay);
    expect(overlayQueryModel(model, overlay)).toBe(q);
    const after = resolveAddressIn(root, q, '/refs/t1');
    expect(after.status).toBe('resolved');
    expect(locationOfResult(root, q, after.result!).location).toEqual(
      rangesLocation('/data/gene', [100, 400, 960, 1422]),
    );
    // a form after the ref reads the node's items (the item's system is not a prefix of location forms yet)
    const item = resolveAddressIn(root, q, '/refs/t1/(item :id T1)');
    expect(item.result).toMatchObject({
      type: 'elements',
      elements: { keys: [{ owner: `${S}/(node :key T1)`, key: 'T1' }] },
    });
    // and the way back names the node by the ref
    expect(addressFor(root, q, after.result!)).toBe('/refs/t1');
  });

  it('is added as one patch, whether or not its node is held', () => {
    const root = spec();
    const r = applyPatch(root, addRefPatch(root, { id: 't2', to: `${S}/(node :key T2)` }), { language: SeqLanguage });
    expect(r.ok).toBe(true);
  });
});

describe('a layer drawing through a ref', () => {
  it('is placed with the ref’s target as its data', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTACGT' });
    const r = b.data.ranges('r', { target: s, rows: { start: [0], end: [4] } });
    b.refs.ref('feature', { to: r });
    b.views.view('v', (v) => {
      v.section('sec', { of: s }).track('t').layer('l', { type: 'ranges', data: '/refs/feature' });
    });
    const view = resolveView(resolveSpec(b.build()), '/views/v')!;
    const layer = view.sections[0].tracks[0].layers[0];
    expect(layer.placed).toBe(true);
    expect(layer.data).toMatchObject({ kind: 'address', address: '/data/r', system: '/data/s' });
    expect(view.issues).toEqual([]);
  });
});

describe('resolving many addresses', () => {
  const root = (() => {
    const s = p53Spec();
    let r = s;
    const refs: [string, string][] = [
      ['protein', '/data/p53'],
      ['tad', '/data/domains/(item :label TAD)'],
      ['start', '/refs/protein/(range :from 1 :to 10)'],
    ];
    for (const [id, to] of refs) {
      const p = applyPatch(r, addRefPatch(r, { id, to }), { language: SeqLanguage });
      if (!p.ok) throw new Error('the refs did not apply');
      r = p.root;
    }
    return r;
  })();
  const model = resolveSpec(root);
  const addresses = Array.from({ length: 1000 }, (_, i) => {
    switch (i % 5) {
      case 0:
        return '/data/p53';
      case 1:
        return `/data/p53/(range :from ${1 + (i % 20)} :to ${10 + (i % 20)})`;
      case 2:
        return '/refs/start';
      case 3:
        return `/data/domains/(item :label TAD)`;
      default:
        return `/refs/protein/(at :label ${1 + (i % 30)} :numbering ordinal)`;
    }
  });
  const all = () => addresses.map((a) => resolveAddressIn(root, model, a));

  it('resolves ids, forms and refs mixed', () => {
    const r = all();
    expect(r.every((x) => x.status === 'resolved')).toBe(true);
    expect(r[2].via).toEqual(['/refs/start', '/refs/protein']);
  });

  itWithin('1,000 addresses (ids, forms and refs mixed) on the p53 spec', (within) => {
    all(); // warm
    // measured: 6 ms median, 5.6 to 12.9 over 20 runs (an M-series laptop); one sample, so 40 ms leaves room for a
    // collection and a parallel heavy run
    const r = within(40, all);
    expect(r.length).toBe(1000);
  });
});
