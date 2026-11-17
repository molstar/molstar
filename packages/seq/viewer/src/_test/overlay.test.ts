/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { type SeqData, SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { nodeAddress, resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import { StoreModel } from '../models/store.js';

// The viewer holds tier 2: the store model gives a frame the overlay of the store's dynamic nodes (one a node
// version), opens a dynamic annotation's rows and the nested sources the view's templates reach, and lets them go
// with their node.

const assets = { readText: () => Promise.reject(new Error('no network in tests')) };
const Genes = '/data/genes';

/** A gene region whose one answer brings a transcript and a nested source (whose coverage lists nothing). */
function genes(): SeqData {
  const transcript: Node = {
    kind: 'ranges',
    params: { target: '/data/gene', rows: { start: [10], end: [40], id: ['T1'] } },
  };
  const nested: Node = {
    kind: 'source',
    params: { format: 'seq', value: seqToValue({ coverage: { systems: [], shapes: [], levels: [1] }, pairs: [] }) },
  };
  return {
    coverage: { systems: [{ name: 'g', length: 100 }], shapes: [{ system: 'g', shape: 'letters' }], levels: [1] },
    pairs: [
      {
        request: { system: 'g', start: 0, end: 100, level: 1, shapes: ['letters'] },
        answer: {
          ok: true,
          tables: { letters: { form: 'letters', start: 0, letters: 'A'.repeat(100) } },
          nodes: [
            { key: 'T1', node: transcript },
            { key: 'T1-more', node: nested },
          ],
        },
      },
    ],
  };
}

function spec(depth: number) {
  const b = createSeqSpec();
  const src = b.data.source('genes', { format: 'seq', value: seqToValue(genes()) });
  const gene = b.data.sequence('gene', { alphabet: 'dna', source: src, name: 'g' });
  b.views.view('main', (v) =>
    v
      .section('gene', { of: gene })
      .track('t', { for_each: { source: src, node_kind: 'ranges', depth } })
      .layer({ type: 'ranges' }),
  );
  const model = resolveSpec(b.build(), { runtimes: [SeqRuntime] });
  return { model, view: resolveViews(model)[0] };
}

const ctx = { runtimes: [SeqRuntime], assets };
/** Lets the openings' promises settle (no timer). */
const flush = () => new Promise((r) => setImmediate(r));

describe('StoreModel: the overlay', () => {
  it('is one a node version; it opens the nodes’ rows, and a nested source only for a template that reaches it', async () => {
    const m = new StoreModel();
    const { model, view } = spec(0);
    await m.sync(model, assets);
    const empty = m.overlay(model, view, ctx);
    expect(empty.nodes).toEqual([]);
    const region = [{ source: Genes, system: 'g', start: 0, end: 100, level: 1 }];
    m.store.want(region);
    await m.store.whenSettled();
    const o = m.overlay(model, view, ctx);
    expect(o).not.toBe(empty);
    expect(m.overlay(model, view, ctx)).toBe(o);
    expect(o.nodes.map((n) => n.key)).toEqual(['T1', 'T1-more']);
    await flush();
    expect(m.store.coverage(nodeAddress(Genes, 'T1'))?.systems[0].length).toBe(40);
    expect(m.store.coverage(nodeAddress(Genes, 'T1-more'))).toBeUndefined();

    const deeper = spec(1);
    // another spec: its source opened again, its region wanted again (a frame's wants)
    await m.sync(deeper.model, assets);
    m.store.want(region);
    await m.store.whenSettled();
    m.overlay(deeper.model, deeper.view, ctx);
    await flush();
    expect(m.store.coverage(nodeAddress(Genes, 'T1-more'))).toBeDefined();
    expect(m.issues.value).toEqual([]);

    // the nodes go with their region: their sources are let go at the next frame
    m.store.want([]);
    m.store.detach(Genes);
    m.overlay(deeper.model, deeper.view, ctx);
    expect(m.store.coverage(nodeAddress(Genes, 'T1'))).toBeUndefined();
    expect(m.store.coverage(nodeAddress(Genes, 'T1-more'))).toBeUndefined();
    m.retire();
  });
});
