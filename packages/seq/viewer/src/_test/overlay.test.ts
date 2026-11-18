/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf, type Node } from '@molstar/core/repr/tree';
import {
  type FormatRuntime,
  type OpenedSource,
  type SeqData,
  SeqRuntime,
  SourceOpenError,
  seqToValue,
} from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { nodeAddress, ResolveIssues, resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it, vi } from 'vitest';
import { StoreModel } from '../models/store.js';

// The viewer holds tier 2: the store model gives a frame the overlay of the store's dynamic nodes (one a node
// version), opens a dynamic annotation's rows and the nested sources the view's templates reach, and lets them go
// with their node.

const assets = { readText: () => Promise.reject(new Error('no network in tests')) };
const Genes = '/data/genes';

/** A gene region whose one answer brings a transcript and a nested source (whose coverage lists nothing). */
function genes(nestedFormat = 'seq'): SeqData {
  const transcript: Node = {
    kind: 'ranges',
    params: { target: '/data/gene', rows: { start: [10], end: [40], id: ['T1'] } },
  };
  const nested: Node = {
    kind: 'source',
    params: {
      format: nestedFormat,
      value: seqToValue({ coverage: { systems: [], shapes: [], levels: [1] }, pairs: [] }),
    },
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

function spec(depth: number, nestedFormat = 'seq', config: { maxNestedDepth?: number } = {}) {
  const b = createSeqSpec();
  const src = b.data.source('genes', { format: 'seq', value: seqToValue(genes(nestedFormat)) });
  const gene = b.data.sequence('gene', { alphabet: 'dna', source: src, name: 'g' });
  b.views.view('main', (v) =>
    v
      .section('gene', { of: gene })
      .track('t', { for_each: { source: src, node_kind: 'ranges', depth } })
      .layer('ranges', { type: 'ranges' }),
  );
  const model = resolveSpec(b.build(), { runtimes: [SeqRuntime], config });
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

  it('retired with a nested source opening: the late opening is disposed, nothing attached, no timer left', async () => {
    vi.useFakeTimers();
    let finish: (o: OpenedSource) => void = () => {};
    let disposed = false;
    const slow: FormatRuntime = {
      name: 'slow',
      formats: ['slow'],
      open: () => new Promise<OpenedSource>((r) => (finish = r)),
    };
    const m = new StoreModel();
    const { model, view } = spec(1, 'slow');
    await m.sync(model, assets);
    m.store.want([{ source: Genes, system: 'g', start: 0, end: 100, level: 1 }]);
    await m.store.whenSettled();
    m.overlay(model, view, { runtimes: [SeqRuntime, slow], assets });
    const store = m.store;
    m.retire();
    finish({
      coverage: { systems: [], shapes: [], levels: [1] },
      read: () => Promise.reject(new Error('no')),
      dispose: () => {
        disposed = true;
      },
    });
    await Promise.resolve();
    await Promise.resolve();
    expect(disposed).toBe(true);
    expect(store.coverage(nodeAddress(Genes, 'T1-more'))).toBeUndefined();
    expect(store.coverage(Genes)).toBeUndefined();
    expect(vi.getTimerCount()).toBe(0);
    vi.useRealTimers();
  });

  it('a nested source that cannot be opened is an issue while its node is held, cleared when it goes', async () => {
    const broken: FormatRuntime = {
      name: 'broken',
      formats: ['broken'],
      open: () =>
        Promise.reject(new SourceOpenError([issueOf(ResolveIssues, 'seq.resolve.source.runtime', 'Cannot open')])),
    };
    const m = new StoreModel();
    const { model, view } = spec(1, 'broken');
    await m.sync(model, assets);
    m.store.want([{ source: Genes, system: 'g', start: 0, end: 100, level: 1 }]);
    await m.store.whenSettled();
    const c = { runtimes: [SeqRuntime, broken], assets };
    m.overlay(model, view, c);
    await flush();
    await flush();
    expect(m.issues.value.map((i) => i.address)).toEqual([nodeAddress(Genes, 'T1-more')]);
    m.store.detach(Genes);
    m.overlay(model, view, c);
    await flush();
    expect(m.issues.value).toEqual([]);
    m.retire();
  });

  it('reads the model’s configuration and the context’s runtimes: either changed, another overlay', async () => {
    const m = new StoreModel();
    const { model, view } = spec(1, 'other');
    await m.sync(model, assets);
    m.store.want([{ source: Genes, system: 'g', start: 0, end: 100, level: 1 }]);
    await m.store.whenSettled();
    const needs = new Map([[Genes, 1]]);
    const a = m.overlay(model, view, ctx);
    expect(a.nested(needs).map((n) => n.runtime?.name)).toEqual([undefined]);
    const other: FormatRuntime = { ...SeqRuntime, name: 'other', formats: ['other'] };
    const b = m.overlay(model, view, { runtimes: [SeqRuntime, other], assets });
    expect(b).not.toBe(a);
    expect(b.nested(needs).map((n) => n.runtime?.name)).toEqual(['other']);
    // a model resolved with another nesting limit: its overlay reads it
    const shallow = spec(1, 'other', { maxNestedDepth: 0 });
    const c = m.overlay(shallow.model, shallow.view, { runtimes: [SeqRuntime, other], assets });
    expect(c.nested(needs)).toEqual([]);
    expect(c.issues.map((i) => i.code)).toContain('seq.resolve.node.depth');
    m.retire();
  });
});
