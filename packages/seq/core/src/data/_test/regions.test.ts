/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn, issueOf, type Node } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { asyncSeqSource, delayedRuntime } from '../async-source.js';
import { SourceIssues } from '../issues.js';
import { dataKey } from '../keys.js';
import type { RegionAnswer, Tables } from '../protocol.js';
import { type SeqData, seqToValue } from '../seq-format.js';
import { SeqRuntime } from '../seq-runtime.js';
import { DataStore, type StoreSettings } from '../store.js';
import { binValues, coarsenBins, withSummaries } from '../summaries.js';
import { FakeClock, flush } from './synthetic.js';

// The synthetic asynchronous `seq` source on a fake clock: every state an entry takes (pending, loaded, empty,
// failed, cancelled, evicted), the nodes an answer brings and how long they stay, and the value pyramid.

const S = '/data/genes';
const ranges = (start: number): Tables => ({
  ranges: {
    form: 'list',
    start: Int32Array.of(start),
    end: Int32Array.of(start + 10),
    id: categoryColumn([`r${start}`]),
  },
});
const transcript = (key: string): Node => ({ kind: 'ranges', params: { key } });

/** `chr1` of 4,096 positions in four pairs of 1,024: data, empty, a failure, data with a node. */
function genes(): SeqData {
  const req = (start: number) => ({ system: 'chr1', start, end: start + 1024, level: 1, shapes: ['ranges' as const] });
  const failed: RegionAnswer = { ok: false, kind: 'failed', issue: issueOf(SourceIssues, 'seq.source.failed', 'Down') };
  return {
    coverage: { systems: [{ name: 'chr1', length: 4096 }], shapes: [{ system: 'chr1', shape: 'ranges' }], levels: [1] },
    pairs: [
      { request: req(0), answer: { ok: true, tables: ranges(10) } },
      { request: req(1024), answer: { ok: true, tables: {} } },
      { request: req(2048), answer: failed },
      {
        request: req(3072),
        answer: { ok: true, tables: ranges(3100), nodes: [{ key: 'T1', node: transcript('T1') }] },
      },
    ],
  };
}

function setup(settings: Partial<StoreSettings> = {}) {
  const clock = new FakeClock();
  const store = new DataStore({
    clock,
    settings: { retry: { attempts: 1, baseMs: 10, factor: 2, cooldownMs: 1000 }, ...settings },
  });
  store.attach(S, asyncSeqSource(genes(), { clock, delayMs: 50 }));
  return { clock, store };
}

const want = (start: number, end: number) => ({ source: S, system: 'chr1', start, end, level: 1 });
const status = (store: DataStore, tile: number) => store.get(dataKey(S, 'chr1', 1, tile))?.status;

describe('the synthetic asynchronous source', () => {
  it('answers after its delay: pending first, then loaded, empty or failed by the pair', async () => {
    const { clock, store } = setup();
    store.want([want(0, 4096)]);
    expect([0, 1, 2, 3].map((t) => status(store, t))).toEqual(['pending', 'pending', 'pending', 'pending']);
    await clock.advance(49);
    expect(status(store, 0)).toBe('pending');
    await clock.advance(1);
    expect([0, 1, 2, 3].map((t) => status(store, t))).toEqual(['loaded', 'empty', 'failed', 'loaded']);
    expect(store.settled.value).toBe(true);
    expect(store.issues.value.map((i) => i.code)).toEqual(['seq.source.failed']);
  });

  it('cancels what left the view after the grace: the read is aborted, the entry gone, never failed', async () => {
    const { clock, store } = setup({ graceMs: 20 });
    store.want([want(0, 1024)]);
    store.want([]);
    expect(status(store, 0)).toBe('pending');
    await clock.advance(20);
    expect(status(store, 0)).toBeUndefined();
    expect(store.inFlight).toBe(0);
    await clock.advance(100);
    expect(status(store, 0)).toBeUndefined();
    expect(store.issues.value).toEqual([]);
  });

  it('evicts what is not wanted over the byte bound, and its nodes go with it', async () => {
    const { clock, store } = setup({ maxBytes: 1 });
    store.want([want(3072, 4096)]);
    await clock.advance(50);
    expect(status(store, 3)).toBe('loaded');
    expect(store.getNode(S, 'T1')).toBeDefined();
    store.want([want(0, 10)]);
    await clock.advance(50);
    expect(status(store, 3)).toBeUndefined();
    expect(store.getNode(S, 'T1')).toBeUndefined();
  });

  it('a held node survives its region’s eviction, and goes when let go', async () => {
    const { clock, store } = setup({ maxBytes: 1 });
    const v0 = store.nodeVersion.value;
    store.want([want(3072, 4096)]);
    await clock.advance(50);
    expect(store.nodeVersion.value).toBeGreaterThan(v0);
    expect([...store.nodes.keys()]).toEqual([`${S}#T1`]);
    const release = store.hold(S, 'T1');
    store.want([]);
    await flush();
    expect(status(store, 3)).toBeUndefined();
    expect(store.getNode(S, 'T1')?.node.kind).toBe('ranges');
    const v1 = store.nodeVersion.value;
    release();
    expect(store.getNode(S, 'T1')).toBeUndefined();
    expect(store.nodeVersion.value).toBeGreaterThan(v1);
    // a region answered again brings the node back under the same key
    store.want([want(3072, 4096)]);
    await clock.advance(50);
    expect(store.getNode(S, 'T1')?.key).toBe('T1');
  });

  it('as a runtime: a delayed format opened from a spec’s source entry', async () => {
    const clock = new FakeClock();
    const runtime = delayedRuntime(SeqRuntime, { clock, delayMs: 30, formats: ['seq-async'] });
    expect(runtime.formats).toEqual(['seq-async']);
    const signal = new AbortController().signal;
    const opened = await runtime.open(
      { address: S, format: 'seq-async', value: seqToValue(genes()) },
      { assets: { readText: () => Promise.reject(new Error('none')) }, signal },
    );
    let answered = false;
    const p = opened.read({ system: 'chr1', start: 0, end: 1024, level: 1, shapes: ['ranges'] }, signal);
    void p.then(() => (answered = true));
    await clock.advance(29);
    expect(answered).toBe(false);
    await clock.advance(1);
    expect(answered).toBe(true);
    const aborting = new AbortController();
    const q = opened.read({ system: 'chr1', start: 0, end: 1024, level: 1, shapes: ['ranges'] }, aborting.signal);
    aborting.abort();
    await expect(q).rejects.toMatchObject({ name: 'AbortError' });
    expect(clock.pending).toBe(0);
  });
});

describe('value summaries', () => {
  const values = (n: number) => Float32Array.from({ length: n }, (_, i) => i % 7);

  it('bins per-position values and coarsens bins by count, NaN not counted', () => {
    const v = Float32Array.of(1, 5, Number.NaN, 3, 2, 2, 8, Number.NaN);
    const b = binValues({ form: 'positions', start: 0, value: v }, 4, 8);
    expect(Array.from(b.min)).toEqual([1, 2]);
    expect(Array.from(b.max)).toEqual([5, 8]);
    expect(Array.from(b.count)).toEqual([3, 3]);
    expect(b.mean[0]).toBeCloseTo(3);
    const c = coarsenBins(b, 2);
    expect([c.min[0], c.max[0], c.count[0]]).toEqual([1, 8, 6]);
    expect(c.mean[0]).toBeCloseTo(21 / 6);
    const none = binValues({ form: 'positions', start: 0, value: Float32Array.of(Number.NaN) }, 4, 4);
    expect([none.count[0], Number.isNaN(none.max[0])]).toEqual([0, true]);
  });

  it('a level-1 source gets a pyramid: coarse levels answered from the level below, the source read once', async () => {
    const length = 100_000;
    let reads = 0;
    const inner = asyncSeqSource({
      coverage: { systems: [{ name: 's', length }], shapes: [{ system: 's', shape: 'values' }], levels: [1] },
      pairs: [
        {
          request: { system: 's', start: 0, end: length, level: 1, shapes: ['values'] },
          answer: { ok: true, tables: { values: { form: 'positions', start: 0, value: values(length) } } },
        },
      ],
    });
    const counted = {
      ...inner,
      read: (r: Parameters<typeof inner.read>[0], s: AbortSignal) => (reads++, inner.read(r, s)),
    };
    const src = withSummaries(counted, { factor: 16, maxLength: 1 << 20 });
    expect(src.coverage.levels).toEqual([1, 16, 256, 4096, 65536]);
    const signal = new AbortController().signal;
    const a = await src.read({ system: 's', start: 0, end: 4096 * 1024, level: 4096, shapes: ['values'] }, signal);
    const b = await src.read({ system: 's', start: 16384, end: 32768, level: 16, shapes: ['values'] }, signal);
    expect(reads).toBe(1);
    if (!a.ok || a.tables.values?.form !== 'bins' || !b.ok || b.tables.values?.form !== 'bins') throw new Error('bins');
    expect(a.tables.values.count.length).toBe(Math.ceil(length / 4096));
    expect(a.tables.values.min[0]).toBe(0);
    expect(a.tables.values.max[0]).toBe(6);
    expect(b.tables.values.start).toBe(16384);
    expect(b.tables.values.count[0]).toBe(16);
    // a source with its own levels, or without values, is left as it is
    const own = { ...inner, coverage: { ...inner.coverage, levels: [1, 10] } };
    expect(withSummaries(own)).toBe(own);
    const noValues = { ...inner, coverage: { ...inner.coverage, shapes: [] } };
    expect(withSummaries(noValues)).toBe(noValues);
  });

  it('the store gives a level-1 source its pyramid on attach (a setting; false: none)', () => {
    const data: SeqData = {
      coverage: {
        systems: [{ name: 's', length: 1_000_000 }],
        shapes: [{ system: 's', shape: 'values' }],
        levels: [1],
      },
      pairs: [],
    };
    const on = new DataStore({ clock: new FakeClock() });
    on.attach(S, asyncSeqSource(data));
    expect(on.coverage(S)?.levels).toEqual([1, 16, 256, 4096, 65536]);
    const off = new DataStore({ clock: new FakeClock(), settings: { summaries: false } });
    off.attach(S, asyncSeqSource(data));
    expect(off.coverage(S)?.levels).toEqual([1]);
  });
});
