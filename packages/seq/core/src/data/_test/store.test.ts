/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive } from '@molstar/core/reactive';
import { issueOf } from '@molstar/core/repr/tree';
import { describe, expect, it } from 'vitest';
import { SourceIssues } from '../issues.js';
import { dataKey } from '../keys.js';
import type { RegionAnswer, RegionRequest } from '../protocol.js';
import { DataStore, type StoreSettings, type Want } from '../store.js';
import { FakeClock, flush, rangesAnswer, SyntheticSource, syntheticCoverage } from './synthetic.js';

const S = '/data/genes';

function setup(
  options: {
    settings?: Partial<StoreSettings>;
    delayMs?: number;
    answer?: (r: RegionRequest, n: number) => RegionAnswer | Error;
  } = {},
) {
  const clock = new FakeClock();
  const source = new SyntheticSource(syntheticCoverage(), clock, options);
  const store = new DataStore({ clock, settings: options.settings });
  store.attach(S, source);
  return { clock, source, store };
}

const view = (start: number, end: number, extra: Partial<Want> = {}): Want => ({
  source: S,
  system: 'chr1',
  start,
  end,
  level: 1,
  ...extra,
});
const key = (tile: number, level = 1) => dataKey(S, 'chr1', level, tile);
const failed = (): RegionAnswer => ({
  ok: false,
  kind: 'failed',
  issue: issueOf(SourceIssues, 'seq.source.failed', 'Down'),
});

describe('entries (EC-169)', () => {
  it('are pending, then loaded or empty, with versions that move', async () => {
    const { clock, store } = setup({ answer: (r) => (r.start === 0 ? rangesAnswer(r) : { ok: true, tables: {} }) });
    const v0 = store.version.value;
    store.want([view(0, 2048)]);
    const pending = store.get(key(0))!;
    expect(pending.status).toBe('pending');
    expect(store.version.value).toBeGreaterThan(v0);
    expect(store.settled.value).toBe(false);
    await clock.advance(100);
    const loaded = store.get(key(0))!;
    expect(loaded.status).toBe('loaded');
    expect(loaded.version).toBeGreaterThan(pending.version);
    expect(loaded.region).toMatchObject({ system: 'chr1', start: 0, end: 1024, tile: 0 });
    expect(store.get(key(1))!.status).toBe('empty');
    expect(store.settled.value).toBe(true);
  });

  it('tell only what reads them', async () => {
    const { clock, store } = setup({ delayMs: 10 });
    let a = 0;
    let b = 0;
    store.entry(key(0)).listen(() => a++);
    derive(store.entry(key(5)), (e) => e?.status).listen(() => b++);
    store.want([view(0, 1024)]);
    await clock.advance(10);
    expect(a).toBe(2);
    expect(b).toBe(0);
  });

  it('are one per tile: overlapping wants and aliases share a request', async () => {
    const { clock, source, store } = setup();
    store.want([view(0, 1500), view(1000, 2000, { system: '1', priority: 'margin' })]);
    expect(source.reads.map((r) => r.request.start)).toEqual([0, 1024]);
    expect(store.keysOf(view(1000, 2000, { system: '1' }))).toEqual([key(0), key(1)]);
    await clock.advance(100);
    store.want([view(0, 2048)]);
    expect(source.reads.length).toBe(2);
  });

  it('snap an unlisted level to the coarsest listed no wider (EC-170)', () => {
    const { source, store } = setup();
    store.want([view(0, 50_000, { level: 37 })]);
    expect(source.reads.map((r) => r.request.level)).toEqual([10, 10, 10, 10, 10]);
  });
});

describe('requests (EC-172)', () => {
  it('go view, margin, prefetch, at most the bound in flight', async () => {
    const { clock, source, store } = setup({ settings: { maxInFlight: 3 } });
    store.want([
      view(10 * 1024, 11 * 1024, { priority: 'prefetch' }),
      view(5 * 1024, 7 * 1024, { priority: 'margin' }),
      view(0, 2 * 1024),
    ]);
    expect(source.reads.map((r) => r.request.start / 1024)).toEqual([0, 1, 5]);
    expect(store.inFlight).toBe(3);
    await clock.advance(100);
    expect(source.reads.map((r) => r.request.start / 1024)).toEqual([0, 1, 5, 6, 10]);
    await clock.advance(100);
    expect(store.settled.value).toBe(true);
  });
});

describe('cancellation (EC-173)', () => {
  it('waits a grace: a jitter back keeps the request, a move away aborts it', async () => {
    const { clock, source, store } = setup({ delayMs: 1000, settings: { graceMs: 300 } });
    store.want([view(0, 1024)]);
    store.want([view(4096, 5120)]);
    await clock.advance(200);
    store.want([view(0, 1024)]);
    await clock.advance(200);
    expect(source.reads[0].aborted).toBe(false);
    store.want([view(4096, 5120)]);
    await clock.advance(300);
    expect(source.reads[0].aborted).toBe(true);
    expect(store.get(key(0))).toBeUndefined();
    expect(store.issues.value).toEqual([]);
    await clock.advance(1000);
    expect(store.get(key(4))?.status).toBe('loaded');
  });

  it('drops what was not yet asked for at once, and the oldest stale goes first', async () => {
    const { clock, source, store } = setup({ delayMs: 1000, settings: { maxInFlight: 2, graceMs: 300 } });
    store.want([view(0, 3 * 1024)]);
    expect(source.reads.length).toBe(2);
    expect(store.get(key(2))?.status).toBe('pending');
    await clock.advance(100);
    store.want([view(0, 1024)]);
    expect(store.get(key(2))).toBeUndefined();
    await clock.advance(100);
    store.want([]);
    await clock.advance(200);
    expect(source.reads.map((r) => r.aborted)).toEqual([false, true]);
    await clock.advance(100);
    expect(source.reads.map((r) => r.aborted)).toEqual([true, true]);
    expect(source.reads.length).toBe(2);
  });

  it('drops a late answer to an aborted read', async () => {
    const clock = new FakeClock();
    const late = new SyntheticSource(syntheticCoverage(), clock, { delayMs: 500 });
    // A source that ignores aborts: it answers anyway.
    const deaf = {
      coverage: late.coverage,
      dispose() {},
      read: (r: RegionRequest) => late.read(r, new AbortController().signal),
    };
    const store = new DataStore({ clock, settings: { graceMs: 100 } });
    store.attach(S, deaf);
    store.want([view(0, 1024)]);
    store.want([]);
    await clock.advance(1000);
    expect(late.reads[0].answered).toBe(true);
    expect(store.get(key(0))).toBeUndefined();
    expect(store.inFlight).toBe(0);
  });
});

describe('failures (EC-174)', () => {
  it('are retried with backoff while pending, then fail, cool down and are tried again; reported once', async () => {
    const { clock, source, store } = setup({
      delayMs: 10,
      settings: { retry: { attempts: 3, baseMs: 500, factor: 2, cooldownMs: 30_000 } },
      answer: (r, n) => (n === 5 ? rangesAnswer(r) : failed()),
    });
    store.want([view(0, 1024)]);
    await clock.advance(10);
    expect(store.get(key(0))?.status).toBe('pending');
    expect(store.settled.value).toBe(false);
    await clock.advance(500);
    expect(source.reads.length).toBe(2);
    await clock.advance(1010);
    expect(source.reads.length).toBe(3);
    await clock.advance(10);
    expect(store.get(key(0))).toMatchObject({ status: 'failed', kind: 'failed' });
    expect(store.settled.value).toBe(true);
    await clock.advance(30_000);
    expect(source.reads.length).toBe(4);
    expect(store.get(key(0))?.status).toBe('failed');
    await clock.advance(10);
    expect(store.get(key(0))?.status).toBe('failed');
    await clock.advance(510);
    expect(source.reads.length).toBe(5);
    expect(store.get(key(0))?.status).toBe('loaded');
    expect(store.issues.value.map((i) => [i.code, i.address])).toEqual([['seq.source.failed', S]]);
  });

  it('a read that rejects is a failure; not served is never retried', async () => {
    const { clock, source, store } = setup({
      delayMs: 10,
      answer: (r) =>
        r.start === 0
          ? new Error('socket')
          : { ok: false, kind: 'not-served', issue: issueOf(SourceIssues, 'seq.source.not_served', 'No') },
    });
    store.want([view(0, 2048)]);
    await clock.advance(10);
    expect(store.get(key(1))).toMatchObject({ status: 'failed', kind: 'not-served' });
    expect(store.get(key(0))?.status).toBe('pending');
    await clock.advance(60_000);
    expect(source.reads.filter((r) => r.request.start === 1024).length).toBe(1);
    expect(source.reads.filter((r) => r.request.start === 0).length).toBeGreaterThan(3);
    expect(store.issues.value.map((i) => i.code)).toEqual(['seq.source.failed', 'seq.source.not_served']);
  });

  it('a malformed answer is a failure with its issue', async () => {
    const malformed: RegionAnswer = {
      ok: true,
      tables: { ranges: { form: 'rows', start: Int32Array.of(1, 2), end: Int32Array.of(3) } },
    };
    const { clock, store } = setup({
      delayMs: 10,
      settings: { retry: { attempts: 1, baseMs: 1, factor: 1, cooldownMs: 1e6 } },
      answer: () => malformed,
    });
    store.want([view(0, 1024)]);
    await clock.advance(10);
    expect(store.get(key(0))).toMatchObject({
      status: 'failed',
      kind: 'malformed',
      issue: { code: 'seq.source.malformed' },
    });
  });
});

describe('eviction', () => {
  // Each tile of the synthetic answer: 11 ranges of 4 + 4 bytes, ids 'r…' in a category.
  const tileBytes = (store: DataStore) => (store.get(key(0)) as { bytes: number }).bytes;

  it('drops the least recently used outside the view, never what is in view', async () => {
    const probe = setup({ delayMs: 1 });
    probe.store.want([view(0, 1024)]);
    await probe.clock.advance(1);
    const b = tileBytes(probe.store);
    const { clock, store } = setup({ delayMs: 1, settings: { maxBytes: 3 * b + b / 2 } });
    store.want([view(0, 1024)]);
    await clock.advance(1);
    store.want([view(1024, 2048)]);
    await clock.advance(1);
    store.want([view(2048, 3072), view(3072, 4096, { priority: 'prefetch' })]);
    await clock.advance(1);
    expect(store.get(key(0))).toBeUndefined();
    expect(store.get(key(1))?.status).toBe('loaded');
    store.want([view(4096, 5120), view(3072, 4096, { priority: 'prefetch' })]);
    await clock.advance(1);
    expect(store.get(key(1))).toBeUndefined();
    expect([2, 3, 4].map((t) => store.get(key(t))?.status)).toEqual(['loaded', 'loaded', 'loaded']);
    expect(store.bytes).toBeLessThanOrEqual(store.settings.maxBytes);
    store.want([view(0, 5 * 1024)]);
    await clock.advance(1);
    expect(store.keysOf(view(0, 5 * 1024)).every((k) => store.get(k)?.status === 'loaded')).toBe(true);
    expect(store.overBudget.value).toBe(true);
    expect(store.issues.value.map((i) => i.code)).toEqual(['seq.store.over_budget']);
    store.want([view(0, 1024)]);
    expect(store.overBudget.value).toBe(false);
  });
});

describe('settle', () => {
  it('resolves when nothing wanted is pending, or says it did not within a timeout', async () => {
    const { clock, store } = setup({ delayMs: 100 });
    store.want([view(0, 1024)]);
    const settled = store.whenSettled();
    const timedOut = store.whenSettled(50);
    await clock.advance(50);
    expect(await timedOut).toBe(false);
    await clock.advance(50);
    expect(await settled).toBe(true);
    expect(await store.whenSettled()).toBe(true);
  });
});

describe('nodes', () => {
  const T1 = { kind: 'sequence', id: 'T1', params: { length: 300 } };
  const withNode =
    (node = T1) =>
    (r: RegionRequest): RegionAnswer => {
      const a = rangesAnswer(r);
      return a.ok ? { ...a, nodes: [{ key: 'T1', node }] } : a;
    };

  it('are held by the regions that carried them and leave quietly with the last', async () => {
    const { clock, store } = setup({ delayMs: 1, answer: withNode() });
    store.want([view(0, 2048)]);
    await clock.advance(1);
    expect(store.getNode(S, 'T1')?.node).toBe(T1);
    expect(store.get(key(0))).toMatchObject({ nodes: [`${S}#T1`] });
    const release = store.hold(S, 'T1');
    store.detach(S);
    expect(store.getNode(S, 'T1')?.node).toBe(T1);
    release();
    release();
    expect(store.getNode(S, 'T1')).toBeUndefined();
    expect(store.issues.value).toEqual([]);
  });

  it('keep the held node on a conflict, with an issue; return under the same key', async () => {
    let node = T1;
    const { clock, store } = setup({ delayMs: 1, settings: { maxBytes: 1 }, answer: (r) => withNode(node)(r) });
    store.want([view(0, 1024)]);
    await clock.advance(1);
    node = { ...T1, params: { length: 301 } };
    store.want([view(0, 2048)]);
    await clock.advance(1);
    expect(store.getNode(S, 'T1')?.node).toBe(T1);
    expect(store.issues.value.map((i) => i.code)).toContain('seq.store.node_conflict');
    store.want([]);
    expect(store.getNode(S, 'T1')).toBeUndefined();
    store.want([view(0, 1024)]);
    await clock.advance(1);
    expect(store.getNode(S, 'T1')?.key).toBe('T1');
  });

  it('known up front are held while their source is attached', () => {
    const clock = new FakeClock();
    const store = new DataStore({ clock });
    const source = new SyntheticSource({ ...syntheticCoverage(), nodes: [{ key: 'T1', node: T1 }] }, clock);
    store.attach(S, source);
    expect(store.node(S, 'T1').value?.node).toBe(T1);
    store.detach(S);
    expect(store.getNode(S, 'T1')).toBeUndefined();
  });
});

describe('the store', () => {
  it('lets its sources go on detach and retire: reads aborted, entries gone, sources disposed', async () => {
    const { clock, source, store } = setup({ delayMs: 100 });
    store.want([view(0, 2048)]);
    store.retire();
    await flush();
    expect(source.reads.every((r) => r.aborted)).toBe(true);
    expect(source.disposed).toBe(true);
    expect(store.get(key(0))).toBeUndefined();
    expect(store.inFlight).toBe(0);
    expect(clock.pending).toBe(0);
  });
});
