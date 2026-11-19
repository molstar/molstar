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

  it('a malformed answer fails at once with a cooldown, never retried before it', async () => {
    const malformed: RegionAnswer = {
      ok: true,
      tables: { ranges: { form: 'list', start: Int32Array.of(1, 2), end: Int32Array.of(3) } },
    };
    const { clock, source, store } = setup({ delayMs: 10, answer: () => malformed });
    store.want([view(0, 1024)]);
    await clock.advance(10);
    expect(store.get(key(0))).toMatchObject({
      status: 'failed',
      kind: 'malformed',
      issue: { code: 'seq.source.malformed' },
    });
    await clock.advance(29_000);
    expect(source.reads.length).toBe(1);
    await clock.advance(1_000);
    expect(source.reads.length).toBe(2);
  });

  it('a result that is no promise, or an answer that throws when read, is malformed; nothing stays in flight', async () => {
    const clock = new FakeClock();
    const throwing = {
      ok: true,
      get tables(): never {
        throw new Error('boom');
      },
    };
    let n = 0;
    const odd = {
      coverage: syntheticCoverage(),
      dispose() {},
      read: () => (n++ === 0 ? null : Promise.resolve(throwing)) as unknown as Promise<RegionAnswer>,
    };
    const store = new DataStore({ clock });
    store.attach(S, odd);
    store.want([view(0, 2048)]);
    await flush();
    expect(store.get(key(0))).toMatchObject({ status: 'failed', kind: 'malformed' });
    expect(store.get(key(1))).toMatchObject({ status: 'failed', kind: 'malformed' });
    expect(store.inFlight).toBe(0);
    expect(store.settled.value).toBe(true);
  });

  it('too large is never retried', async () => {
    const tooLarge: RegionAnswer = {
      ok: false,
      kind: 'too-large',
      issue: issueOf(SourceIssues, 'seq.source.too_large', 'Too many'),
    };
    const { clock, source, store } = setup({ delayMs: 10, answer: () => tooLarge });
    store.want([view(0, 1024)]);
    await clock.advance(120_000);
    expect(source.reads.length).toBe(1);
    expect(store.get(key(0))).toMatchObject({ status: 'failed', kind: 'too-large' });
  });

  it('an AbortError the source raised itself is a failure', async () => {
    const own = new Error('gone');
    own.name = 'AbortError';
    const { clock, store } = setup({ delayMs: 10, answer: () => own });
    store.want([view(0, 1024)]);
    await clock.advance(10);
    expect(store.issues.value.map((i) => i.code)).toEqual(['seq.source.failed']);
    expect(store.get(key(0))?.status).toBe('pending');
    expect(store.inFlight).toBe(0);
  });

  it('keep their backoff across a move away and back: not asked at once, attempts kept', async () => {
    const { clock, source, store } = setup({ delayMs: 10, answer: () => failed() });
    store.want([view(0, 1024)]);
    await clock.advance(10);
    store.want([view(8192, 9216)]);
    store.want([view(0, 1024)]);
    expect(source.reads.filter((r) => r.request.start === 0).length).toBe(1);
    await clock.advance(500);
    expect(source.reads.filter((r) => r.request.start === 0).length).toBe(2);
    await clock.advance(10);
    store.want([]);
    await clock.advance(1000);
    store.want([view(0, 1024)]);
    // The second backoff (1,000) is over while away: asked again now, and its failure is the third (failed).
    expect(source.reads.filter((r) => r.request.start === 0).length).toBe(3);
    await clock.advance(10);
    expect(store.get(key(0))).toMatchObject({ status: 'failed', kind: 'failed' });
  });

  it('settle after a cooldown: a failed entry asked again is settled, and stays current', async () => {
    const { clock, store } = setup({
      delayMs: 10,
      settings: { retry: { attempts: 1, baseMs: 1, factor: 1, cooldownMs: 1000 } },
      answer: (r, n) => (n === 1 ? failed() : rangesAnswer(r)),
    });
    const settles: boolean[] = [];
    store.settled.listen((s) => settles.push(s));
    store.want([view(0, 1024)]);
    await clock.advance(10);
    expect(store.get(key(0))?.status).toBe('failed');
    await clock.advance(1000);
    expect(store.inFlight).toBe(1);
    expect(store.settled.value).toBe(true);
    await clock.advance(10);
    expect(store.get(key(0))?.status).toBe('loaded');
    expect(settles).toEqual([false, true]);
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
    expect(store.issues.value).toEqual([]);
  });

  it('never drops what is wanted: prefetch then margin become evicted, not asked again, without a loop', async () => {
    const probe = setup({ delayMs: 1 });
    probe.store.want([view(0, 1024)]);
    await probe.clock.advance(1);
    const b = tileBytes(probe.store);
    const { clock, source, store } = setup({ delayMs: 1, settings: { maxBytes: 2 * b + b / 2 } });
    const wants = [
      view(0, 1024),
      view(1024, 2048, { priority: 'margin' }),
      view(2048, 3072, { priority: 'prefetch' }),
      view(3072, 4096, { priority: 'margin' }),
    ];
    // A planner that wants again on every change.
    const stop = store.version.listen(() => store.want(wants));
    store.want(wants);
    await clock.advance(10);
    expect([0, 1, 2, 3].map((t) => store.get(key(t))?.status)).toEqual(['loaded', 'evicted', 'evicted', 'loaded']);
    expect(source.reads.length).toBe(4);
    await clock.advance(100);
    expect(source.reads.length).toBe(4);
    expect(store.settled.value).toBe(true);
    // Leaving and coming back asks again.
    stop();
    store.want([view(0, 1024)]);
    expect(store.get(key(2))).toBeUndefined();
  });

  it('never touches a request in flight', async () => {
    const { clock, source, store } = setup({ delayMs: 100, settings: { maxBytes: 1, graceMs: 1000 } });
    store.want([view(0, 1024)]);
    await clock.advance(100);
    store.want([view(1024, 2048)]);
    expect(store.get(key(0))).toBeUndefined();
    store.want([view(4096, 5120)]);
    await clock.advance(50);
    expect(source.reads[1].aborted).toBe(false);
    expect(store.get(key(1))?.status).toBe('pending');
  });
});

describe('want sets', () => {
  it('of an export merge with the live view’s and settle on their own (EC-125)', async () => {
    const { clock, source, store } = setup({ delayMs: 100 });
    store.want([view(0, 1024)]);
    await clock.advance(100);
    const exporter = store.wanter();
    exporter.want([view(0, 2048, { priority: 'prefetch' })]);
    expect(source.reads.length).toBe(2);
    expect([...exporter.keys]).toEqual([key(0), key(1)]);
    const done = exporter.whenSettled();
    await clock.advance(100);
    expect(await done).toBe(true);
    exporter.dispose();
    expect(store.get(key(0))?.status).toBe('loaded');
    expect(store.get(key(1))?.status).toBe('loaded');
    store.want([]);
    expect(store.get(key(1))?.status).toBe('loaded');
  });

  it('are bounded: lengths clip, a want over the tile cap is an issue and not asked', () => {
    const { source, store } = setup({ settings: { maxTilesPerWant: 4 } });
    store.want([view(999_000, 2_000_000)]);
    expect(source.reads.map((r) => [r.request.start, r.request.end])).toEqual([
      [998_400, 999_424],
      [999_424, 1_000_000],
    ]);
    store.want([view(0, 10 * 1024)]);
    expect(source.reads.length).toBe(2);
    expect(store.issues.value.map((i) => [i.code, i.address])).toEqual([['seq.store.too_many_tiles', S]]);
    // its keys are not built either: a reader of a long system at a fine level stops at the bound
    expect(store.keysOf(view(0, 4 * 1024))).toHaveLength(4);
    expect(store.keysOf(view(0, 10 * 1024))).toEqual([]);
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

describe('cancellation when the view waits (EC-173)', () => {
  it('lets the oldest stale request go at once, not after its grace', async () => {
    const { clock, source, store } = setup({ delayMs: 1000, settings: { maxInFlight: 2, graceMs: 5000 } });
    store.want([view(0, 2048)]);
    await clock.advance(10);
    store.want([view(1024, 2048)]);
    await clock.advance(10);
    store.want([view(8192, 9216)]);
    expect(source.reads.map((r) => [r.request.start, r.aborted])).toEqual([
      [0, true],
      [1024, false],
      [8192, false],
    ]);
    store.want([view(8192, 10240)]);
    expect(source.reads.map((r) => r.aborted)).toEqual([true, true, false, false]);
  });

  it('and a detach during a grace leaves no timer behind', async () => {
    const { clock, source, store } = setup({ delayMs: 1000 });
    store.want([view(0, 1024)]);
    store.want([]);
    store.detach(S);
    expect(clock.pending).toBe(0);
    await clock.advance(5000);
    expect(source.reads[0].aborted).toBe(true);
    expect(store.inFlight).toBe(0);
  });
});

describe('a detach with tiles queued and in flight', () => {
  it('leaves nothing pending, and a fresh attach and want settle promptly', async () => {
    const { clock, source, store } = setup({ delayMs: 100, settings: { maxInFlight: 1 } });
    const other = store.wanter();
    other.want([view(0, 3072)]);
    store.want([view(0, 3072)]);
    expect(source.reads.length).toBe(1);
    store.detach(S);
    // another set's wants (planned before, kept) merged again while it is gone make no entry nothing can answer
    store.want([view(0, 3072)]);
    expect([0, 1, 2].map((t) => store.get(key(t)))).toEqual([undefined, undefined, undefined]);
    expect(store.settled.value).toBe(true);
    expect(store.inFlight).toBe(0);
    expect(await other.whenSettled(0)).toBe(true);
    expect(clock.pending).toBe(0);

    const again = new SyntheticSource(syntheticCoverage(), clock, { delayMs: 100 });
    store.attach(S, again);
    store.want([view(0, 3072)]);
    other.want([view(0, 3072)]);
    expect(store.get(key(0))?.status).toBe('pending');
    const settled = other.whenSettled(1000);
    await clock.advance(300);
    expect(await settled).toBe(true);
    expect([0, 1, 2].map((t) => store.get(key(t))?.status)).toEqual(['loaded', 'loaded', 'loaded']);
    expect(store.settled.value).toBe(true);
  });
});

describe('coverage', () => {
  it('is checked on attach: a malformed one is an issue and the source unavailable', () => {
    const clock = new FakeClock();
    const store = new DataStore({ clock });
    const bad = (coverage: object) => new SyntheticSource(coverage as never, clock);
    expect(store.attach(S, bad({ ...syntheticCoverage(), levels: [10, 1] }))).toBe(false);
    expect(store.attach('/b', bad({ ...syntheticCoverage(), levels: [] }))).toBe(false);
    expect(store.attach('/c', bad({ ...syntheticCoverage(), levels: [0.5] }))).toBe(false);
    expect(store.attach('/d', bad({ ...syntheticCoverage(), systems: [{ name: 'chr1', length: -1 }] }))).toBe(false);
    expect(store.issues.value.map((i) => [i.code, i.address])).toEqual([
      ['seq.source.malformed', S],
      ['seq.source.malformed', '/b'],
      ['seq.source.malformed', '/c'],
      ['seq.source.malformed', '/d'],
    ]);
    expect(store.coverage(S)).toBeUndefined();
    store.want([view(0, 1024)]);
    expect(store.get(key(0))).toBeUndefined();
  });
});

describe('versions', () => {
  it('move only for the entries that changed', async () => {
    const { clock, store } = setup();
    store.want([view(0, 1024)]);
    await clock.advance(100);
    const first = store.get(key(0))!;
    const told: number[] = [];
    store.entry(key(0)).listen((e) => told.push(e!.version));
    store.want([view(0, 2048)]);
    await clock.advance(100);
    expect(store.get(key(1))?.status).toBe('loaded');
    expect(store.get(key(0))).toBe(first);
    expect(told).toEqual([]);
    expect(store.version.value).toBe(store.get(key(1))!.version);
  });
});

describe('issues', () => {
  it('are deduped, capped, and cleared with their source', async () => {
    const { clock, store } = setup({
      delayMs: 10,
      settings: { maxIssues: 1 },
      answer: (r) =>
        r.start === 0
          ? failed()
          : { ok: false, kind: 'not-served', issue: issueOf(SourceIssues, 'seq.source.not_served', 'No') },
    });
    store.want([view(0, 3072)]);
    await clock.advance(2000);
    expect(store.issues.value.map((i) => i.code)).toEqual(['seq.source.failed', 'seq.store.limit_exceeded']);
    store.detach(S);
    expect(store.issues.value).toEqual([]);
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
