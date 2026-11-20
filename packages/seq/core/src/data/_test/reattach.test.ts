/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A new reading of an attached source (`reattach`, _spec/README.md "Data", Store) and a source's values channels
// (`RegionRequest.field`, Store, Channels): node identity kept by key, keys gone removed, tiles read again; a field's
// tiles, summaries and entries its own.

import { describe, expect, it } from 'vitest';
import { channelAddress, channelOf, channelSource, fieldChannel, runtimeFor, runtimesByFormat } from '../channels.js';
import { dataKey } from '../keys.js';
import type { Coverage, FormatRuntime, RegionAnswer, RegionRequest } from '../protocol.js';
import { DataStore } from '../store.js';
import { FakeClock, rangesAnswer, SyntheticSource, syntheticCoverage } from './synthetic.js';

const S = '/data/structure';
const unit = (length: number) => ({ kind: 'sequence', id: 'u', params: { alphabet: 'protein', length } });
const withNodes = (nodes: Record<string, number>): Coverage => ({
  ...syntheticCoverage(),
  nodes: Object.entries(nodes).map(([key, length]) => ({ key, node: unit(length) })),
});
const view = (source = S) => ({ source, system: 'chr1', start: 0, end: 1024, level: 1 });

function setup(nodes: Record<string, number>) {
  const clock = new FakeClock();
  const store = new DataStore({ clock });
  const first = new SyntheticSource(withNodes(nodes), clock, { delayMs: 1 });
  store.attach(S, first);
  return { clock, store, first };
}

describe('reattach: a new reading of a source in place', () => {
  it('keeps a node’s identity by key and takes its new content, a new version; an unchanged one keeps its entry (EC-176)', () => {
    const { clock, store } = setup({ A: 10, B: 20 });
    const a = store.getNode(S, 'A')!;
    const b = store.getNode(S, 'B')!;
    const nodeVersion = store.nodeVersion.value;
    expect(store.reattach(S, new SyntheticSource(withNodes({ A: 10, B: 21 }), clock))).toBe(true);
    expect(store.getNode(S, 'A')).toBe(a);
    expect(store.getNode(S, 'B')?.node.params?.length).toBe(21);
    expect(store.getNode(S, 'B')!.version).toBeGreaterThan(b.version);
    expect(store.nodeVersion.value).toBeGreaterThan(nodeVersion);
    // a new reading is no conflict: EC-176's rule is for two answers of one reading
    expect(store.issues.value).toEqual([]);
  });

  it('removes a key the new reading lacks, even while something holds it; the hold’s release is then a no-op', () => {
    const { clock, store } = setup({ A: 10, B: 20 });
    const release = store.hold(S, 'B');
    store.reattach(S, new SyntheticSource(withNodes({ A: 10 }), clock));
    expect(store.getNode(S, 'B')).toBeUndefined();
    release();
    expect(store.getNode(S, 'A')).toBeDefined();
    // the key back in a later reading is held by it alone
    store.reattach(S, new SyntheticSource(withNodes({ A: 10, B: 22 }), clock));
    expect(store.getNode(S, 'B')?.node.params?.length).toBe(22);
    store.detach(S);
    expect(store.nodes.size).toBe(0);
  });

  it('drops every tile and asks the wanted ones of the new reading; the old one disposed, its reads aborted', async () => {
    const { clock, store, first } = setup({ A: 10 });
    store.want([view()]);
    await clock.advance(1);
    expect(store.get(dataKey(S, 'chr1', 1, 0))?.status).toBe('loaded');
    store.want([{ ...view(), start: 1024, end: 2048 }, view()]);
    const second = new SyntheticSource(withNodes({ A: 10 }), clock, { delayMs: 1 });
    store.reattach(S, second);
    expect(first.disposed).toBe(true);
    expect(first.reads.some((r) => r.aborted)).toBe(true);
    expect(store.get(dataKey(S, 'chr1', 1, 0))?.status).toBe('pending');
    await clock.advance(1);
    expect(second.reads.map((r) => r.request.start)).toEqual([1024, 0]);
    expect(store.get(dataKey(S, 'chr1', 1, 0))?.status).toBe('loaded');
    expect(store.settled.value).toBe(true);
  });

  it('a node an old tile brought goes when the new reading’s tiles no longer bring it', async () => {
    const clock = new FakeClock();
    const store = new DataStore({ clock });
    let brings = true;
    const answer = (r: RegionRequest): RegionAnswer => {
      const a = rangesAnswer(r);
      return a.ok && brings ? { ...a, nodes: [{ key: 'T', node: unit(5) }] } : a;
    };
    store.attach(S, new SyntheticSource(syntheticCoverage(), clock, { delayMs: 1, answer }));
    store.want([view()]);
    await clock.advance(1);
    const release = store.hold(S, 'T');
    brings = false;
    store.reattach(S, new SyntheticSource(syntheticCoverage(), clock, { delayMs: 1, answer }));
    await clock.advance(1);
    expect(store.getNode(S, 'T')).toBeUndefined();
    release();
  });

  it('a held node an old tile brought keeps its identity while it is served again: kept as it was, or new content (EC-176)', async () => {
    const clock = new FakeClock();
    const store = new DataStore({ clock });
    let length = 5;
    const answer = (r: RegionRequest): RegionAnswer => {
      const a = rangesAnswer(r);
      return a.ok ? { ...a, nodes: [{ key: 'T', node: unit(length) }] } : a;
    };
    store.attach(S, new SyntheticSource(syntheticCoverage(), clock, { delayMs: 1, answer }));
    store.want([view()]);
    await clock.advance(1);
    const release = store.hold(S, 'T');
    const before = store.getNode(S, 'T');
    store.reattach(S, new SyntheticSource(syntheticCoverage(), clock, { delayMs: 1, answer }));
    // between the drop and the answer it stays: held, never gone for a frame
    expect(store.getNode(S, 'T')).toBe(before);
    await clock.advance(1);
    expect(store.getNode(S, 'T')).toBe(before);
    length = 6;
    store.reattach(S, new SyntheticSource(syntheticCoverage(), clock, { delayMs: 1, answer }));
    await clock.advance(1);
    expect(store.getNode(S, 'T')?.node.params?.length).toBe(6);
    expect(store.issues.value).toEqual([]);
    release();
  });

  it('reading a channel’s coverage makes nothing: the source’s own, until a tile of the channel is read', () => {
    const { store } = setup({ A: 10 });
    const plddt = channelAddress(S, 'plddt');
    expect(store.coverage(plddt)).toBe(store.coverage(S));
    expect(store.keysOf({ ...view(plddt) })).toHaveLength(1);
    expect(store.summaryBytes).toBe(0);
  });

  it('refuses a malformed coverage with an issue and keeps the old reading', () => {
    const { clock, store, first } = setup({ A: 10 });
    const bad = new SyntheticSource({ ...syntheticCoverage(), levels: [] }, clock);
    expect(store.reattach(S, bad)).toBe(false);
    expect(store.coverage(S)).toBe(first.coverage);
    expect(first.disposed).toBe(false);
    expect(store.issues.value.map((i) => i.code)).toEqual(['seq.source.malformed']);
  });

  it('of a source not attached attaches it', () => {
    const clock = new FakeClock();
    const store = new DataStore({ clock });
    expect(store.reattach(S, new SyntheticSource(withNodes({ A: 1 }), clock))).toBe(true);
    expect(store.getNode(S, 'A')).toBeDefined();
  });
});

describe('values channels (RegionRequest.field)', () => {
  const values: Coverage = {
    systems: [{ name: 'A', length: 8 }],
    shapes: [{ system: 'A', shape: 'values' }],
    levels: [1],
    nodes: [{ key: 'A', node: unit(8) }],
  };
  const answer = (r: RegionRequest): RegionAnswer => ({
    ok: true,
    tables: {
      values: {
        form: 'positions',
        start: r.start,
        value: new Float32Array(r.end - r.start).fill(r.field === 'plddt' ? 90 : 1),
      },
    },
    nodes: [{ key: 'A', node: unit(8) }],
  });
  const want = (source: string) => ({ source, system: 'A', start: 0, end: 8, level: 1 });

  it('name a source and a field by an address no spec address holds, read back', () => {
    const c = channelAddress(S, 'plddt');
    expect(channelOf(c)).toEqual({ source: S, field: 'plddt' });
    expect(channelSource(c)).toBe(S);
    expect(channelOf(S)).toBeUndefined();
    expect(channelSource(S)).toBe(S);
  });

  it('read one field each from a source opened once: tiles their own, nodes the source’s alone', async () => {
    const clock = new FakeClock();
    const store = new DataStore({ clock });
    const source = new SyntheticSource(values, clock, { delayMs: 1, answer });
    store.attach(S, source);
    const plddt = channelAddress(S, 'plddt');
    store.want([want(S), want(plddt)]);
    await clock.advance(1);
    expect(source.reads.map((r) => r.request.field)).toEqual([undefined, 'plddt']);
    const read = (address: string) => store.get(dataKey(address, 'A', 1, 0));
    const v = (address: string) => {
      const e = read(address);
      return e?.status === 'loaded' && e.tables.values?.form === 'positions' ? e.tables.values.value[0] : undefined;
    };
    expect([v(S), v(plddt)]).toEqual([1, 90]);
    expect([...store.nodes.keys()]).toEqual([`${S}#A`]);
    expect(store.coverage(plddt)?.systems).toEqual(values.systems);
    // gone with its source
    store.detach(S);
    expect(store.coverage(plddt)).toBeUndefined();
    expect(read(plddt)).toBeUndefined();
  });

  it('a channel’s failure is reported at its source’s address', async () => {
    const clock = new FakeClock();
    const store = new DataStore({
      clock,
      settings: { retry: { attempts: 1, baseMs: 1, factor: 1, cooldownMs: 1000 } },
    });
    const fails = (r: RegionRequest) => (r.field ? new Error('no such field') : answer(r));
    store.attach(S, new SyntheticSource(values, clock, { delayMs: 1, answer: fails }));
    store.want([want(channelAddress(S, 'nope'))]);
    await clock.advance(1);
    expect(store.issues.value.map((i) => [i.code, i.address])).toEqual([['seq.source.failed', S]]);
  });

  it('a channel strips nodes from its answers and disposes nothing of its source', async () => {
    const clock = new FakeClock();
    const source = new SyntheticSource(values, clock, { delayMs: 0, answer });
    const channel = fieldChannel(source, 'b_factor');
    expect(channel.coverage.nodes).toBeUndefined();
    const p = channel.read(
      { system: 'A', start: 0, end: 8, level: 1, shapes: ['values'] },
      new AbortController().signal,
    );
    await clock.advance(0);
    const a = await p;
    expect(a).toMatchObject({ ok: true });
    expect(a.ok && a.nodes).toBeUndefined();
    channel.dispose();
    expect(source.disposed).toBe(false);
  });
});

describe('runtimes by format', () => {
  const rt = (name: string, ...formats: string[]): FormatRuntime => ({
    name,
    formats,
    open: () => Promise.reject(new Error('unused')),
  });

  it('bind one runtime a format; several bind none, in any order', () => {
    const a = rt('a', 'x', 'y');
    const b = rt('b', 'y');
    const c = rt('c', 'y');
    for (const list of [
      [a, b, c],
      [c, b, a],
      [b, a, c],
    ]) {
      expect(runtimeFor(list, 'x')).toBe(a);
      expect(runtimeFor(list, 'y')).toBe('conflict');
      expect(runtimeFor(list, 'z')).toBeUndefined();
      expect([...(runtimesByFormat(list).conflicts.get('y') ?? [])].sort()).toEqual(['a', 'b', 'c']);
    }
    expect(runtimeFor([a, a], 'y')).toBe(a);
  });
});
