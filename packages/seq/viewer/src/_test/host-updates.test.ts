/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A host's updates (_spec/README.md "StoreModel", Host updates): an opened source's `changes` reopen it in place a
// frame later, changes inside one frame once; node identity kept by key, a key gone removed though held (EC-176); a
// reopen that fails is an issue and the reading before it stays; a source gone stops listening.

import type { Node } from '@molstar/core/repr/tree';
import type { FormatRuntime, OpenedSource } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { DomHost } from '@molstar/ui/dom';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { StoreModel } from '../models/store.js';

const assets = {
  readText: () => Promise.reject(new Error('no network in tests')),
  readBytes: () => Promise.reject(new Error('no network in tests')),
};
const unit = (length: number): Node => ({ kind: 'sequence', id: 'u', params: { alphabet: 'protein', length } });

/**
 * A host holding units by key, telling its edits (one since a reading was opened, at its subscription); `fail` makes
 * its next opening fail, `gate` holds an opening after it read the units, until it resolves.
 */
function host() {
  const state = {
    units: { A: 10, B: 20 } as Record<string, number>,
    fail: false,
    opens: 0,
    listeners: 0,
    version: 0,
    gate: undefined as Promise<void> | undefined,
  };
  const listeners = new Set<() => void>();
  const runtime: FormatRuntime = {
    name: 'host',
    formats: ['host'],
    async open() {
      state.opens++;
      if (state.fail) throw new Error('the host lost its model');
      const at = state.version;
      const nodes = Object.entries(state.units).map(([key, length]) => ({ key, node: unit(length) }));
      await state.gate;
      const opened: OpenedSource = {
        coverage: { systems: [], shapes: [], levels: [1], nodes },
        read: () => Promise.resolve({ ok: true, tables: {} }),
        changes(fn) {
          listeners.add(fn);
          state.listeners = listeners.size;
          if (state.version !== at) queueMicrotask(() => listeners.has(fn) && fn());
          return () => {
            listeners.delete(fn);
            state.listeners = listeners.size;
          };
        },
        dispose() {},
      };
      return opened;
    },
  };
  const edit = (units: Record<string, number>) => {
    state.units = units;
    state.version++;
    for (const fn of [...listeners]) fn();
  };
  return { runtime, state, edit };
}

/** The next frames, run by hand. */
function frames() {
  const queued: ((t: number) => void)[] = [];
  vi.spyOn(DomHost, 'requestAnimationFrame').mockImplementation((fn) => void queued.push(fn));
  return {
    get queued() {
      return queued.length;
    },
    async run() {
      for (const fn of queued.splice(0)) fn(0);
      await new Promise((r) => setImmediate(r));
    },
  };
}

function model(runtime: FormatRuntime, withSource = true) {
  const b = createSeqSpec();
  if (withSource) b.data.source('zf', { format: 'host', uri: 'zf.cif' });
  b.data.sequence('p', { alphabet: 'protein', letters: 'MEEP' });
  return resolveSpec(b.build(), { runtimes: [runtime] });
}

afterEach(() => vi.restoreAllMocks());

describe('a host’s updates', () => {
  it('reopen the source a frame later, once for changes inside one frame; identity kept by key (EC-176)', async () => {
    const h = host();
    const frame = frames();
    const m = new StoreModel();
    await m.sync(model(h.runtime), assets);
    const a = m.store.getNode('/data/zf', 'A');
    const release = m.store.hold('/data/zf', 'B');
    h.edit({ A: 10, B: 21 });
    h.edit({ A: 10, C: 5 });
    h.edit({ A: 10, C: 6 });
    expect(frame.queued).toBe(1);
    expect(h.state.opens).toBe(1);
    await frame.run();
    expect(h.state.opens).toBe(2);
    expect(m.store.getNode('/data/zf', 'A')).toBe(a);
    // gone from the new reading, though the selection held it: as when a patch removes it
    expect(m.store.getNode('/data/zf', 'B')).toBeUndefined();
    expect(m.store.getNode('/data/zf', 'C')?.node.params?.length).toBe(6);
    release();
    // the new reading is listened to, the old one not
    expect(h.state.listeners).toBe(1);
    h.edit({ A: 11 });
    await frame.run();
    expect(m.store.getNode('/data/zf', 'A')?.node.params?.length).toBe(11);
    expect(m.issues.value).toEqual([]);
    m.retire();
    expect(h.state.listeners).toBe(0);
  });

  it('a reopen that fails is an issue at the source, the reading before it kept; the next change tries again', async () => {
    const h = host();
    const frame = frames();
    const m = new StoreModel();
    await m.sync(model(h.runtime), assets);
    h.state.fail = true;
    h.edit({ A: 1 });
    await frame.run();
    expect(m.issues.value.map((i) => [i.code, i.address])).toEqual([['seq.source.unreadable', '/data/zf']]);
    expect(m.store.getNode('/data/zf', 'B')).toBeDefined();
    h.state.fail = false;
    h.edit({ A: 2 });
    await frame.run();
    expect(m.issues.value).toEqual([]);
    expect(m.store.getNode('/data/zf', 'B')).toBeUndefined();
    m.retire();
  });

  it('a source gone stops listening: a later change reopens nothing', async () => {
    const h = host();
    const frame = frames();
    const m = new StoreModel();
    await m.sync(model(h.runtime), assets);
    h.edit({ A: 1 });
    await m.sync(model(h.runtime, false), assets);
    expect(h.state.listeners).toBe(0);
    await frame.run();
    h.edit({ A: 2 });
    expect(frame.queued).toBe(0);
    expect(h.state.opens).toBe(1);
    m.retire();
  });

  it('an edit while the source is being opened, or reopened, is not lost: the latest reading is attached', async () => {
    const h = host();
    const frame = frames();
    const m = new StoreModel();
    let open!: () => void;
    h.state.gate = new Promise<void>((r) => (open = r));
    const syncing = m.sync(model(h.runtime), assets);
    await Promise.resolve();
    h.edit({ A: 1 });
    open();
    await syncing;
    await new Promise((r) => setImmediate(r));
    expect(m.store.getNode('/data/zf', 'A')?.node.params?.length).toBe(10);
    await frame.run();
    expect(m.store.getNode('/data/zf', 'A')?.node.params?.length).toBe(1);
    // during a reopen
    h.state.gate = new Promise<void>((r) => (open = r));
    h.edit({ A: 2 });
    await frame.run();
    h.edit({ A: 3 });
    open();
    await new Promise((r) => setImmediate(r));
    h.state.gate = undefined;
    await frame.run();
    expect(m.store.getNode('/data/zf', 'A')?.node.params?.length).toBe(3);
    m.retire();
  });
});
