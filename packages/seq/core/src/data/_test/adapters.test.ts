/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Seq data ports (_spec/README.md "Data", Protocol; "Seq io"): the default adapter (the `seq` runtime reading a file
// through the assets port) and a host's (a test structure host serving a description from memory, through
// `structureNodes`) pass one conformance suite, give equal nodes for one structure, and the host's run reads no asset;
// the host's edit is told by `changes` and reattached in place with node identity kept (EC-176, EC-185).

import type { Issue } from '@molstar/core/repr/tree';
import { type StructureDescription, structureNodes } from '@molstar/seqx-core/io/structure';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { buildOverlay, resolveSpec } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { Source, zincFinger, zincFingerNodes } from '../../io/_test/fixtures.js';
import { channelAddress } from '../channels.js';
import { dataKey } from '../keys.js';
import {
  type AssetReader,
  abortError,
  type FormatRuntime,
  type OpenedSource,
  type RegionAnswer,
  type SourceEntry,
} from '../protocol.js';
import { encodeSeq } from '../seq-format.js';
import { SeqRuntime } from '../seq-runtime.js';
import { DataStore } from '../store.js';
import { checkAdapter } from './conformance.js';
import { FakeClock } from './synthetic.js';

const File = 'entries/zf.seq';

/** The assets port with the structure's `seq` file, counting its reads. */
function fileAssets() {
  const text = encodeSeq({
    coverage: {
      systems: [],
      shapes: [],
      levels: [1],
      nodes: zincFingerNodes().map((node) => ({ key: node.id!, node })),
    },
    pairs: [],
  });
  const reads: string[] = [];
  const assets: AssetReader = {
    readText: (uri) => (reads.push(uri), uri === File ? Promise.resolve(text) : Promise.reject(new Error('none'))),
    readBytes: (uri) => (
      reads.push(uri),
      uri === File ? Promise.resolve(new TextEncoder().encode(text)) : Promise.reject(new Error('none'))
    ),
  };
  return { assets, reads };
}

/**
 * A test host holding one structure in memory (a description it edits): it serves the sources whose `uri` or `ref`
 * names it, the units' values by `field`, tells its edits; any other source goes to `fallback`.
 */
function testHost(fallback: FormatRuntime) {
  let description = zincFinger();
  const listeners = new Set<() => void>();
  const serve = (entry: SourceEntry): OpenedSource => {
    const read = structureNodes(description, { source: entry.address });
    const values = (field: string) => read.values(field);
    return {
      coverage: {
        systems: read.systems,
        shapes: read.systems.map((s) => ({ system: s.name, shape: 'values' as const })),
        levels: [1],
        nodes: read.nodes,
      },
      issues: read.issues as Issue[],
      read(request, signal): Promise<RegionAnswer> {
        if (signal.aborted) return Promise.reject(abortError());
        const all = values(request.field ?? 'b_factor')?.get(request.system);
        if (!all) return Promise.resolve({ ok: true, tables: {} });
        const value = all.slice(request.start, request.end);
        return Promise.resolve({ ok: true, tables: { values: { form: 'positions', start: request.start, value } } });
      },
      changes(fn) {
        listeners.add(fn);
        return () => listeners.delete(fn);
      },
      dispose() {},
    };
  };
  const runtime: FormatRuntime = {
    name: 'test-host',
    formats: ['seq'],
    open: (entry, ctx) =>
      entry.uri === File || entry.ref === 'zf' ? Promise.resolve(serve(entry)) : fallback.open(entry, ctx),
  };
  const edit = (next: StructureDescription = description) => {
    description = next;
    for (const fn of [...listeners]) fn();
  };
  return { runtime, edit };
}

/** The zinc finger without chain B (its residues and values with it). */
function withoutB(d: StructureDescription): StructureDescription {
  const rows = [...d.residues.chain.keys()].filter((i) => d.residues.chain[i] !== 1);
  const chain = (c: number) => (c > 1 ? c - 1 : c);
  const drop = <T>(xs: readonly T[]) => xs.filter((_, i) => i !== 1);
  return {
    ...d,
    chains: {
      label_asym_id: drop(d.chains.label_asym_id),
      auth_asym_id: drop(d.chains.auth_asym_id!),
      entity: d.chains.entity.filter((_, i) => i !== 1),
      operator: drop(d.chains.operator!),
    },
    residues: {
      chain: Int32Array.from(rows, (i) => chain(d.residues.chain[i])),
      position: Int32Array.from(rows, (i) => d.residues.position[i]),
      auth_seq_id: Int32Array.from(rows, (i) => d.residues.auth_seq_id![i]),
      ins_code: rows.map((i) => d.residues.ins_code![i]),
    },
    links: { ...d.links!, a_chain: Int32Array.of(0), b_chain: Int32Array.of(1) },
    values: { b_factor: Float32Array.from(rows, (i) => d.values!.b_factor[i]) },
  };
}

const entry: SourceEntry = { address: Source, format: 'seq', uri: File };

describe('the data port: one conformance suite', () => {
  it('the default adapter (the seq runtime, reading the file through the assets port) passes', async () => {
    await checkAdapter({ runtime: SeqRuntime, entry, assets: fileAssets().assets });
  });

  it('a host’s adapter (a description in memory, through structureNodes) passes, its edits told', async () => {
    const host = testHost(SeqRuntime);
    await checkAdapter({ runtime: host.runtime, entry, assets: fileAssets().assets, edit: () => host.edit() });
  });

  it('both give equal nodes, keys and content; the host’s run reads no asset', async () => {
    const signal = new AbortController().signal;
    const limits = { maxRows: 1e6, maxBytes: 1e8 };
    const file = fileAssets();
    const byFile = await SeqRuntime.open(entry, { assets: file.assets, signal, limits });
    const none = fileAssets();
    const byHost = await testHost(SeqRuntime).runtime.open(entry, { assets: none.assets, signal, limits });
    expect(byHost.coverage.nodes!.map((n) => n.key)).toEqual(byFile.coverage.nodes!.map((n) => n.key));
    for (let i = 0; i < byFile.coverage.nodes!.length; i++)
      expect(byHost.coverage.nodes![i].node).toEqual(byFile.coverage.nodes![i].node);
    expect(file.reads).toEqual([File]);
    expect(none.reads).toEqual([]);
  });

  it('a host hands a source it does not hold to its fallback, and a `ref` alone reaches it (EC-185)', async () => {
    const signal = new AbortController().signal;
    const limits = { maxRows: 1e6, maxBytes: 1e8 };
    const host = testHost(SeqRuntime);
    const other = fileAssets();
    await expect(
      host.runtime.open(
        { address: '/data/x', format: 'seq', uri: 'other.seq' },
        { assets: other.assets, signal, limits },
      ),
    ).rejects.toThrow();
    expect(other.reads).toEqual(['other.seq']);
    const byRef = await host.runtime.open(
      { address: Source, format: 'seq', ref: 'zf' },
      { assets: other.assets, signal, limits },
    );
    expect(byRef.coverage.nodes?.length).toBe(4);
  });
});

describe('a host’s edit through the store', () => {
  async function attached() {
    const clock = new FakeClock();
    const store = new DataStore({ clock });
    const host = testHost(SeqRuntime);
    const ctx = {
      assets: fileAssets().assets,
      signal: new AbortController().signal,
      limits: { maxRows: 1e6, maxBytes: 1e8 },
    };
    store.attach(Source, await host.runtime.open(entry, ctx));
    return { clock, store, host, reopen: () => host.runtime.open(entry, ctx) };
  }

  it('a chain removed: the source reattached, every other node kept as it was, the complex’s content new (EC-176)', async () => {
    const { store, host, reopen } = await attached();
    const sequence = store.getNode(Source, '1');
    const complex = store.getNode(Source, 'complex')!;
    host.edit(withoutB(zincFinger()));
    expect(store.reattach(Source, await reopen())).toBe(true);
    expect(store.getNode(Source, '1')).toBe(sequence);
    const units = store.getNode(Source, 'complex')!.node.children!.filter((c) => c.kind === 'unit');
    expect(units.map((u) => u.id)).toEqual(['A-1_555', 'C-1_555']);
    expect(store.getNode(Source, 'complex')!.version).toBeGreaterThan(complex.version);
    expect(store.coverage(Source)?.systems.map((s) => s.name)).toEqual(['A-1_555', 'C-1_555']);
    expect(store.issues.value).toEqual([]);
  });

  it('values by field from the one opened source: B-factors and another channel, read again after an edit', async () => {
    const { clock, store, host, reopen } = await attached();
    const bf = channelAddress(Source, 'b_factor');
    const want = { source: bf, system: 'A-1_555', start: 0, end: 5, level: 1 };
    store.want([want]);
    await clock.advance(0);
    const values = () => {
      const e = store.get(dataKey(bf, 'A-1_555', 1, 0));
      return e?.status === 'loaded' && e.tables.values?.form === 'positions'
        ? Array.from(e.tables.values.value)
        : e?.status;
    };
    expect(values()).toEqual([Number.NaN, 20, 21, 22, 23]);
    const d = zincFinger();
    host.edit({ ...d, values: { b_factor: d.values!.b_factor.map((v) => v + 1) } });
    store.reattach(Source, await reopen());
    expect(values()).toBe('pending');
    await clock.advance(0);
    expect(values()).toEqual([Number.NaN, 21, 22, 23, 24]);
  });
});

describe('a host’s malformed node', () => {
  it('is an issue at its address and left out; the host’s other nodes kept, nothing thrown', () => {
    const store = new DataStore({ clock: new FakeClock() });
    const b = createSeqSpec();
    b.data.source('zf', { format: 'seq', uri: File });
    const model = resolveSpec(b.build(), { runtimes: [SeqRuntime] });
    const read = structureNodes(zincFinger(), { source: Source });
    const bad = { key: 'bad', node: { kind: 'sequence', id: 'bad', params: { alphabet: 'klingon', letters: 'Q' } } };
    store.attach(Source, {
      coverage: { systems: [], shapes: [], levels: [1], nodes: [read.nodes[0], bad] },
      read: () => Promise.resolve({ ok: true, tables: {} }),
      dispose() {},
    });
    const overlay = buildOverlay(model, store);
    expect(overlay.nodes.map((n) => n.key)).toEqual(['1']);
    expect(overlay.issues.map((i) => [i.code, i.address])).toEqual([
      ['seq.resolve.node.invalid', `${Source}/(node :key bad)`],
    ]);
  });
});
