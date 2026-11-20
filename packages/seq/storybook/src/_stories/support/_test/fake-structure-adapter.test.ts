/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import {
  type AssetReader,
  DefaultOpenLimits,
  type FormatRuntime,
  type OpenedSource,
  type SourceEntry,
} from '@molstar/seqx-core/data';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { readChainTable } from '@molstar/seqx-examples/p53/structure-record';
import { MmcifRuntime } from '@molstar/seqx-core/runtimes/mmcif';
import { describe, expect, it } from 'vitest';
import type { ChainTable } from '../fake-structure/chain-table.js';
import { FakeStructureViewer } from '../fake-structure/fake-viewer.js';
import { descriptionOf, entryOfUri, fakeStructureAdapter } from '../fake-structure-adapter.js';

// The fake host's adapter (seq-data-ports, The fake host's adapter): the viewer's table served through `structureNodes`
// gives the nodes the mmCIF runtime reads from the entry's file, keys and content, and reads no asset; its values by
// `field` are the file's; a chain dropped is told by `changes` and the other nodes are as they were; what it does not
// hold goes to its fallback. The default adapter's file is read here for the comparison (2OCJ: four chains and four
// zincs; 1FIN: two chain types in two copies, an ATP each).

const table = (file: string): ChainTable => readChainTable(readFileSync(`${PackageFolder}${file}`, 'utf8'));
const signal = new AbortController().signal;

/** The assets port reading the examples' fixtures, listing what it read. */
function assets() {
  const reads: string[] = [];
  const reader: AssetReader = {
    readText: async (uri) => (reads.push(uri), readFileSync(`${PackageFolder}${uri}`, 'utf8')),
    readBytes: async (uri) => (reads.push(uri), new Uint8Array(readFileSync(`${PackageFolder}${uri}`))),
  };
  return { reader, reads };
}

const entryOf = (id: string): SourceEntry => ({
  address: `/data/${id.toLowerCase()}`,
  format: 'mmcif',
  uri: `fixtures/structures/${id}.cif.gz`,
});

async function open(runtime: FormatRuntime, entry: SourceEntry, a = assets()) {
  const opened = await runtime.open(entry, { assets: a.reader, signal, limits: DefaultOpenLimits });
  return { opened, reads: a.reads };
}

const read = async (s: OpenedSource, system: string, field?: string) => {
  const r = await s.read(
    { system, start: 0, end: 1e6, level: 1, shapes: ['values'], ...(field ? { field } : {}) },
    signal,
  );
  return r.ok && r.tables.values?.form === 'positions' ? Array.from(r.tables.values.value) : r;
};

describe('the fake host’s adapter: one path', () => {
  const p53 = table('fixtures/p53/2OCJ-chains.json');
  const fake = new FakeStructureViewer(p53);
  const host = fakeStructureAdapter(fake, MmcifRuntime);

  it('2OCJ: the nodes the fake serves are those the file gives (keys and content), and the fake reads no asset', async () => {
    const byFile = await open(MmcifRuntime, entryOf('2OCJ'));
    const byHost = await open(host, entryOf('2OCJ'), assets());
    expect(byFile.reads).toEqual(['fixtures/structures/2OCJ.cif.gz']);
    expect(byHost.reads).toEqual([]);
    const a = byFile.opened.coverage;
    const b = byHost.opened.coverage;
    expect(b.nodes!.map((n) => n.key)).toEqual(a.nodes!.map((n) => n.key));
    for (let i = 0; i < a.nodes!.length; i++) expect(b.nodes![i].node).toEqual(a.nodes![i].node);
    expect(b.systems).toEqual(a.systems);
    // the channels it serves are the file's (it holds no occupancy: the file lists one more)
    expect(b.shapes.map((x) => [x.system, x.shape])).toEqual(a.shapes.map((x) => [x.system, x.shape]));
    for (const [i, x] of b.shapes.entries()) for (const f of x.fields!) expect(a.shapes[i].fields).toContainEqual(f);
    expect(b.levels).toEqual(a.levels);
    expect(byHost.opened.issues).toEqual(byFile.opened.issues);
  });

  it('values by field: the file’s B-factors per unit, the first channel the table has the default', async () => {
    const byFile = await open(MmcifRuntime, entryOf('2OCJ'));
    const byHost = await open(host, entryOf('2OCJ'));
    for (const system of ['A', 'D']) {
      const want = (await read(byFile.opened, system, 'b_factor')) as number[];
      const got = (await read(byHost.opened, system, 'b_factor')) as number[];
      expect(got).toHaveLength(want.length);
      // the table keeps the mean as the file's double, the runtime as a float32: the same to float32's digits
      got.forEach((v, i) => (Number.isNaN(want[i]) ? expect(v).toBeNaN() : expect(v).toBeCloseTo(want[i], 3)));
      expect(await read(byHost.opened, system)).toEqual(got);
    }
    const missing = await read(byHost.opened, 'A', 'plddt');
    expect(missing).toMatchObject({ ok: false, kind: 'not-served' });
    expect(await read(byHost.opened, 'Q')).toMatchObject({ ok: false, kind: 'not-served' });
  });

  it('a pLDDT column of a table is served by name, and is the default channel before B-factors', async () => {
    const t: ChainTable = {
      entry: 'AF-TEST',
      chains: [
        {
          label_asym_id: 'A',
          auth_asym_id: 'A',
          entity: '1',
          operators: ['1_555'],
          residues: [1, 2, 3].map((n) => ({
            label: n,
            auth: String(n),
            ins: '',
            comp: 'ALA',
            observed: true,
            b: n,
            plddt: 90 - n,
          })),
        },
      ],
    };
    const adapter = fakeStructureAdapter(new FakeStructureViewer(t), MmcifRuntime);
    const { opened } = await open(adapter, { address: '/data/af', format: 'mmcif', ref: 'AF-TEST' });
    expect(opened.coverage.shapes[0].fields?.map((f) => f.name)).toEqual(['plddt', 'b_factor']);
    expect(await read(opened, 'A')).toEqual([89, 88, 87]);
    expect(await read(opened, 'A', 'b_factor')).toEqual([1, 2, 3]);
    expect(opened.coverage.nodes!.map((n) => [n.key, n.node.kind])).toEqual([
      ['1', 'sequence'],
      ['complex', 'complex'],
    ]);
  });

  it('a source it does not hold, one asking for an assembly, goes to the fallback; a `ref` alone is held', async () => {
    const fallback: string[] = [];
    const spy: FormatRuntime = {
      name: 'spy',
      formats: ['mmcif'],
      open: (e) => (fallback.push(e.address), Promise.reject(new Error('fallback'))),
    };
    const h = fakeStructureAdapter(new FakeStructureViewer(p53), spy);
    const ctx = { assets: assets().reader, signal, limits: DefaultOpenLimits };
    await expect(
      h.open({ address: '/data/x', format: 'mmcif', uri: 'fixtures/structures/1FIN.cif.gz' }, ctx),
    ).rejects.toThrow('fallback');
    await expect(h.open({ ...entryOf('2OCJ'), assembly: '1' } as SourceEntry, ctx)).rejects.toThrow('fallback');
    await expect(h.open({ address: '/data/x', format: 'mmcif', value: {} as never }, ctx)).rejects.toThrow('fallback');
    expect(fallback).toEqual(['/data/x', '/data/2ocj', '/data/x']);
    const byRef = await h.open({ address: '/data/y', format: 'mmcif', ref: '2ocj' }, ctx);
    expect(byRef.coverage.nodes!.map((n) => n.key)).toContain('complex');
    expect(entryOfUri('https://files.rcsb.org/download/2OCJ.cif.gz?x=1')).toBe('2ocj');
    expect(entryOfUri('a/b/2ocj.bcif')).toBe('2ocj');
  });

  it('a malformed table is the converter’s issue, never a throw: a chain with no entity row left out', async () => {
    const t = table('fixtures/p53/2OCJ-chains.json');
    const bad: ChainTable = { ...t, model: { ...t.model!, entities: t.model!.entities.filter((e) => e.id !== '2') } };
    const { opened } = await open(fakeStructureAdapter(new FakeStructureViewer(bad), MmcifRuntime), entryOf('2OCJ'));
    // the zinc's entity is not described: its chains are left out, the protein stays
    expect(opened.coverage.systems.map((s) => s.name)).toEqual(['A', 'B', 'C', 'D']);
  });
});

describe('the fake host’s edits', () => {
  it('dropping a chain: told by `changes` once an edit; the others as they were, the chain and its zinc gone', async () => {
    const fake = new FakeStructureViewer(table('fixtures/p53/2OCJ-chains.json'));
    const host = fakeStructureAdapter(fake, MmcifRuntime);
    const first = (await open(host, entryOf('2OCJ'))).opened;
    let told = 0;
    const off = first.changes!(() => told++);
    expect(first.coverage.systems.map((s) => s.name)).toEqual(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']);
    expect(fake.dropChain('B')).toBe(true);
    expect(told).toBe(1);
    expect(fake.dropChain('B')).toBe(false);
    expect(told).toBe(1);
    const second = (await open(host, entryOf('2OCJ'))).opened;
    expect(second.coverage.systems.map((s) => s.name)).toEqual(['A', 'C', 'D', 'E', 'G', 'H']);
    // every node that is not the complex is as it was; the complex lost the units of B and F, and their links
    const was = new Map(first.coverage.nodes!.map((n) => [n.key, n.node]));
    for (const n of second.coverage.nodes!) if (n.key !== 'complex') expect(n.node).toEqual(was.get(n.key));
    const units = (o: OpenedSource) =>
      o.coverage.nodes!.find((n) => n.key === 'complex')!.node.children!.filter((c) => c.kind === 'unit');
    expect(units(second).map((u) => u.id)).toEqual(['A', 'C', 'D', 'E', 'G', 'H']);
    for (const u of units(second)) expect(u).toEqual(units(first).find((x) => x.id === u.id));
    off();
    fake.dropChain('C');
    expect(told).toBe(1);
  });

  it('a reading opened before an edit and listened to after it is told once, soon after', async () => {
    const fake = new FakeStructureViewer(table('fixtures/p53/2OCJ-chains.json'));
    const { opened } = await open(fakeStructureAdapter(fake, MmcifRuntime), entryOf('2OCJ'));
    fake.dropChain('D');
    let told = 0;
    const off = opened.changes!(() => told++);
    expect(told).toBe(0);
    await Promise.resolve();
    expect(told).toBe(1);
    off();
  });

  it('the viewer drops the chain of the selection, else the last polymer’s, and forgets what pointed at it', () => {
    const fake = new FakeStructureViewer(table('fixtures/p53/2OCJ-chains.json'));
    expect(fake.chainToDrop()).toBe('D');
    fake.click({ row: 1, index: 10 });
    expect(fake.chainToDrop()).toBe('B');
    fake.dropChain('B');
    expect(fake.selection.value).toEqual([]);
    expect(fake.rows.map((r) => r.asym)).toEqual(['A', 'C', 'D', 'E', 'G', 'H']);
    expect(fake.rows.map((r) => r.index)).toEqual([0, 1, 2, 3, 4, 5]);
    expect(descriptionOf(fake.table).chains.label_asym_id).toEqual(['A', 'C', 'D', 'E', 'G', 'H']);
  });
});

describe('the fake host’s adapter on 1FIN', () => {
  it('two chain types, copies and a ligand each: the nodes (keys and content) are the file’s, no asset read', async () => {
    const fake = new FakeStructureViewer(table('fixtures/complex/1FIN-chains.json'));
    const byFile = await open(MmcifRuntime, entryOf('1FIN'));
    const byHost = await open(fakeStructureAdapter(fake, MmcifRuntime), entryOf('1FIN'));
    expect(byHost.reads).toEqual([]);
    const a = byFile.opened.coverage;
    const b = byHost.opened.coverage;
    expect(b.nodes!.map((n) => n.key)).toEqual(a.nodes!.map((n) => n.key));
    for (let i = 0; i < a.nodes!.length; i++) expect(b.nodes![i].node).toEqual(a.nodes![i].node);
    expect(b.systems).toEqual(a.systems);
    expect(byHost.opened.issues).toEqual(byFile.opened.issues);
  });
});
