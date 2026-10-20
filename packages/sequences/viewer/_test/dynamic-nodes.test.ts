/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1b-a: dynamic nodes (seqspec:sources#dynamic-nodes; seqspec:viewer#updates-and-status,
// seqspec:sources#the-viewers-side), headless over the streamed genome. A manifest's nodes are in the set from the
// start; a patch's arrive with its region, are checked as document nodes (a bad one left out with its issue), held once
// (the same key again a no-op, other params a conflict), and reach the set through `#assemble.stream` once a frame; a
// track naming only dynamic nodes not held is hidden until one arrives; eviction takes the nodes nobody holds with
// their regions, quietly, and a node that returns gets its id back; the document,
// its history and the person's edits never see any of it; and data mode validates views naming dynamic nodes.

import type { StreamRequest } from '@molstar/seq-core/format';
import { KeyKind, type KeyTarget, OrderedSet, type SeqSet } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { GRCH38 } from '../../core/_test/format/stream-provider.js';
import type { Streams } from '../src/streams.js';
import {
  CacheBytes,
  dynamicDoc,
  GENES,
  GENOME,
  MaxInFlight,
  openGenome,
  PEAK_10MB,
  peakNode,
  type StreamedGenome,
  withPeaks,
} from './stream-fixtures.js';

const flush = () => new Promise((ok) => setTimeout(ok, 0));
const THREE = GRCH38.slice(0, 3);
const CHR1 = 0;
const t0: KeyTarget = { kind: -1, id: -1, row: -1 };
const idOf = (set: SeqSet | undefined, key: string) => (set?.byName(key, t0) ? t0.id : -1);
const setOf = (t: StreamedGenome) => t.ctx.graph.assembled.value?.set;

/** Block b's window [start, stop), its regions answered, and the frame that applies their nodes run. */
async function visit(t: StreamedGenome, b: number, start: number, stop: number) {
  t.vp(GENOME).zoomTo(b, start, stop);
  t.layout(GENOME);
  t.manager.pass();
  for (let k = 0; k < 4; k++) {
    await flush();
    t.manager.pass();
  }
  t.scheduler.advance(0);
  await t.ctx.graph.update();
}
/** Whether the view shows a track (not hidden), by ref. */
function shown(t: StreamedGenome, track: string): boolean {
  const v = t.ctx.graph.assembled.value?.view;
  const rv = v?.resolved;
  for (let k = 0; v && rv && k < rv.tracks.count; k++)
    if (v.set.strings.get(rv.tracks.ref[k]) === track) return rv.tracks.hidden[k] === 0;
  return false;
}
const open = (o: Parameters<typeof openGenome>[0] = {}) =>
  openGenome({ chromosomes: THREE, doc: dynamicDoc(), nodes: [GENES], answer: withPeaks(), ...o });
const codes = (t: StreamedGenome) => t.ctx.issues.current.value.map((i) => i.code);

describe('dynamic nodes in, from a manifest and with a patch', () => {
  it("the manifest's are in the set up front; a patch's arrive with its region, its named track shown then", async () => {
    const t = await open();
    const set0 = setOf(t);
    expect(set0?.byName('phylop_src/genes', t0) && t0.kind).toBe(KeyKind.annotation);
    expect([shown(t, 'genes'), shown(t, 'peak')]).toEqual([true, false]);
    // Pending, never unplaced.
    expect(codes(t)).not.toContain('svs.layer.unplaced');
    const runs = t.ctx.graph.runs;
    await visit(t, CHR1, 10_000_000, 10_010_000);
    const set = setOf(t);
    expect(set?.lineage).toBe(set0?.lineage);
    expect(idOf(set, PEAK_10MB)).toBeGreaterThanOrEqual(0);
    expect(idOf(set, 'phylop_src/genes')).toBe(idOf(set0, 'phylop_src/genes'));
    expect(shown(t, 'peak')).toBe(true);
    const v = t.ctx.graph.assembled.value?.view;
    const l = v
      ? [...Array(v.resolved.layers.count).keys()].find((k) => v.set.strings.get(v.resolved.layers.ref[k]) === 'peak.l')
      : -1;
    expect(v && l !== undefined && l >= 0 && v.resolved.layers.placementCount[l]).toBe(1);
    // One graph update for every node of the visit (its regions' peaks), not one a patch.
    expect((t.manager as Streams).nodes.entries.size).toBeGreaterThan(1);
    expect(t.ctx.graph.ran.has('#assemble.stream')).toBe(true);
    expect(t.ctx.graph.runs - runs).toBeLessThan(40);
    t.ctx.dispose();
  });

  it('held once: the same key again is a no-op, with other params a conflict and the held node kept', async () => {
    let other = false;
    const t = await open({
      answer: withPeaks(() => (other ? [peakNode('chr1', 10_000_000, { label: 'another' })] : [])),
    });
    await visit(t, CHR1, 10_000_000, 10_010_000);
    const set = setOf(t);
    const m = t.manager as Streams;
    const version = m.nodes.version;
    // Every region of the window carries the same node again: equal, nothing changes.
    other = false;
    await visit(t, CHR1, 10_020_000, 10_030_000);
    await visit(t, CHR1, 10_000_000, 10_010_000);
    expect(setOf(t)?.annotations[idOf(set, PEAK_10MB)]).toBe(set?.annotations[idOf(set, PEAK_10MB)]);
    expect(codes(t)).not.toContain('svs.stream.conflict');
    other = true;
    await visit(t, CHR1, 12_000_000, 12_010_000);
    expect(codes(t)).toContain('svs.stream.conflict');
    expect(m.nodes.entries.get(PEAK_10MB)?.node.params.label).toBe('peak 2441');
    expect(m.nodes.version).toBeGreaterThan(version);
    t.ctx.dispose();
  });

  it('checked as document nodes: a bad param or a ref to nothing leaves the node out, a ref into the same patch is in', async () => {
    const bad = (r: StreamRequest) => [
      { kind: 'annotation' as const, key: 'typo', params: { geometry: 'ranges', target: 'hg38/chr1', valeu: {} } },
      {
        kind: 'annotation' as const,
        key: 'lost',
        params: { geometry: 'ranges', target: 'hg38/chr9', value: { start: [0], stop: [5] } },
      },
      { kind: 'sequence' as const, key: `p-${r.start}`, params: { value: 'MKVLAAGG' } },
      {
        kind: 'annotation' as const,
        key: `on-p-${r.start}`,
        params: { geometry: 'ranges', target: `phylop_src/p-${r.start}/unit`, value: { start: [1], stop: [4] } },
      },
    ];
    const t = await open({ answer: withPeaks(bad) });
    await visit(t, CHR1, 10_000_000, 10_010_000);
    const set = setOf(t);
    const issues = t.ctx.issues.current.value;
    expect(issues.find((i) => i.code === 'svs.params.unknown')?.message).toContain('phylop_src');
    expect(issues.find((i) => i.code === 'svs.stream.ref')?.message).toContain('"hg38/chr9"');
    expect([idOf(set, 'phylop_src/typo'), idOf(set, 'phylop_src/lost')]).toEqual([-1, -1]);
    const start = Math.floor(10_000_000 / 4096) * 4096;
    expect(set?.byName(`phylop_src/p-${start}`, t0) && t0.kind).toBe(KeyKind.molecule);
    expect(set?.byName(`phylop_src/on-p-${start}`, t0) && t0.kind).toBe(KeyKind.annotation);
    t.ctx.dispose();
  });
});

describe('dynamic nodes out: eviction', () => {
  it('nodes nobody holds go with their regions, quietly; a shown one stays; one that returns gets its id back', async () => {
    const t = await open();
    t.ctx.config.set(MaxInFlight, 100);
    await visit(t, CHR1, 10_000_000, 10_010_000);
    const other = 'phylop_src/peak-chr1-2442';
    const set = setOf(t);
    const ids = [idOf(set, PEAK_10MB), idOf(set, other)];
    expect(ids.every((id) => id >= 0)).toBe(true);
    const before = t.ctx.issues.current.value.length;
    t.ctx.config.set(CacheBytes, 150_000);
    for (const at of [60, 110, 160]) await visit(t, 1, at * 1e6, at * 1e6 + 10_000);
    const out = setOf(t);
    // The shown peak is held by its placement; its neighbor, nobody's, left with its region and nothing was reported.
    expect([idOf(out, PEAK_10MB), idOf(out, other)]).toEqual([ids[0], -1]);
    expect(idOf(out, 'phylop_src/genes')).toBeGreaterThanOrEqual(0);
    expect(t.ctx.issues.current.value.slice(before).map((i) => i.code)).toEqual([]);
    expect(shown(t, 'peak')).toBe(true);
    t.ctx.config.set(CacheBytes, 64 * 2 ** 20);
    await visit(t, CHR1, 10_000_000, 10_010_000);
    expect([idOf(setOf(t), PEAK_10MB), idOf(setOf(t), other)]).toEqual(ids);
    t.ctx.dispose();
  });

  it('a node the selection names stays while selected, and goes with the next eviction once it is not', async () => {
    const t = await open();
    t.ctx.config.set(MaxInFlight, 100);
    await visit(t, CHR1, 10_000_000, 10_010_000);
    const other = 'phylop_src/peak-chr1-2442';
    const set = setOf(t) as SeqSet;
    const id = idOf(set, other);
    const items = { kind: 'item-loci' as const, set, elements: [{ annotation: id, items: OrderedSet.ofSingleton(0) }] };
    t.ctx.selection.set(items as never, 'api');
    t.ctx.config.set(CacheBytes, 150_000);
    for (const at of [60, 110]) await visit(t, 1, at * 1e6, at * 1e6 + 10_000);
    expect(idOf(setOf(t), other)).toBe(id);
    t.ctx.selection.clear('api');
    t.scheduler.advance(0);
    await t.ctx.graph.update();
    expect(idOf(setOf(t), other)).toBe(-1);
    t.ctx.dispose();
  });

  /** Views visited until `n` node patches have landed, under a 400 KB cache: the document and edits unchanged, the cache bounded. */
  async function nodePatches(n: number) {
    const t = await open();
    t.ctx.config.set(CacheBytes, 400_000);
    const doc = {
      base: t.ctx.document.base.value,
      shown: t.ctx.document.shown.value,
      local: t.ctx.document.local.value,
      revision: t.ctx.document.revision.value,
    };
    const m = t.manager as Streams;
    let patches = 0;
    const seen = new Set<string>();
    for (let k = 0; patches < n; k++) {
      const at = (k * 1_000_003) % 40_000_000;
      await visit(t, k % 3, at, at + 10_000);
      for (const ref of m.nodes.entries.keys()) seen.add(ref);
      patches = t.manager.stats.current.fine + t.manager.stats.current.coarse;
    }
    expect(seen.size).toBeGreaterThan(n / 2);
    expect({
      base: t.ctx.document.base.value,
      shown: t.ctx.document.shown.value,
      local: t.ctx.document.local.value,
      revision: t.ctx.document.revision.value,
    }).toEqual(doc);
    expect(t.ctx.document.shown.value).toBe(doc.shown);
    expect(t.ctx.history.canUndo.value).toBe(false);
    // Held under the bound: what the set has of the stream is what the manager holds, plus the manifest's.
    const set = setOf(t);
    let dynamic = 0;
    for (const a of set?.annotations ?? []) if (a && set?.strings.get(a.key).startsWith('phylop_src/')) dynamic++;
    expect(dynamic).toBe(m.nodes.entries.size + 1);
    expect(m.stats.current.bytes).toBeLessThanOrEqual(400_000 + 3 * 32_768 + 6_144);
    t.ctx.dispose();
  }

  it("two hundred node patches leave the document, its history and the person's edits untouched", () =>
    nodePatches(200));

  it(
    "a thousand node patches leave the document, its history and the person's edits untouched",
    { tags: ['heavy'] },
    () => nodePatches(1000),
  );
});

describe('data mode', () => {
  it("views naming dynamic nodes of the data's streamed sources validate; other refs do not", async () => {
    const t = await open();
    const set = setOf(t) as SeqSet;
    const ok = await t.ctx.setData(
      { set },
      { views: (b) => b.view({ ref: 'v', label: 'v' }).track().ranges({ data: 'phylop_src/peak-chr2-7' }) },
    );
    expect(ok.issues.filter((i) => i.severity === 'error')).toEqual([]);
    const no = await t.ctx.setData(
      { set },
      { views: (b) => b.view({ ref: 'v', label: 'v' }).track().ranges({ data: 'elsewhere/x' }) },
    );
    expect(no.issues.map((i) => i.code)).toContain('svs.ref.unresolved');
    t.ctx.dispose();
  });
});
