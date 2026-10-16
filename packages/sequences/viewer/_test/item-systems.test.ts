/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST7: item coordinate systems (seqspec:sources#data-patches; seqspec:sources#dynamic-nodes). A streamed transcript
// whose rows carry a `map` becomes, when its region arrives, a coordinate system of its own (`genes_src/T1`: a
// letterless molecule and unit) mapped onto chr2; a protein mapped onto it shows its domains along the transcript's
// exons in a genome view, and a section of the transcript's own residues shows them with the genome's data.

import { ItemSystemKind, ItemSystemProp, type StreamNode } from '@molstar/seq-core/format';
import { OrderedSet, type SeqLoci, type SeqSet } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import {
  CODING_CDS,
  GENE_LEVELS,
  GRCH38,
  hg38Streams,
  itemsManifest,
  type Scripted,
  scripted,
  TRANSCRIPT,
  transcriptDoc,
  transcriptItems,
  VOCABULARY,
  variantAt,
} from '../../core/_test/format/stream-provider.js';
import { DynamicNodes, heldNodeRefs } from '../src/dynamic.js';
import { ChannelMarks, createSeqViewer, Hit, MarkBit } from '../src/index.js';
import { DefaultSeqViewerSpec } from '../src/default-spec.js';
import { itemSystemNodes } from '../src/item-systems.js';
import type { Features } from '../src/prepare.js';
import type { Streams } from '../src/streams.js';
import { composed } from './composer-fixtures.js';
import { CacheBytes } from './stream-fixtures.js';

const flush = () => new Promise((ok) => setTimeout(ok, 0));
const THREE = GRCH38.slice(0, 3);
const KEYS = THREE.map((c) => c[0]);

type DocOptions = NonNullable<Parameters<typeof transcriptDoc>[0]>;
async function opened(
  o: DocOptions & {
    map?: Record<string, unknown> | null;
    cds?: [number, number][];
    section?: string;
    nodes?: StreamNode[];
  } = {},
) {
  const s = hg38Streams({ chromosomes: THREE, gated: true, letters: true, nodes: o.nodes });
  const genes: Scripted = scripted(itemsManifest(THREE, GENE_LEVELS, VOCABULARY.genes), (r) =>
    transcriptItems(r, KEYS, o.map ?? null, o.cds),
  );
  genes.gated = true;
  const stream = (src: { stream: string | null }) => (src.stream === 'genes' ? genes.provider : s.stream(src));
  const loaders = { load: () => Promise.reject(new Error('none')), stream };
  const t = await composed(transcriptDoc(o), {
    width: 1000,
    height: 400,
    loaders,
    section: o.section ?? 'tx.genome',
  });
  const m = t.ctx.streams as Streams;
  const drain = async () => {
    for (let k = 0; k < 8; k++) {
      for (const sec of t.ctx.viewports.sections.value) t.layout(sec);
      m.pass();
      for (const p of [s.genome, s.phylop, genes, s.variants]) p.release();
      await flush();
      await t.ctx.graph.update();
    }
  };
  await drain();
  return { t, s, genes, m, drain };
}

/** Layer `ref`'s prepared features (a `ranges` layer's), and the resolved view. */
function featuresOf(t: Awaited<ReturnType<typeof opened>>['t'], ref: string) {
  const v = t.ctx.graph.assembled.value?.view;
  const r = v?.resolved;
  for (let k = 0; v && r && k < r.layers.count; k++)
    if (v.set.strings.get(r.layers.ref[k]) === ref) return { f: v.prepared[k].state as Features, k, v };
  throw new Error(`no layer ${ref}`);
}
/** The chr2 positions of P1's residues [a, e), through its CDS (plus strand), as runs. */
function exonRuns(a: number, e: number): [number, number][] {
  const out: [number, number][] = [];
  let base = 0;
  for (const [s, t] of CODING_CDS) {
    const lo = Math.max(3 * a, base);
    const hi = Math.min(3 * e, base + t - s);
    if (hi > lo) out.push([s + lo - base, s + hi - base]);
    base += t - s;
  }
  return out;
}

const unitOf = (set: SeqSet, key: string) => Array.from(set.units.key).indexOf(set.strings.find(key));
/** A key's id in the current set; -1: none. */
function idOf(t: Awaited<ReturnType<typeof opened>>['t'], ref: string): number {
  const key = { kind: -1, id: -1, row: -1 };
  return t.ctx.graph.assembled.value?.set.byName(ref, key) ? key.id : -1;
}
function lengthOf(t: Awaited<ReturnType<typeof opened>>['t'], ref: string): number {
  const set = t.ctx.graph.assembled.value?.set;
  const id = idOf(t, ref);
  return id >= 0 ? (set?.molecules[id]?.length ?? -1) : -1;
}
const codes = (t: Awaited<ReturnType<typeof opened>>['t']) => t.ctx.issues.current.value.map((i) => i.code);
/** P1's domains' fragments on chr2's block of the genome section, sorted. */
function domainsOnChr2(t: Awaited<ReturnType<typeof opened>>['t']): [number, number][] {
  const { f, v } = featuresOf(t, 'tx.ranges-2');
  const chr2 = v.resolved.blocks.unit.indexOf(v.set.units.key.indexOf(v.set.strings.find('hg38/chr2')));
  const frags: [number, number][] = [];
  for (let q = 0; f && q < f.count; q++) if (f.block[q] === chr2) frags.push([f.start[q], f.stop[q]]);
  return frags.sort((a, b) => a[0] - b[0]);
}

describe('item coordinate systems', () => {
  it("a transcript's map becomes its own coordinate system, mapped onto chr2", async () => {
    const g = await opened();
    const set = g.t.ctx.graph.assembled.value?.set;
    const key = { kind: -1, id: -1, row: -1 };
    expect(set?.byName(`genes_src/${TRANSCRIPT}`, key)).toBe(true);
    const mol = set?.molecules[key.id];
    expect([mol?.length, mol?.lazy, mol ? set?.strings.get(mol.label) : '']).toEqual([1334, true, TRANSCRIPT]);
    expect(mol?.custom.props.get(ItemSystemProp)).toBe(true);
    expect(set?.byName(`genes_src/${TRANSCRIPT}.map`, key)).toBe(true);
    expect(g.t.ctx.issues.current.value.map((i) => i.code)).toEqual([]);
    g.t.ctx.dispose();
  });

  it("a protein mapped onto it shows its domains along the transcript's exons on chr2", async () => {
    const g = await opened();
    const { f, v } = featuresOf(g.t, 'tx.ranges-2');
    const chr2 = v.resolved.blocks.unit.indexOf(v.set.units.key.indexOf(v.set.strings.find('hg38/chr2')));
    const frags: [number, number][] = [];
    for (let q = 0; q < f.count; q++) if (f.block[q] === chr2) frags.push([f.start[q], f.stop[q]]);
    frags.sort((a, b) => a[0] - b[0]);
    // Alpha (20–200) in the first exon; Beta (480–620) across the first intron (two parts); Gamma (880–1000) across the
    // second.
    expect(frags).toEqual([...exonRuns(20, 200), ...exonRuns(480, 620), ...exonRuns(880, 1000)]);
    g.t.ctx.dispose();
  });

  it("a section of the transcript's own residues: the protein's domains, the genome's variants through its map", async () => {
    const g = await opened({ item: true, section: 'tx.item' });
    const v = g.t.ctx.graph.assembled.value?.view;
    const r = v?.resolved;
    const s = r ? r.sections.ref.indexOf(v.set.strings.find('tx.item')) : -1;
    expect(r?.sections.blockCount[s]).toBe(1);
    const unit = r?.blocks.unit[r.sections.blockStart[s]] ?? -1;
    expect(v?.set.strings.get(v.set.units.key[unit])).toBe(`genes_src/${TRANSCRIPT}/unit`);
    // The domains on the item's residues as on P1's (base for base).
    const { f } = featuresOf(g.t, 'tx.ranges-3');
    const own: [number, number][] = [];
    for (let q = 0; q < f.count; q++) own.push([f.start[q], f.stop[q]]);
    expect(own.sort((a, b) => a[0] - b[0])).toEqual([
      [20, 200],
      [480, 620],
      [880, 1000],
    ]);
    // Variants through the item's codon map, on the residues their positions code.
    const st = g.m.get('tx.variants-2');
    const out = new Int32Array(7);
    const got = new Set<number>();
    for (let k = 0; st?.items; k++) {
      const x = st.items(0, 0, 1334, k, out);
      if (x === null) break;
      for (let q = 0; x.features && q < x.features.count; q++) got.add(x.features.start[q]);
    }
    const want = new Set<number>();
    let base = 0;
    for (const [a, e] of CODING_CDS) {
      for (let x = a; x < e; x++) if (variantAt(1, x)) want.add(Math.floor((base + x - a) / 3));
      base += e - a;
    }
    expect([...got].sort((a, b) => a - b)).toEqual([...want].sort((a, b) => a - b));
    expect(g.t.ctx.issues.current.value.map((i) => i.code)).toEqual([]);
    g.t.ctx.dispose();
  });

  it('hovered, a transcript residue is numbered in its own positions; its codon is read through its map', async () => {
    const g = await opened({ item: true, section: 'tx.item' });
    g.t.vp?.zoomTo(0, 480, 520);
    await g.drain();
    // The label of residue 500 (its codon split by the first intron), read where the pointer would be.
    const h = Hit.create();
    const l = g.t.layout('tx.item');
    let label: Record<string, string> = {};
    for (let y = 1; y < l.height && !label.codon; y += 2)
      for (let x = l.dataLeft; x < l.width && !label.codon; x++) {
        g.t.ctx.layout.hitTest('tx.item', x, y, h);
        if (h.block !== 0 || h.index !== 500) continue;
        g.t.ctx.hover.setFromHit(h, 'pointer');
        const lines = g.t.ctx.labels.compute(g.t.ctx.hover.live.current)?.lines ?? [];
        if (lines.some((q) => q.key === 'codon' && q.value !== 'loading…'))
          label = Object.fromEntries(lines.map((q) => [q.key, q.value]));
        else await g.drain();
      }
    expect(g.t.ctx.labels.compute(g.t.ctx.hover.live.current)?.title).toBe('T1 501');
    expect([label.codon, label['codon at']]).toEqual(['T|CA', 'chr2 100011501 | 100090001–100090002 (+)']);
    g.t.ctx.dispose();
  });

  it("comes and goes with its regions, its id kept; the document's mapping onto it and P1's domains come back", async () => {
    const g = await opened();
    const T1 = `genes_src/${TRANSCRIPT}`;
    const id = idOf(g.t, T1);
    // Far from the gene, under a small cache: its regions evicted, the item goes (nothing shows it), with p1.t1.
    g.t.ctx.config.set(CacheBytes, 1);
    g.t.vp?.zoomTo(0, 10_000_000, 10_100_000);
    await g.drain();
    expect([idOf(g.t, T1), idOf(g.t, 'p1.t1')]).toEqual([-1, -1]);
    g.t.ctx.config.set(CacheBytes, 64 << 20);
    g.t.vp?.zoomTo(1, CODING_CDS[0][0] - 5000, CODING_CDS[2][1] + 5000);
    await g.drain();
    expect([idOf(g.t, T1), idOf(g.t, 'p1.t1') >= 0]).toEqual([id, true]);
    expect(domainsOnChr2(g.t)).toEqual([...exonRuns(20, 200), ...exonRuns(480, 620), ...exonRuns(880, 1000)]);
    expect(codes(g.t)).toEqual([]);
    g.t.ctx.dispose();
  });

  it("made only when something names the item: unnamed, the row's map stays a map", async () => {
    const g = await opened({ named: false });
    expect([idOf(g.t, `genes_src/${TRANSCRIPT}`), idOf(g.t, `genes_src/${TRANSCRIPT}.map`)]).toEqual([-1, -1]);
    // Held in the manager (counted in its cache), not read into the set.
    expect(g.m.nodes.entries.has(`genes_src/${TRANSCRIPT}`)).toBe(true);
    g.t.ctx.dispose();
  });

  it('a storm of items with maps nothing names reruns no section, layer or view, and the set is the same', async () => {
    const g = await opened({ named: false });
    const set0 = g.t.ctx.graph.assembled.value?.set;
    const reran = new Set<string>();
    for (let k = 0; k < 300; k++) {
      const out: StreamNode[] = [];
      const t = { id: [`X${k}`], map: [{ method: 'codon', cds: [[100 * k, 100 * k + 300]] }], label: [`X${k}`] };
      itemSystemNodes('genes_src', 'hg38/chr2', t, out, () => {});
      g.m.nodes.add('genes_src', out, true);
      await g.t.ctx.graph.update();
      for (const r of g.t.ctx.graph.ran) reran.add(r);
    }
    expect([...reran].sort()).toEqual(['#assemble.stream']);
    expect(g.t.ctx.graph.assembled.value?.set).toBe(set0);
    g.t.ctx.dispose();
  });

  it("a map's faults: once per source as warnings; no item coordinate system for them", async () => {
    for (const [map, cds, said] of [
      [{ method: 'codon' }, undefined, 'a codon map needs cds'],
      [
        null,
        [
          [100_010_000, 100_011_500],
          [100_090_000, 100_091_201],
        ],
        'not a multiple of 3',
      ],
      [{ method: 'segments', segments: [[0, 100_010_000, null]] }, undefined, "each segment's length"],
      // A reversed CDS range, as core's codon map refuses it (VS-27).
      [null, [[100_011_500, 100_010_000]], 'is empty or reversed'],
    ] as const) {
      const g = await opened({
        map: map as Record<string, unknown> | null,
        cds: cds as [number, number][] | undefined,
      });
      const issues = g.t.ctx.issues.current.value.filter((i) => i.code === 'svs.stream.item_map');
      expect([issues.length, issues[0]?.severity, issues[0]?.message.includes(said)]).toEqual([1, 'warning', true]);
      expect(idOf(g.t, `genes_src/${TRANSCRIPT}`)).toBe(-1);
      g.t.ctx.dispose();
    }
  });

  it('from_range gives the length; a segments map; the minus strand', async () => {
    const r = await opened({ map: { method: 'codon', cds: CODING_CDS, from_range: [0, 1334] } });
    expect(lengthOf(r.t, `genes_src/${TRANSCRIPT}`)).toBe(1334);
    r.t.ctx.dispose();
    // Segments: the transcript's first 100 bases, base for base: P1's Alpha (20–200) over it as its first 100.
    const s = await opened({ map: { method: 'segments', segments: [[0, CODING_CDS[0][0], 100]] } });
    expect(lengthOf(s.t, `genes_src/${TRANSCRIPT}`)).toBe(100);
    expect(domainsOnChr2(s.t)).toEqual([[CODING_CDS[0][0] + 20, CODING_CDS[0][0] + 100]]);
    s.t.ctx.dispose();
    // Minus: the same exons read downward: Alpha at the top of the last exon.
    const minus = await opened({ map: { method: 'codon', cds: CODING_CDS.slice().reverse(), strand: '-' } });
    const top = CODING_CDS[2][1];
    expect(domainsOnChr2(minus.t)[0 + 2]).toEqual([top - 600, top - 60]);
    minus.t.ctx.dispose();
  });

  it("a mapping naming the item's unit, and an annotation targeting the item, wait for it and join", async () => {
    const u = await opened({ to: `genes_src/${TRANSCRIPT}/unit` });
    expect([idOf(u.t, 'p1.t1') >= 0, domainsOnChr2(u.t).length, codes(u.t)]).toEqual([true, 5, []]);
    u.t.ctx.dispose();
    const a = await opened({ target: true });
    expect([domainsOnChr2(a.t).length, codes(a.t)]).toEqual([5, []]);
    a.t.ctx.dispose();
  });

  it('a mapping retargeted between the item and a document sequence, both ways, stays in the set', async () => {
    const g = await opened();
    const retarget = async (to: string) => {
      await g.t.ctx.applyPatch({
        format: 'seq_view_spec_patch',
        ops: [{ op: 'set', ref: 'p1.t1', params: { to } }],
      } as never);
      await g.drain();
      return idOf(g.t, 'p1.t1') >= 0;
    };
    expect(await retarget('p1copy')).toBe(true);
    expect(await retarget(`genes_src/${TRANSCRIPT}`)).toBe(true);
    expect(domainsOnChr2(g.t).length).toBe(5);
    g.t.ctx.dispose();
  });

  it('shown, the item keeps its map when the regions that carried them are evicted: its section keeps its tracks', async () => {
    const g = await opened({ item: true, section: 'tx.genome' });
    g.t.ctx.config.set(CacheBytes, 1);
    g.t.vp?.zoomTo(0, 10_000_000, 10_100_000);
    await g.drain();
    expect([idOf(g.t, `genes_src/${TRANSCRIPT}`) >= 0, idOf(g.t, `genes_src/${TRANSCRIPT}.map`) >= 0]).toEqual([
      true,
      true,
    ]);
    const v = g.t.ctx.graph.assembled.value?.view;
    const r = v?.resolved;
    let k = 0;
    while (r && k < r.layers.count && v?.set.strings.get(r.layers.ref[k]) !== 'tx.variants-2') k++;
    expect(r?.layers.placementCount[k]).toBeGreaterThan(0);
    g.t.ctx.dispose();
  });

  it("a provider's manifest may not add an item coordinate system", async () => {
    const node = { kind: ItemSystemKind as StreamNode['kind'], key: 'Z', params: { length: 9 } };
    const g = await opened({ nodes: [node] });
    expect([codes(g.t).includes('svs.node.kind'), idOf(g.t, 'phylop_src/Z')]).toEqual([true, -1]);
    g.t.ctx.dispose();
  });

  it("templates: item_system lists the source's item coordinate systems; sequence does not", async () => {
    const g = await opened({ templates: true });
    const v = g.t.ctx.graph.assembled.value?.view;
    const tracks = v ? Array.from(v.resolved.tracks.ref, (x) => v.set.strings.get(x)) : [];
    const items = tracks.filter((x) => x.includes(`genes_src/${TRANSCRIPT}`));
    expect(items.length).toBe(1);
    g.t.ctx.dispose();
  });

  it("an item's regions released one by one: it stays until the last; a conflicting map is said; sources apart", () => {
    const said: string[] = [];
    const nodes = new DynamicNodes(
      () => {},
      (code) => said.push(code),
      heldNodeRefs(createSeqViewer(DefaultSeqViewerSpec()).registry.svs),
    );
    const out: StreamNode[] = [];
    const t = { id: ['T1'], map: [{ method: 'codon', cds: [[0, 30]] }], label: ['T1'] };
    itemSystemNodes('genes_src', 'hg38/chr2', t, out, () => {});
    const a = nodes.add('genes_src', out, true) ?? [];
    const b = nodes.add('genes_src', out, true) ?? [];
    nodes.release(a);
    expect(nodes.entries.has('genes_src/T1')).toBe(true);
    nodes.release(b);
    expect(nodes.entries.has('genes_src/T1')).toBe(false);
    const other: StreamNode[] = [];
    itemSystemNodes('other_src', 'hg38/chr2', t, other, () => {});
    nodes.add('genes_src', out, true);
    nodes.add('other_src', other, true);
    expect([nodes.entries.has('genes_src/T1'), nodes.entries.has('other_src/T1')]).toEqual([true, true]);
    const changed: StreamNode[] = [];
    itemSystemNodes('genes_src', 'hg38/chr2', { ...t, map: [{ method: 'codon', cds: [[0, 33]] }] }, changed, () => {});
    nodes.add('genes_src', changed, true);
    expect(said).toContain('svs.stream.conflict');
  });

  it("a selection of the transcript's residues lights their codons on chr2, across its map", async () => {
    const g = await opened({ item: true });
    const set = g.t.ctx.graph.assembled.value?.set as SeqSet;
    const unit = Array.from(set.units.key).indexOf(set.strings.find(`genes_src/${TRANSCRIPT}/unit`));
    const loci = {
      kind: 'seq-loci',
      set,
      elements: [{ unit, indices: OrderedSet.ofRange(499, 501) }],
    } as unknown as SeqLoci;
    g.t.ctx.selection.set(loci, 'api');
    const ch = new ChannelMarks(g.t.ctx);
    ch.sync();
    const v = g.t.ctx.graph.assembled.value?.view;
    const r = v?.resolved;
    const sec = r ? r.sections.ref.indexOf(v.set.strings.find('tx.genome')) : -1;
    const chr2 = r ? r.blocks.unit.indexOf(unitOf(set, 'hg38/chr2')) - (r.sections.blockStart[sec] ?? 0) : -1;
    const m = ch.marks(sec)?.runs(MarkBit.selection);
    const runs: number[][] = [];
    for (let q = m ? m.from(chr2) : 0; m && q < m.to(chr2); q++) runs.push([m.start[q], m.stop[q]]);
    // Residues 499 and 500: 500's codon split by the first intron.
    expect(runs).toEqual(exonRuns(499, 501));
    g.t.ctx.dispose();
  });

  it('an annotation whose rows target the item one by one is said not to be supported yet', async () => {
    const g = await opened({ rows: true });
    const said = g.t.ctx.issues.current.value.find((i) => i.message.includes('per-row targets'));
    // An unsupported case, not an unknown ref: the node is unavailable, a warning (GV-10).
    expect([said?.code, said?.severity]).toEqual(['svs.provider.unavailable', 'warning']);
    g.t.ctx.dispose();
  });

  it("a held-back mapping's issue stays while it is in, across later stream versions", async () => {
    const g = await opened();
    await g.t.ctx.applyPatch({
      format: 'seq_view_spec_patch',
      ops: [{ op: 'set', ref: 'p1.t1', params: { segments: [[0, 2000, 5]] } }],
    } as never);
    await g.drain();
    const has = () => g.t.ctx.issues.current.value.some((i) => i.message.includes('"p1.t1"'));
    expect(has()).toBe(true);
    // More arrivals (a pan): new stream versions, the issue kept.
    g.t.vp?.zoomTo(1, CODING_CDS[0][0] - 50_000, CODING_CDS[2][1] + 50_000);
    await g.drain();
    expect(has()).toBe(true);
    g.t.ctx.dispose();
  });

  it("a held node names what core's `nodeRefs` reads: a dictionary-coded table target holds its node (VS-24)", () => {
    const nodes = new DynamicNodes(
      () => {},
      () => {},
      heldNodeRefs(createSeqViewer(DefaultSeqViewerSpec()).registry.svs),
    );
    const seq: StreamNode = { kind: 'sequence', key: 's1', params: { value: 'MKV' } };
    nodes.add('src', [seq]);
    const ann: StreamNode = {
      kind: 'annotation',
      key: 'a1',
      params: {
        geometry: 'ranges',
        value: { target: { dict: ['src/s1'], codes: [0, 0] }, start: [0, 1], stop: [1, 2] },
      },
    };
    nodes.add('src', [ann]);
    expect(nodes.entries.get('src/s1')?.referrers).toBe(1);
    expect(nodes.entries.get('src/a1')?.names).toEqual(['src/s1']);
  });

  it("a provider's patch may not add an item coordinate system itself", () => {
    const said: string[] = [];
    const nodes = new DynamicNodes(
      () => {},
      (code) => said.push(code),
      heldNodeRefs(createSeqViewer(DefaultSeqViewerSpec()).registry.svs),
    );
    const node = { kind: ItemSystemKind as StreamNode['kind'], key: 'X', params: { length: 5 } };
    expect([nodes.add('genes_src', [node]), said, nodes.entries.size]).toEqual([null, ['svs.node.kind'], 0]);
    expect([nodes.add('genes_src', [node], true), nodes.entries.size]).toEqual([['genes_src/X'], 1]);
  });
});
