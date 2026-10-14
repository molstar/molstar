/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST7: streamed data through a path (seqspec:sources#requests; seqspec:model#paths). A protein's block shows its
// genome's streamed genes, variants and phyloP through its codon map: the regions asked are those its residues' exons
// reach (three exons over 220 kb: three regions of items, never the introns between them), read back on the block's
// indices (a variant on its residue, a gene's exon clipped to its coding part and labeled, values over a residue's
// codon), and a streamed item hovered as its range on the protein. ST7's review: an item starting between two exons'
// regions, a map edited under held items, the minus strand, a source with levels finer than a codon, the hover at a
// codon split across two regions, case 9's coverage on its chains through their rows.

import { createSeqViewSpec } from "@molstar/seq-core/builder";
import type { DataPatch, Loaders, SeqViewSpecDoc, StreamRequest } from "@molstar/seq-core/format";
import { StreamNotServed } from "@molstar/seq-core/format";
import { Paths } from "@molstar/seq-core/model";
import { tilesStream } from "@molstar/seq-core/stream/tiles";
import { tileServer } from "@molstar/seq-core/stream/tiles/server";
import { describe, expect, it } from "vitest";
import { case9 } from "../../core/_test/format/cases.ts";
import { fixtureLoaders } from "../../core/_test/format/load-fixtures.ts";
import {
  CODING_CDS,
  codingDoc,
  codingProtein,
  GENE_LEVELS,
  GRCH38,
  hg38Streams,
  itemsManifest,
  letterAt,
  type Scripted,
  sampledItems,
  scripted,
  VOCABULARY,
  variantAt,
} from "../../core/_test/format/stream-provider.ts";
import { KernelRenderers, type LayerRenderer, RegionState } from "../src/index.ts";
import type { StreamItems, Streams } from "../src/streams.ts";
import { composed, recording } from "./composer-fixtures.ts";
import { driven } from "./input-fixtures.ts";
import { MaxInFlight } from "./stream-fixtures.ts";

const flush = () => new Promise((ok) => setTimeout(ok, 0));
const THREE = GRCH38.slice(0, 3);
const SECTION = "coding.protein";
const RESIDUES = codingProtein().length;
/** The codons the introns split. */
const SPLIT = [500, 900];

const KEYS = THREE.map((c) => c[0]);
/** Sources in place of the fixture's (by stream id). */
type Over = Partial<Record<"genes" | "variants", Scripted>>;
async function protein(o: { doc?: SeqViewSpecDoc; over?: Over; dpr?: number } = {}) {
  const base = hg38Streams({ chromosomes: THREE, gated: true });
  const s = { ...base, ...o.over };
  const stream = (src: { stream: string | null }) =>
    src.stream && src.stream in (o.over ?? {}) ? (o.over?.[src.stream as "genes"]?.provider ?? null) : base.stream(src);
  const loaders = { load: () => Promise.reject(new Error("none")), stream };
  const doc = o.doc ?? codingDoc();
  const t = await composed(doc, { width: 1000, height: 300, dpr: o.dpr ?? 1, loaders, section: SECTION });
  t.ctx.config.set(MaxInFlight, 256);
  const m = t.ctx.streams as Streams;
  const ref = (type: string) => {
    const v = t.ctx.graph.assembled.value?.view;
    for (let k = 0; v && k < v.resolved.layers.count; k++)
      if (v.set.strings.get(v.resolved.layers.type[k]) === type) return v.set.strings.get(v.resolved.layers.ref[k]);
    throw new Error(`no ${type} layer`);
  };
  const show = (start: number, stop: number) => {
    t.vp?.zoomTo(0, start, stop);
    t.layout(SECTION);
    m.pass();
  };
  const drain = async () => {
    for (let k = 0; k < 6; k++) {
      for (const p of [s.genome, s.phylop, s.genes, s.variants]) p.release();
      await flush();
      m.pass();
    }
  };
  return { t, s, m, ref, show, drain };
}

/** Whether [a, e) of chr2 overlaps the CDS. */
const coding = (a: number, e: number) => CODING_CDS.some(([s, t]) => a < t && s < e);

type J = { kind: string; ref?: string; params?: Record<string, unknown>; children?: J[] };
const byRef = (n: J, ref: string): J | undefined =>
  n.ref === ref ? n : n.children?.map((c) => byRef(c, ref)).find(Boolean);
/** Case 9 with a streamed coverage over its columns (`seqstar-tiles/1`, value i / 320), a layer of it on the chains. */
function coverageCase9(): { doc: SeqViewSpecDoc; loaders: Loaders } {
  const doc = case9();
  const COV = "https://example.org/pfam/PF00042/coverage";
  const data = (doc.root as J).children?.find((c) => c.kind === "data");
  data?.children?.push(
    { kind: "source", ref: "cov.src", params: { uri: COV, protocol: "seqstar-tiles/1" } },
    {
      kind: "annotation",
      ref: "cov",
      params: { geometry: "values", source_ref: "cov.src", target: "pf00042", length: 320 },
    },
  );
  const layer = { kind: "layer", ref: "cov.a", params: { type: "values", data: "cov" } };
  byRef(doc.root as J, "msa.chains")?.children?.push({ kind: "track", ref: "cov.chains", children: [layer] });
  const cov = tileServer({ rowCount: 1, columnCount: 320, length: 320, values: (i) => i / 320 });
  const base = fixtureLoaders();
  const loaders: Loaders = {
    load: base.load,
    stream: (source) =>
      source.uri === COV ? tilesStream(() => cov.transport)(source) : (base.stream?.(source) ?? null),
  };
  return { doc, loaders };
}

/** The residue position x of chr2 codes, on either strand (the CDS read in transcript order). */
function residueOf(x: number, strand: "+" | "-" = "+"): number {
  const cds = strand === "+" ? CODING_CDS : CODING_CDS.slice().reverse();
  let base = 0;
  for (const [a, e] of cds) {
    if (x >= a && x < e) return Math.floor((base + (strand === "+" ? x - a : e - 1 - x)) / 3);
    base += e - a;
  }
  return -1;
}
/** The residues the variants layer's regions hold features on (block 0), and the expected ones on `strand`. */
function variantResidues(g: Awaited<ReturnType<typeof protein>>, strand: "+" | "-" = "+") {
  const st = g.m.get(g.ref("variants"));
  const out = new Int32Array(7);
  const got = new Set<number>();
  for (let k = 0; st?.items; k++) {
    const r = st.items(0, 0, RESIDUES, k, out);
    if (r === null) break;
    for (let q = 0; r.features && q < r.features.count; q++) got.add(r.features.start[q]);
  }
  const want = new Set<number>();
  for (const [a, e] of CODING_CDS) for (let x = a; x < e; x++) if (variantAt(1, x)) want.add(residueOf(x, strand));
  return { got: [...got].sort((a, b) => a - b), want: [...want].sort((a, b) => a - b) };
}
/** Per item id the regions that draw a fragment of it (the drawn-once rule: `first` at or past `out[3]`). */
function drawnBy(st: ReturnType<Streams["get"]>, out: Int32Array) {
  const by = new Map<string, number[]>();
  for (let k = 0; st?.items; k++) {
    const r = st.items(0, 0, RESIDUES, k, out);
    if (r === null) break;
    const f = r.features;
    const ids = r.table?.id;
    for (let q = 0; f && ids && q < f.count; q++) {
      if (r.first[q] < out[3]) continue;
      const id = String(ids[f.row[q]]);
      by.set(id, [...new Set([...(by.get(id) ?? []), k])]);
    }
  }
  return by;
}
describe("streamed data through a path", () => {
  it("places the genome's layers on the protein's block, with no issue", async () => {
    const { t } = await protein();
    const v = t.ctx.graph.assembled.value?.view;
    expect(t.ctx.issues.current.value.map((i) => i.code)).toEqual([]);
    expect(v?.resolved.layers.count).toBe(5);
    const paths = v?.resolved.paths ?? [];
    // The path through chr2 holds the CDS (4,002 positions), not the chromosome's 242 Mb.
    expect([paths.length > 0, Math.max(...paths.map((p) => p.heldTo.length))]).toEqual([true, 4002]);
  });

  it("asks only the regions the residues' exons reach, on chr2, each once", async () => {
    const g = await protein();
    g.show(0, RESIDUES);
    await g.drain();
    for (const p of [g.s.genes, g.s.variants, g.s.phylop]) {
      expect(p.requests.length).toBeGreaterThan(0);
      expect(p.requests.every((r) => r.system === "chr2" && coding(r.start, r.stop))).toBe(true);
    }
    // Items: one data region of 65,536 an exon (none for the 200 kb of introns between them).
    expect(g.s.genes.requests.map((r) => [r.level, r.start / 65_536]).sort((a, b) => a[1] - b[1])).toEqual(
      CODING_CDS.map(([a]) => [-1, Math.floor(a / 65_536)]),
    );
    // Zoomed to residues 600–700 (the second exon): nothing more of the first or the third.
    const before = g.s.variants.requests.length;
    g.show(600, 700);
    await g.drain();
    expect(g.s.variants.requests.slice(before).every((r) => coding(r.start, r.stop))).toBe(true);
  });

  it("reads back on the block's indices: variants on their residues, exons clipped, values over codons", async () => {
    const g = await protein();
    g.show(0, RESIDUES);
    await g.drain();
    const st = g.m.get(g.ref("variants"));
    const out = new Int32Array(7);
    const residues: number[] = [];
    for (let k = 0; st?.items; k++) {
      const r = st.items(0, 0, RESIDUES, k, out);
      if (r === null) break;
      // A region between two exons' (never asked) reaches no residue: an empty extent, pending.
      const intronic = out[1] === out[2];
      expect([intronic || out[4] !== RegionState.pending, out[1] >= 0, out[2] <= RESIDUES]).toEqual([true, true, true]);
      const f = r.features;
      for (let q = 0; f && q < f.count; q++) residues.push(f.start[q]);
    }
    // Every variant of the CDS, on the residue its position codes.
    const expected: number[] = [];
    let base = 0;
    for (const [a, e] of CODING_CDS) {
      for (let x = a; x < e; x++) if (variantAt(1, x)) expected.push(Math.floor((base + x - a) / 3));
      base += e - a;
    }
    expect(residues.length).toBeGreaterThan(0);
    expect([...new Set(residues)].sort((a, b) => a - b)).toEqual([...new Set(expected)]);
    // A gene's exons land as their coding parts, inside the protein.
    const genes = g.m.get(g.ref("ranges"));
    for (let k = 0; genes?.items; k++) {
      const r = genes.items(0, 0, RESIDUES, k, out);
      if (r === null) break;
      const f = r.features;
      for (let q = 0; f && q < f.count; q++) expect(f.start[q] >= 0 && f.stop[q] <= RESIDUES).toBe(true);
    }
    // Values over a residue's codon, its three positions at the data level, both parts of a codon an intron splits.
    const values = g.m.get(g.ref("values"));
    const v = new Float64Array(4);
    for (const i of [0, SPLIT[0] - 1, SPLIT[0], SPLIT[1], RESIDUES - 1]) {
      v[0] = 0;
      expect([i, values?.values?.(0, i, i + 1, v), v[3], Number.isFinite(v[0])]).toEqual([i, true, 3, true]);
    }
    expect(values?.stateAt?.(0, 10)).toBe(RegionState.loaded);
  });
  it("draws the genome's tracks on the protein; an item's range on the protein is the hull of its parts", async () => {
    const g = await protein();
    g.show(0, RESIDUES);
    await g.drain();
    const frame = g.t.frame(recording());
    expect(frame.opCount).toBeGreaterThan(0);
    // A variant's range read back on the block: its residue.
    const st = g.m.get(g.ref("variants"));
    const out = new Int32Array(7);
    let r: StreamItems | null | undefined = null;
    for (let k = 0; !r?.features?.count; k++) r = st?.items?.(0, 0, RESIDUES, k, out);
    const f = r?.features;
    const span = new Int32Array(2);
    const end = Math.max(r?.ends?.[0] ?? 0, (r?.first[0] ?? 0) + 1);
    expect(st?.toBlock?.(0, r?.first[0] ?? 0, end, span)).toBe(true);
    expect([span[0], span[1]]).toEqual([f?.start[0], f?.stop[0]]);
  });

  it("labels an item on its first fragment on the protein", async () => {
    const g = await protein();
    g.show(0, RESIDUES);
    await g.drain();
    const st = g.m.get(g.ref("ranges"));
    const out = new Int32Array(7);
    const labels: string[] = [];
    for (let k = 0; st?.items; k++) {
      const r = st.items(0, 0, RESIDUES, k, out);
      if (r === null) break;
      for (const x of r.features?.label ?? []) if (x) labels.push(x);
    }
    expect(labels).toEqual(["G2705"]);
  });

  it("draws an item starting between two exons' regions once, from the region it reaches", async () => {
    // A gene from 100,150,000 (a region no exon reaches, never asked) into the third exon.
    const genes = scripted(itemsManifest(THREE, GENE_LEVELS, VOCABULARY.genes), (r: StreamRequest): DataPatch => {
      const p = sampledItems(r, KEYS, "genes", GENE_LEVELS);
      if (p.kind !== "items" || r.start > 100_230_100 || r.stop <= 100_150_000) return p;
      const t = p.items as Record<string, unknown[]>;
      for (const [k, v] of Object.entries({ start: 100_150_000, stop: 100_230_100, id: "bridge", label: "B" }))
        t[k] = [...(t[k] ?? []), v];
      t.category = [...(t.category ?? []), "lncRNA"];
      return p;
    });
    genes.gated = true;
    const g = await protein({ over: { genes } });
    g.show(0, RESIDUES);
    await g.drain();
    const by = drawnBy(g.m.get(g.ref("ranges")), new Int32Array(7));
    expect(by.get("bridge")?.length).toBe(1);
  });

  it("places held items again through an edited map: the strand flipped on the same exons", async () => {
    const g = await protein();
    g.show(0, RESIDUES);
    await g.drain();
    const plus = variantResidues(g);
    expect(plus.got).toEqual(plus.want);
    const cds = CODING_CDS.slice().reverse();
    await g.t.ctx.applyPatch({
      format: "seq_view_spec_patch",
      ops: [{ op: "set", ref: "p1.codon", params: { cds, strand: "-" } }],
    } as never);
    g.show(0, RESIDUES);
    await g.drain();
    const minus = variantResidues(g, "-");
    expect([minus.got, minus.got.length > 0]).toEqual([minus.want, true]);
  });

  it("reads the minus strand: variants on the residues their positions code there", async () => {
    const g = await protein({ doc: codingDoc({ strand: "-" }) });
    g.show(0, RESIDUES);
    await g.drain();
    const r = variantResidues(g, "-");
    expect([r.got, r.got.length > 0]).toEqual([r.want, true]);
  });

  it("a source whose levels are finer than a codon: the draw reads the level the pass asked", async () => {
    // Levels of 1, 2 and 4 positions a bin: 1,334 residues over 2,000 device px is 2 positions a device px.
    const variants = scripted(itemsManifest(THREE, [1, 2, 4], VOCABULARY.variants), (r) =>
      sampledItems(r, KEYS, "variants", [1, 2, 4]),
    );
    variants.gated = true;
    const g = await protein({ over: { variants }, dpr: 2 });
    g.t.vp?.fit();
    g.t.layout(SECTION);
    g.m.pass();
    await g.drain();
    const st = g.m.get(g.ref("variants"));
    const v = g.t.ctx.layout.axis(SECTION);
    const out = new Int32Array(7);
    out[0] = Math.floor(1 / ((v?.cell ?? 1) * 2));
    const states: number[] = [];
    for (let k = 0; st?.items; k++) {
      if (st.items(0, 0, RESIDUES, k, out) === null) break;
      states.push(out[4]);
    }
    expect(states.length).toBeGreaterThan(0);
    expect(states.every((x) => x !== RegionState.pending)).toBe(true);
  });

  it("hovers a gene at a codon split across two regions, as its range on the protein", async () => {
    const s = hg38Streams({ chromosomes: THREE, gated: true });
    const t = await driven(codingDoc(), {
      renderers: [...KernelRenderers] as LayerRenderer[],
      loaders: { load: () => Promise.reject(new Error("none")), stream: s.stream },
      size: { [SECTION]: [1000, 300] },
    });
    t.vp(SECTION).zoomTo(0, 450, 550);
    await t.ctx.graph.update();
    for (let k = 0; k < 4; k++) {
      t.layout(SECTION);
      (t.ctx.streams as Streams).pass();
      for (const p of [s.genome, s.genes, s.variants, s.phylop]) p.release();
      await flush();
    }
    const l = t.layout(SECTION);
    const m = t.ctx.graph.assembled.value?.view;
    const h = t.ctx.hover.live.current;
    let hit = false;
    // Residue 500's codon: one base in the first exon's region, two in the second's (where G2705 is).
    for (let x = 0; x < 1000 && !hit; x += 2)
      for (let y = 1; y < 300 && !hit; y += 2) {
        t.pointer("move", SECTION, l.dataLeft + x, y);
        const genes = h.track >= 0 && m?.set.strings.get(m.resolved.tracks.label[h.track]) === "Genes";
        hit = genes && h.index === 500 && h.part >= 0;
      }
    expect([hit, h.span[0], h.span[1]]).toEqual([true, 500, 901]);
    t.ctx.dispose();
  });
});

describe("the codon label through a path", () => {
  /** A section in the input harness at block 0's [a, e), driven by passes; the label at an index of it. */
  async function opened(
    doc: SeqViewSpecDoc,
    stream: Loaders["stream"],
    held: Scripted[],
    section: string,
    a = 0,
    e = 0,
  ) {
    const t = await driven(doc, {
      renderers: [...KernelRenderers] as LayerRenderer[],
      loaders: { load: () => Promise.reject(new Error("none")), stream },
      size: { [section]: [1000, 300] },
    });
    if (e > a) t.vp(section).zoomTo(0, a, e);
    await t.ctx.graph.update();
    const settle = async () => {
      for (let k = 0; k < 4; k++) {
        t.layout(section);
        (t.ctx.streams as Streams).pass();
        for (const p of held) p.release();
        await flush();
      }
    };
    await settle();
    /** The label at index i of block 0 (the pointer over it, on the first row that hits the block). */
    const at = (i: number) => {
      const l = t.layout(section);
      const h = t.ctx.hover.live.current;
      let y = 1;
      for (; y < 300; y += 2) {
        t.pointer("move", section, l.dataLeft + 2, y);
        if (h.block === 0 && h.index >= 0) break;
      }
      for (let x = 0; x < 1000 && h.index !== i; x++) t.pointer("move", section, l.dataLeft + x, y);
      const lines = t.ctx.labels.compute(h)?.lines ?? [];
      return Object.fromEntries(lines.map((x) => [x.key, x.value]));
    };
    return { t, settle, at };
  }
  /** The protein's section at residues [a, e), its letters held or not; `genome` in place of the letters stream. */
  async function hovered(strand: "+" | "-", a: number, e: number, genome?: Scripted["provider"]) {
    const s = hg38Streams({ chromosomes: THREE, gated: true, letters: true });
    const stream = (src: { stream: string | null }) => (src.stream === "genome" && genome ? genome : s.stream(src));
    const g = await opened(
      codingDoc({ strand }),
      stream as Loaders["stream"],
      [s.genome, s.genes, s.variants, s.phylop],
      SECTION,
      a,
      e,
    );
    return { ...g, s };
  }
  /** Residue i's three positions in reading order, and its codon read there. */
  function codonOf(i: number, strand: "+" | "-") {
    const xs: number[] = [];
    for (const [a, e] of strand === "+" ? CODING_CDS : CODING_CDS.slice().reverse())
      if (strand === "+") for (let x = a; x < e; x++) xs.push(x);
      else for (let x = e - 1; x >= a; x--) xs.push(x);
    const at = xs.slice(3 * i, 3 * i + 3);
    const C: Record<string, string> = { A: "T", T: "A", C: "G", G: "C", N: "N" };
    const bases = at.map((x) => letterAt(1, x).toUpperCase()).map((c) => (strand === "+" ? c : C[c]));
    return { at, bases };
  }

  it("shows a residue's codon and its positions; a codon an intron splits in two parts", async () => {
    const g = await hovered("+", 480, 520);
    // The first hover asks for the letters; they arrive, and the label is made again.
    g.at(490);
    g.at(500);
    await g.settle();
    const one = codonOf(490, "+");
    const label = g.at(490);
    expect([label.codon, label["codon at"]]).toEqual([
      one.bases.join(""),
      `chr2 ${one.at[0] + 1}–${one.at[2] + 1} (+)`,
    ]);
    const split = codonOf(500, "+");
    const two = g.at(500);
    expect(two.codon).toBe(`${split.bases[0]}|${split.bases[1]}${split.bases[2]}`);
    expect(two["codon at"]).toBe(`chr2 ${split.at[0] + 1} | ${split.at[1] + 1}–${split.at[2] + 1} (+)`);
    g.t.ctx.dispose();
  });

  it("on the minus strand: the complement, read downward", async () => {
    const g = await hovered("-", 0, 40);
    g.at(20);
    await g.settle();
    const c = codonOf(20, "-");
    const label = g.at(20);
    expect([label.codon, label["codon at"]]).toEqual([c.bases.join(""), `chr2 ${c.at[0] + 1}–${c.at[2] + 1} (−)`]);
    g.t.ctx.dispose();
  });

  it("says it is loading until the letters arrive, and asks for them", async () => {
    const g = await hovered("+", 480, 520);
    const before = g.s.genome.requests.length;
    // Residue 900's codon lies in the third exon: no letters of it held yet.
    g.t.vp(SECTION).zoomTo(0, 880, 920);
    const first = g.at(905);
    expect(first.codon).toBe("loading…");
    await g.settle();
    expect(g.s.genome.requests.length).toBeGreaterThan(before);
    expect(g.at(905).codon).toBe(codonOf(905, "+").bases.join(""));
    g.t.ctx.dispose();
  });

  it("a minus-strand codon an intron splits: the complement in reading order, positions high to low", async () => {
    // On the minus strand the last exon (1,301 bases) is read first: residue 433's codon is split.
    const g = await hovered("-", 420, 450);
    g.at(433);
    await g.settle();
    const c = codonOf(433, "-");
    const label = g.at(433);
    expect(label.codon).toBe(`${c.bases[0]}${c.bases[1]}|${c.bases[2]}`);
    expect(label["codon at"]).toBe(`chr2 ${c.at[0] + 1}–${c.at[1] + 1} | ${c.at[2] + 1} (−)`);
    g.t.ctx.dispose();
  });

  it("a letters region that fails says so (Retry reaches it); one not served says that", async () => {
    for (const [e, said] of [
      [new Error("503"), "failed: Retry in a track's … menu"],
      [new StreamNotServed("not served here"), "not served"],
    ] as const) {
      const base = hg38Streams({ chromosomes: THREE, letters: true });
      const genome: Scripted["provider"] = {
        fetch: (r, sig) => (r.system === null ? base.genome.provider.fetch(r, sig) : Promise.reject(e)),
      };
      const g = await hovered("+", 480, 520, genome);
      g.at(490);
      await g.settle();
      expect(g.at(490).codon).toBe(said);
      if (said.startsWith("failed")) expect(g.t.ctx.streams.retry("coding.values")).toBeGreaterThan(0);
      g.t.ctx.dispose();
    }
  });

  it("a hover that moves on before its letters arrive: only the last codon's region is asked", async () => {
    const g = await hovered("+", 880, 920);
    const before = g.s.genome.requests.length;
    g.at(905);
    g.t.vp(SECTION).zoomTo(0, 480, 520);
    g.at(490);
    await g.settle();
    const asked = g.s.genome.requests.slice(before).map((r) => r.start);
    expect(asked).toEqual([Math.floor(codonOf(490, "+").at[0] / 4096) * 4096]);
    g.t.ctx.dispose();
  });

  it("another document while letters are pending: the old stream is not asked again, and the view settles", async () => {
    const g = await hovered("+", 880, 920);
    g.at(905);
    const other = createSeqViewSpec({ title: "Another" });
    const x = other.data.sequence({ ref: "x", alphabet: "protein", value: "MKVLAAGIV" });
    other.view({ ref: "ov", label: "Another" }).section({ ref: "ov.main", label: "x" }).block({ sequence: x });
    const opened = await g.t.ctx.setDocument(other.getDocument());
    expect([opened.ok, (g.t.ctx.streams as unknown as { asides: unknown[] }).asides.length]).toEqual([true, 0]);
    const before = g.s.genome.requests.length;
    for (let k = 0; k < 4; k++) {
      (g.t.ctx.streams as Streams).pass();
      g.s.genome.release();
      await flush();
    }
    expect(g.s.genome.requests.length).toBe(before);
    await g.t.ctx.streams.settled();
    g.t.ctx.dispose();
  });

  it("a molecule held whole (letters by value, a cyclic one across its origin), and no codon without a codon map", async () => {
    // A plasmid of 30 bases: a CDS across its origin (27–30, then 0–5), a gene item streamed over it.
    const dna = "ATGAAACCCGGGTTTAAACCCGGGTTTGCA";
    const b = createSeqViewSpec({ title: "Plasmid" });
    const plasmid = b.data.sequence({ ref: "pl", alphabet: "dna", value: dna, cyclic: true });
    const pep = b.data.sequence({ ref: "pep", alphabet: "protein", value: "AMK" });
    b.data.mapping({
      ref: "pep.codon",
      from: pep,
      to: plasmid,
      method: "codon",
      cds: [
        [27, 30],
        [0, 6],
      ],
    });
    b.data.source({ ref: "feat_src", stream: "feat", data_kind: "items", coordinate_system: "pl" });
    b.data.annotation({ ref: "feat", geometry: "ranges", source_ref: "feat_src" });
    const pv = b.view({ ref: "pv", label: "Plasmid" });
    const v = pv.section({ ref: "pv.main", label: "pep" });
    v.block({ sequence: pep });
    v.track({ label: "Features" }).ranges({ data: "feat" });
    v.track({ label: "pep" }).letters();
    const cp = b.data.sequence({ ref: "cp", alphabet: "dna", value: dna, cyclic: true });
    b.data.mapping({ ref: "cp.map", from: cp, to: plasmid, method: "segments", segments: [[0, 30, 0]] });
    const w = pv.section({ ref: "pv.copy", label: "copy" });
    w.block({ sequence: cp });
    w.track({ label: "Features" }).ranges({ data: "feat" });
    w.track({ label: "cp" }).letters();
    const feat = scripted(itemsManifest([["pl", 30]], [], {}), (r) => ({
      source: r.source,
      system: r.system,
      start: r.start,
      stop: r.stop,
      row_start: 0,
      row_stop: 0,
      level: r.level,
      kind: "items",
      items: { start: [0], stop: [30], id: ["f"], label: ["F"] },
    }));
    const g = await opened(
      b.getDocument(),
      ((src: { stream: string | null }) => (src.stream === "feat" ? feat.provider : null)) as Loaders["stream"],
      [feat],
      "pv.main",
    );
    // Residue 0: bases 27, 28, 29 (no gap); residue 1: 0, 1, 2; across the origin 29 → 0 is no gap either.
    expect([g.at(0).codon, g.at(0)["codon at"]]).toEqual(["GCA", "pl 28–30 (+)"]);
    expect([g.at(1).codon, g.at(1)["codon at"]]).toEqual(["ATG", "pl 1–3 (+)"]);
    g.t.ctx.dispose();
    const g2 = await opened(
      b.getDocument(),
      ((src: { stream: string | null }) => (src.stream === "feat" ? feat.provider : null)) as Loaders["stream"],
      [feat],
      "pv.copy",
    );
    // A copy of the plasmid mapped base for base (no codon step): its features reach it, its label has no codon.
    const label = g2.at(4);
    expect([label["codon at"] === undefined, Object.keys(label).length > 0]).toEqual([true, true]);
    g2.t.ctx.dispose();
  });
});

describe("case 9's coverage on its chains, through their rows", () => {
  it("asks the columns its residues reach and reads a residue's column's value", async () => {
    const { doc, loaders } = coverageCase9();
    const t = await composed(doc, { width: 1000, height: 400, loaders, section: "msa.chains" });
    const m = t.ctx.streams as Streams;
    for (let k = 0; k < 4; k++) {
      t.layout("msa.chains");
      m.pass();
      await flush();
    }
    // Each chain's block through its row's path (the section's blocks: the layer's own).
    const st = m.get("cov.a");
    const paths = st?.source.paths ?? [];
    const block = paths.findIndex((x) => x !== null);
    const path = paths[block];
    expect([block >= 0, paths.filter((x) => x !== null).length]).toEqual([true, 2]);
    // A residue's column through the row, and its value there (the server's i / 320).
    let checked = 0;
    const v = new Float64Array(4);
    for (let i = 0; i < 400 && checked < 5; i++) {
      const c: number[] = [];
      if (path) Paths.forEachFrom(path, i, i + 1, (a, e, acc) => acc.push(a, e), c);
      if (c.length !== 2) continue;
      v[0] = 0;
      if (!st?.values?.(block, i, i + 1, v) || v[3] === 0) continue;
      expect(v[2] / v[3]).toBeCloseTo(c[0] / 320, 3);
      checked++;
    }
    expect(checked).toBe(5);
    t.ctx.dispose();
  });
});
