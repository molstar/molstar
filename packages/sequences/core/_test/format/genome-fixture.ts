/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A synthetic genome for the Seq* viewer's scale test: 24 chromosomes sized as the human ones, scaled to about
// 3,000,000 nucleotides, with cytobands, genes as transcripts (exons as fragments of one item, UTR and CDS), proteins
// translated from a subset of the CDSs (codon mappings, as case 14), GC content per position, variants and repeats.
// Seeded and generated where it is used (the format's tests, Storybook's genome story); nothing fetched, nothing
// committed.

import { createSeqViewSpec, type NodeHandle, translate } from "@molstar/seq-core/builder";
import type { SeqViewSpecDoc } from "../../src/format/index.ts";

/** GRCh38's chromosome lengths in Mb (chr1–22, X, Y). */
const HUMAN_MB = [
  248.96, 242.19, 198.3, 190.21, 181.54, 170.81, 159.35, 145.14, 138.39, 133.8, 135.09, 133.28, 114.36, 107.04, 101.99,
  90.34, 83.26, 80.37, 58.62, 64.44, 46.71, 50.82, 156.04, 57.23,
];
const NAMES = [...Array.from({ length: 22 }, (_, i) => `chr${i + 1}`), "chrX", "chrY"];

/** mulberry32: a small seeded PRNG, uniform in [0, 1). */
function prng(seed: number) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export interface GenomeOptions {
  /** Total nucleotides (default 3,000,000), split over 24 chromosomes in the human proportions. */
  readonly total?: number;
  readonly genes?: number;
  /** Genes whose CDS is translated into a protein with a codon mapping (default 300). */
  readonly proteins?: number;
  readonly variants?: number;
  readonly repeats?: number;
  readonly seed?: number;
}

export interface Chromosome {
  readonly name: string;
  readonly length: number;
  readonly seq: string;
}
interface Exon {
  start: number;
  stop: number;
}
export interface Gene {
  readonly chr: number;
  readonly name: string;
  readonly exons: readonly Exon[];
  /** The coding ranges on the chromosome, in reading order (strand +), total a multiple of 3. */
  readonly cds: [number, number][];
}
export interface SynthGenome {
  readonly chromosomes: readonly Chromosome[];
  readonly genes: readonly Gene[];
  /** Indices into `genes` of those translated (their CDSs written as clean ORFs, never overlapping each other). */
  readonly coding: readonly number[];
  /** Cytobands per chromosome: [chr, start, stop, stain]. */
  readonly bands: { chr: number[]; start: number[]; stop: number[]; stain: string[] };
  /** GC fraction per position, of the position's 1 kb window, chromosome after chromosome. */
  readonly gc: Float64Array;
  readonly variants: { chr: number[]; start: number[]; alt: string[] };
  readonly repeats: { chr: number[]; start: number[]; stop: number[]; family: string[] };
}

const BASES = "ACGT";
/** Sense codons (no stop) of the standard code, for clean ORFs. */
const SENSE: string[] = [];
for (const a of BASES)
  for (const b of BASES) for (const c of BASES) if (!["TAA", "TAG", "TGA"].includes(a + b + c)) SENSE.push(a + b + c);
const STAINS = ["gneg", "gpos25", "gpos50", "gpos75", "gpos100"];
const FAMILIES = ["LINE", "SINE", "LTR", "DNA", "Simple"];

/** The data, deterministic for a seed. */
export function synthGenome(o: GenomeOptions = {}): SynthGenome {
  const total = o.total ?? 3_000_000;
  const rand = prng(o.seed ?? 20261009);
  const int = (lo: number, hi: number) => lo + Math.floor(rand() * (hi - lo));
  const sum = HUMAN_MB.reduce((a, b) => a + b, 0);
  const lengths = HUMAN_MB.map((mb) => Math.round((mb / sum) * total));
  const bytes = lengths.map((n) => new Uint8Array(n));

  // Nucleotides with GC drifting over isochores of 20–200 kb (0.35–0.60), so the GC track has structure.
  for (const buf of bytes) {
    let gc = 0.41;
    let next = 0;
    for (let i = 0; i < buf.length; i++) {
      if (i === next) {
        gc = 0.35 + rand() * 0.25;
        next = i + int(20_000, 200_000);
      }
      const r = rand();
      buf[i] = r < gc ? (r < gc / 2 ? 67 : 71) : r < (1 + gc) / 2 ? 65 : 84;
    }
  }

  // Cytobands: alternating gneg and a G-positive stain, 1–5 % of the chromosome each, a centromere (acen) pair.
  const bands = { chr: [] as number[], start: [] as number[], stop: [] as number[], stain: [] as string[] };
  lengths.forEach((n, c) => {
    const cen = Math.floor(n * (c < 3 || c === 22 ? 0.45 : 0.3));
    let at = 0;
    let k = 0;
    while (at < n) {
      let stop = Math.min(n, at + int(Math.floor(n / 100), Math.floor(n / 20)));
      if (at < cen && stop > cen) stop = cen;
      const acen = stop === cen || at === cen;
      bands.chr.push(c);
      bands.start.push(at);
      bands.stop.push(stop);
      bands.stain.push(acen ? "acen" : k % 2 === 0 ? "gneg" : STAINS[int(1, STAINS.length)]);
      at = stop;
      k++;
    }
  });

  // Genes, proportional to length: a log-uniform span of 1–30 kb, 2–12 exons of 60–400 nt, UTRs at the ends.
  const geneCount = o.genes ?? 3000;
  const cum = lengths.map(
    (
      (s) => (n: number) =>
        (s += n)
    )(0),
  );
  const genes: Gene[] = [];
  for (let g = 0; g < geneCount; g++) {
    const x = rand() * total;
    const chr = cum.findIndex((c) => x < c);
    const n = lengths[chr];
    const span = Math.min(n - 1, Math.round(Math.exp(Math.log(1000) + rand() * Math.log(30))));
    const start = int(0, Math.max(1, n - span));
    const count = int(2, Math.min(13, Math.floor(span / 300) + 2));
    // Exon starts spread over the span: cut it into `count` slots, an exon at the start of each.
    const exons: Exon[] = [];
    const slot = span / count;
    for (let e = 0; e < count; e++) {
      const s = start + Math.floor(e * slot) + (e === 0 ? 0 : int(0, Math.floor(slot / 3)));
      const len = Math.min(int(60, 400), Math.floor(slot * 0.6));
      exons.push({ start: s, stop: s + Math.max(30, len) });
    }
    exons[count - 1].stop = Math.min(n, Math.max(exons[count - 1].stop, start + span));
    // The CDS: from inside the first exon to inside the last, cut to whole codons.
    const first = exons[0];
    const last = exons[count - 1];
    const cdsStart = first.start + int(5, Math.max(6, Math.floor((first.stop - first.start) / 2)));
    const cdsStop = last.stop - int(5, Math.max(6, Math.floor((last.stop - last.start) / 2)));
    const cds: [number, number][] = [];
    for (const ex of exons) {
      const a = Math.max(ex.start, cdsStart);
      const b = Math.min(ex.stop, cdsStop);
      if (a < b) cds.push([a, b]);
    }
    let len = cds.reduce((t, [a, b]) => t + b - a, 0);
    while (len % 3 !== 0 || cds[cds.length - 1][1] - cds[cds.length - 1][0] < 3) {
      const tail = cds[cds.length - 1];
      tail[1]--;
      len--;
      if (tail[1] <= tail[0]) {
        cds.pop();
        len = cds.reduce((t, [a, b]) => t + b - a, 0);
      }
    }
    genes.push({ chr, name: `G${String(g + 1).padStart(4, "0")}`, exons, cds });
  }

  // Coding genes: spread over the list, never overlapping each other; their CDSs written as ATG, sense codons, a stop.
  const coding: number[] = [];
  const want = Math.min(o.proteins ?? 300, geneCount);
  const taken: [number, number, number][] = [];
  // Every `stride`-th gene first, then the others, until `want` are taken.
  const stride = Math.max(1, Math.floor(geneCount / want));
  const order = [...genes.keys()].sort((x, y) => (x % stride) - (y % stride) || x - y);
  for (const g of order) {
    if (coding.length >= want) break;
    const gene = genes[g];
    const a = gene.exons[0].start;
    const b = gene.exons[gene.exons.length - 1].stop;
    if (taken.some(([c, s, e]) => c === gene.chr && s < b && a < e)) continue;
    taken.push([gene.chr, a, b]);
    coding.push(g);
    const buf = bytes[gene.chr];
    const nts = gene.cds.reduce((t, [x, y]) => t + y - x, 0);
    let codon = "";
    let k = 0;
    for (const [x, y] of gene.cds)
      for (let i = x; i < y; i++, k++) {
        if (k % 3 === 0)
          codon = k === 0 ? "ATG" : k + 3 >= nts ? ["TAA", "TAG", "TGA"][int(0, 3)] : SENSE[int(0, SENSE.length)];
        buf[i] = codon.charCodeAt(k % 3);
      }
  }

  // GC per position: its 1 kb window's fraction (the windows' values, dense, as `values` takes them).
  const gc = new Float64Array(total + 24);
  let off = 0;
  for (const buf of bytes) {
    for (let w = 0; w < buf.length; w += 1000) {
      const end = Math.min(buf.length, w + 1000);
      let n = 0;
      for (let i = w; i < end; i++) if (buf[i] === 67 || buf[i] === 71) n++;
      gc.fill(n / (end - w), off + w, off + end);
    }
    off += buf.length;
  }

  // Variants: single positions, an alternative base other than the reference.
  const variants = { chr: [] as number[], start: [] as number[], alt: [] as string[] };
  for (let k = 0; k < (o.variants ?? 20_000); k++) {
    const x = rand() * total;
    const chr = cum.findIndex((c) => x < c);
    const i = int(0, lengths[chr]);
    const ref = bytes[chr][i];
    let alt = BASES.charCodeAt(int(0, 4));
    if (alt === ref) alt = BASES.charCodeAt((BASES.indexOf(String.fromCharCode(ref)) + 1) % 4);
    variants.chr.push(chr);
    variants.start.push(i);
    variants.alt.push(String.fromCharCode(alt));
  }

  // Repeats: 50–600 nt, five families.
  const repeats = { chr: [] as number[], start: [] as number[], stop: [] as number[], family: [] as string[] };
  for (let k = 0; k < (o.repeats ?? 10_000); k++) {
    const x = rand() * total;
    const chr = cum.findIndex((c) => x < c);
    const len = int(50, 600);
    const s = int(0, lengths[chr] - len);
    repeats.chr.push(chr);
    repeats.start.push(s);
    repeats.stop.push(s + len);
    repeats.family.push(FAMILIES[int(0, FAMILIES.length)]);
  }

  const dec = new TextDecoder("latin1");
  const chromosomes = bytes.map((b, c) => ({ name: NAMES[c], length: b.length, seq: dec.decode(b) }));
  return { chromosomes, genes, coding, bands, gc: gc.subarray(0, off), variants, repeats };
}

// #region genome
/** The document: 24 chromosome blocks on one axis, each track's data one table with a `target` column. */
export function genomeDoc(g: SynthGenome, o: { gcTable?: "one" | "each" } = {}): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: "A synthetic genome of 3,000,000 nucleotides" });
  const chrs = g.chromosomes.map((c) => b.data.sequence({ ref: c.name, label: c.name, alphabet: "dna", value: c.seq }));
  // A string column dictionary-coded (seqspec:format#positions): 24 names, a code per row. Number columns are typed
  // arrays, kept as given (seqspec:format#values): JSON writes them as arrays, MessagePack packs them.
  const target = (chr: ArrayLike<number>) => ({ dict: g.chromosomes.map((c) => c.name), codes: Uint8Array.from(chr) });
  const i32 = (xs: number[]) => Int32Array.from(xs);

  const bands = b.data.annotation({
    ref: "cytobands",
    geometry: "ranges",
    value: { target: target(g.bands.chr), start: i32(g.bands.start), stop: i32(g.bands.stop), category: g.bands.stain },
  });

  // Transcripts: one item per gene (`id`), its exons as fragments, each cut into UTR and CDS parts (`category`).
  const tx = { target: [] as number[], start: [] as number[], stop: [] as number[], id: [] as string[] };
  const cat: string[] = [];
  const label: string[] = [];
  const part = (gene: Gene, s: number, e: number, c: string) => {
    if (e <= s) return;
    tx.target.push(gene.chr);
    tx.start.push(s);
    tx.stop.push(e);
    tx.id.push(gene.name);
    cat.push(c);
    label.push(gene.name);
  };
  for (const gene of g.genes) {
    const c0 = gene.cds[0][0];
    const c1 = gene.cds[gene.cds.length - 1][1];
    for (const ex of gene.exons) {
      part(gene, ex.start, Math.min(ex.stop, c0), "UTR");
      part(gene, Math.max(ex.start, c0), Math.min(ex.stop, c1), "CDS");
      part(gene, Math.max(ex.start, c1), ex.stop, "UTR");
    }
  }
  const genes = b.data.annotation({
    ref: "genes",
    geometry: "ranges",
    value: { ...tx, target: target(tx.target), start: i32(tx.start), stop: i32(tx.stop), category: cat, label },
  });

  // Proteins: each coding gene's CDS translated, with its codon mapping (case 14's `translate`).
  for (const k of g.coding) {
    const gene = g.genes[k];
    translate(b.data, {
      ref: `${gene.name}.p`,
      label: `${gene.name} protein`,
      from: chrs[gene.chr],
      cds: gene.cds,
    });
  }

  // GC: dense values (no `index`), one per position in index order. One table per chromosome (`target` a param), or
  // one table for all with a `target` per row (3M refs, dictionary-coded; slower to check and load).
  const gc: NodeHandle<"annotation">[] = [];
  if (o.gcTable === "one") {
    const codes = new Uint8Array(g.gc.length);
    let at = 0;
    g.chromosomes.forEach((c, k) => {
      codes.fill(k, at, at + c.length);
      at += c.length;
    });
    gc.push(b.data.annotation({ ref: "gc", geometry: "values", value: { target: target(codes), value: g.gc } }));
  } else {
    let at = 0;
    g.chromosomes.forEach((c, k) => {
      const value = g.gc.subarray(at, at + c.length);
      at += c.length;
      gc.push(b.data.annotation({ ref: `gc.${c.name}`, geometry: "values", target: chrs[k], value: { value } }));
    });
  }
  const variants = b.data.annotation({
    ref: "variants",
    geometry: "ranges",
    value: { target: target(g.variants.chr), start: i32(g.variants.start), alt: g.variants.alt },
  });
  const rp = g.repeats;
  const repeats = b.data.annotation({
    ref: "repeats",
    geometry: "ranges",
    value: { target: target(rp.chr), start: i32(rp.start), stop: i32(rp.stop), category: rp.family },
  });

  const v = b.view({ ref: "genome", label: "Genome" });
  for (const c of chrs) v.block({ sequence: c });
  v.track({ label: "" }).ruler();
  v.track({ label: "Cytobands" }).color({
    data: bands,
    color: {
      theme: "field_palette",
      params: {
        palette: {
          kind: "categorical",
          colors: {
            gneg: "#f4f4f4",
            gpos25: "#c8c8c8",
            gpos50: "#969696",
            gpos75: "#5f5f5f",
            gpos100: "#1e1e1e",
            acen: "#c0392b",
          },
        },
      },
    },
  });
  v.track({ label: "Nucleotides" }).letters();
  // One layer for every protein: each over its codons on its chromosome (`of: "mapped"`).
  if (g.coding.length > 0) v.track({ label: "Proteins" }).letters({ of: "mapped" });
  v.track({ label: "Genes", height: 4 }).ranges({ data: genes, labels: true });
  const gcTrack = v.track({ label: "GC content", height: 2 });
  for (const data of gc) gcTrack.values({ data, domain: [0.3, 0.65] });
  v.track({ label: "Variants" }).variants({ data: variants });
  v.track({ label: "Repeats" }).color({ data: repeats });
  return b.getDocument();
}
// #endregion
