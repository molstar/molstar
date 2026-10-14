/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Documents for the prepare tests (LY2), by value and DOM-free (Node and the browser): every R2 type on two chains, and
// format case 14's protein over its coding DNA, with overlapping features packed in lanes.

import { createSeqViewSpec, translate } from "@molstar/seq-core/builder";
import type { Features } from "../src/prepare.js";
import { segWindow } from "../src/pyramids.js";

const PROTEIN = "MSEEKLKELAEKAIRLAKEGNPEEALKLLEEALRLAKESGDPELIARVLELLAEVTKQSVHAYRLDDGKWFLQPGGRVSLEMAPTDLQAAH";
const residues = (n: number, seed = 0) =>
  Array.from({ length: n }, (_, i) => PROTEIN[(i * 7 + seed * 13) % PROTEIN.length]).join("");

/**
 * Every R2 type on two chains: a ruler, letters, labeled ranges, points, emphasis, variants, secondary structure,
 * values, color, and pairs within and between the chains (`arcs` as given).
 */
export function allTypes(
  o: { length?: number; arcs?: "always" | "hover"; side?: "above" | "below"; layout?: "tracks" | "text" } = {},
) {
  const n = o.length ?? 300;
  const m = Math.round(n * 0.7);
  const b = createSeqViewSpec({ title: "All types" });
  const a = b.data.sequence({ ref: "a", label: "Chain A", value: residues(n, 1) });
  const c = b.data.sequence({ ref: "c", label: "Chain C", value: residues(m, 2) });
  const ann = (ref: string, geometry: string, target: unknown, value: unknown) =>
    b.data.annotation({ ref, geometry, target: target as never, value } as never);
  const spread = (k: number, len: number) => (k * 37) % (len - 40);
  const ra = ann("ra", "ranges", a, {
    start: Array.from({ length: 16 }, (_, k) => spread(k, n)),
    stop: Array.from({ length: 16 }, (_, k) => spread(k, n) + 4 + ((k * 7) % 30)),
    label: Array.from({ length: 16 }, (_, k) => `Feature ${k}`),
    category: Array.from({ length: 16 }, (_, k) => `c${k % 4}`),
  });
  const ss = ann("ss", "ranges", a, {
    start: Array.from({ length: 12 }, (_, k) => k * 24),
    stop: Array.from({ length: 12 }, (_, k) => k * 24 + 6 + (k % 3) * 4),
    category: Array.from({ length: 12 }, (_, k) => ["helix", "strand", "turn"][k % 3]),
  });
  const va = ann("var", "ranges", a, {
    start: Array.from({ length: 40 }, (_, k) => (k * 13) % n),
    alt: Array.from({ length: 40 }, (_, k) => "ACDEFGHIKL"[k % 10]),
  });
  const vals = ann("vals", "values", a, { value: Array.from({ length: n }, (_, i) => ((i * 31) % 100) / 100) });
  const pairs = ann("pp", "pairs", null, {
    target: Array.from({ length: 24 }, (_, k) => (k % 3 === 0 ? "c" : "a")),
    start: Array.from({ length: 24 }, (_, k) => (k * 11) % (m - 1)),
    b_target: Array.from({ length: 24 }, () => "a"),
    b_start: Array.from({ length: 24 }, (_, k) => (k * 11 + 7 + k * 3) % n),
    category: Array.from({ length: 24 }, (_, k) => `k${k % 3}`),
  });
  const v = b.view({ ref: "v", label: "All", height: "box", ...(o.layout ? { layout: o.layout } : {}) });
  v.block({ sequence: a }).block({ sequence: c });
  v.track({ label: "" }).ruler();
  v.track({ label: "Letters" }).letters();
  v.track({ label: "Contacts" }).pairs({ data: pairs, arcs: o.arcs ?? "always", side: o.side ?? null } as never);
  v.track({ label: "Features" }).ranges({ data: ra, labels: true });
  v.track({ label: "Sites" }).points({ data: va });
  v.track({ label: "Emphasis" }).emphasis({ data: va });
  v.track({ label: "Variants" }).variants({ data: va });
  v.track({ label: "Secondary" }).secondary({ data: ss });
  v.track({ label: "Values" }).values({ data: vals });
  v.track({ label: "Color" }).color({ data: vals });
  return b.getDocument();
}

/** Format case 14: a protein over its coding DNA, a codon split by the intron; residues `form`. */
export function gene(form: "auto" | "one" | "three" = "auto") {
  const b = createSeqViewSpec({ title: "A protein over its coding DNA" });
  const GENE = [
    "GAATTCGCCACC", // 5′ UTR
    "ATGGCTTCCAAGGAACTGACCCCACGTGACGGCC", // exon 1's coding part
    "GTAAGTCTTTCCTAACTTCTGTTTCCACAG", // the intron
    "ATTGGTACTTCAACAGCCAGGTGGAGTAA", // exon 2's coding part and stop
    "GCGGCCGCATAAAGT", // 3′ UTR
  ].join("");
  const g = b.data.sequence({ ref: "gene", label: "Gene fragment", alphabet: "dna", value: GENE });
  const cds: [number, number][] = [
    [12, 46],
    [76, 105],
  ];
  const protein = translate(b.data, { ref: "protein", label: "Protein", from: g, cds });
  const v = b.view({ ref: "gene-view", label: "Gene" }).block({ sequence: g });
  v.track({ label: "" }).ruler();
  v.track({ label: "DNA" }).letters();
  v.track({ ref: "protein-track", label: "Protein" }).letters({ of: protein, residue_form: form });
  return b.getDocument();
}

/** One chain of `n` residues with ranges at `[start, stop)` (lanes packed in prepare), points at them too. */
export function rangesDoc(ranges: readonly (readonly [number, number])[], n = 60) {
  const b = createSeqViewSpec({ title: "Ranges" });
  const a = b.data.sequence({ ref: "a", label: "Chain A", value: residues(n, 1) });
  const r = b.data.annotation({
    ref: "r",
    geometry: "ranges",
    target: a,
    value: { start: ranges.map((x) => x[0]), stop: ranges.map((x) => x[1]), label: ranges.map((_, k) => `R${k}`) },
  });
  const v = b.view({ ref: "v", label: "Ranges" });
  v.block({ sequence: a });
  v.track({ label: "Letters" }).letters();
  v.track({ ref: "features", label: "Features" }).ranges({ data: r, labels: true });
  v.track({ ref: "sites", label: "Sites" }).points({ data: r });
  v.track({ ref: "vals", label: "Values" }).values({
    data: b.data.annotation({
      ref: "vv",
      geometry: "values",
      target: a,
      value: { value: Array.from({ length: n }, (_, i) => i / n) },
    }),
  });
  return b.getDocument();
}

/**
 * A large section for the index for scale (G1): `blocks` DNA blocks of `length` nucleotides with bands (`color` by
 * category), dense values, variants, ranges in lanes (many per device column zoomed out) and repeats. Seeded and
 * deterministic.
 */
export function scaleDoc(o: { length?: number; blocks?: number; ranges?: number; variants?: number } = {}) {
  const n = o.length ?? 200_000;
  let s = 7;
  const rnd = () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 2 ** 32;
  };
  const int = (lo: number, hi: number) => lo + Math.floor(rnd() * (hi - lo));
  const b = createSeqViewSpec({ title: "Scale" });
  const view = b.view({ ref: "v", label: "Scale", height: "box" });
  const blocks = Array.from({ length: o.blocks ?? 2 }, (_, k) => {
    const seq = Array.from({ length: n }, () => "ACGT"[int(0, 4)]).join("");
    return b.data.sequence({ ref: `chr${k}`, label: `chr${k}`, alphabet: "dna", value: seq });
  });
  const table = (ref: string, geometry: string, rows: number, len: [number, number], field: string, cats: string[]) => {
    const t = { target: [] as string[], start: [] as number[], stop: [] as number[], [field]: [] as string[] };
    for (let r = 0; r < rows; r++) {
      const k = int(0, blocks.length);
      const l = int(len[0], len[1]);
      const at = int(0, n - l - 1);
      t.target.push(`chr${k}`);
      t.start.push(at);
      t.stop.push(at + l);
      (t[field] as string[]).push(cats[int(0, cats.length)]);
    }
    if (len[1] <= 1) delete (t as { stop?: number[] }).stop;
    return b.data.annotation({ ref, geometry, value: t } as never);
  };
  const bands = table("bands", "ranges", 200, [500, 5000], "category", ["gneg", "gpos50", "gpos100"]);
  const genes = table("genes", "ranges", o.ranges ?? 3000, [50, 3000], "category", ["UTR", "CDS"]);
  const variants = table("variants", "ranges", o.variants ?? 4000, [1, 1], "alt", ["A", "C", "G", "T"]);
  for (const c of blocks) view.block({ sequence: c });
  view.track({ label: "" }).ruler();
  view.track({ label: "Bands" }).color({ data: bands });
  view.track({ label: "Letters" }).letters();
  view.track({ label: "Genes", height: 3 }).ranges({ data: genes });
  const values = view.track({ label: "Values" });
  blocks.forEach((c, k) => {
    const value = Float64Array.from({ length: n }, (_, i) => 0.5 + 0.4 * Math.sin(i / (300 + 50 * k)) * rnd());
    values.values({ data: b.data.annotation({ ref: `vals${k}`, geometry: "values", target: c, value: { value } }) });
  });
  view.track({ label: "Variants" }).variants({ data: variants });
  return b.getDocument();
}

/**
 * `color` layers with and without values: "Bands" over two ranges with no `value` column, [10, 20) and [30, 40);
 * "Shade" over dense values i / n.
 */
export function colorDoc(n = 60) {
  const b = createSeqViewSpec({ title: "Color" });
  const a = b.data.sequence({ ref: "a", label: "Chain A", value: residues(n, 1) });
  const value = { start: [10, 30], stop: [20, 40], label: ["Band 1", "Band 2"] };
  const bands = b.data.annotation({ ref: "bands", geometry: "ranges", target: a, value });
  const shade = Array.from({ length: n }, (_, i) => i / n);
  const vals = b.data.annotation({ ref: "shade", geometry: "values", target: a, value: { value: shade } });
  const v = b.view({ ref: "v", label: "Color" });
  v.block({ sequence: a });
  v.track({ label: "Bands" }).color({ data: bands });
  v.track({ label: "Shade" }).color({ data: vals });
  return b.getDocument();
}

/**
 * Transcripts as items: fragments sharing an `id` (exons), two transcripts overlapping in extent and a third clear of
 * both, on one chain; `blocks` chains side by side (the block gap shrinks past 8).
 */
export function transcripts(o: { blocks?: number } = {}) {
  const b = createSeqViewSpec({ title: "Transcripts" });
  const chains = Array.from({ length: o.blocks ?? 1 }, (_, k) =>
    b.data.sequence({ ref: `c${k}`, label: `C${k}`, value: residues(200, k) }),
  );
  const exons = [
    ["t1", 10, 30],
    ["t1", 60, 80],
    ["t1", 120, 140],
    ["t2", 40, 50],
    ["t2", 90, 110],
    ["t3", 150, 160],
    ["t3", 170, 190],
  ] as const;
  const t = b.data.annotation({
    ref: "tx",
    geometry: "ranges",
    target: chains[0],
    value: { id: exons.map((x) => x[0]), start: exons.map((x) => x[1]), stop: exons.map((x) => x[2]) },
  });
  const v = b.view({ ref: "v", label: "Transcripts" });
  for (const c of chains) v.block({ sequence: c });
  v.track({ label: "Genes" }).ranges({ data: t });
  return b.getDocument();
}

/** Four bands as a discrete palette over [0, 100] (AlphaFold DB's pLDDT colors). */
export const BANDS = {
  kind: "discrete",
  mode: "absolute",
  colors: [
    ["#FF7D45", 0],
    ["#FFDB13", 50],
    ["#65CBF3", 70],
    ["#0053D6", 90],
  ],
} as const;
/** Each band's score in `linesDoc`, a run of 10 positions each in turn; fractional, so no value is a small integer. */
export const BAND_SCORES = [20.5, 60.5, 80.5, 95.5] as const;
/**
 * `values` colored by their values on two chains of `n` and 0.7 `n` residues: "Line" a line in four bands over
 * [0, 100], "Area" an area under a continuous palette, "Plain" a line with no color theme; the scores step through
 * the bands (`BAND_SCORES`, 10 positions each) with a missing value every 50; `alternate`: the first and the last band's
 * scores in turn, position by position (any 4 or more average into the second band); `style` circular: on rings.
 */
export function linesDoc(n = 400, o: { style?: "circular"; alternate?: boolean } = {}) {
  const style = o.style;
  const b = createSeqViewSpec({ title: "Colored lines" });
  const chains = [n, Math.round(n * 0.7)].map((len, k) =>
    b.data.sequence({ ref: `c${k}`, label: `Chain ${k}`, value: residues(len, k) }),
  );
  const scores = chains.map((c, k) => {
    const len = k === 0 ? n : Math.round(n * 0.7);
    const value = Float64Array.from({ length: len }, (_, i) =>
      i % 50 === 49 ? Number.NaN : BAND_SCORES[o.alternate ? 3 * (i % 2) : Math.floor(i / 10) % 4],
    );
    return b.data.annotation({ ref: `s${k}`, geometry: "values", target: c, value: { value } });
  });
  const v = b.view({ ref: "v", label: "Lines", ...(style ? { style } : {}) });
  for (const c of chains) v.block({ sequence: c });
  const palette = (p: unknown) => ({ theme: "field_palette", params: { field: "value", palette: p } });
  v.track({ label: "Line" }).values({ data: scores, encoding: "line", domain: [0, 100], color: palette(BANDS) });
  v.track({ label: "Area" }).values({
    data: scores,
    encoding: "area",
    domain: [0, 100],
    color: palette({ kind: "continuous", colors: "Viridis", value_domain: [0, 100] }),
  });
  v.track({ label: "Plain" }).values({ data: scores, encoding: "line", domain: [0, 100] });
  return b.getDocument();
}

const win = new Int32Array(2);
/**
 * The items of block b whose lanes' windows meet [from, to) on display row `row` (column data counts on every row),
 * stopping at `limit`: the segment tables' windows, counted (prepare's tests and guard; no draw reads it). Allocation-free.
 */
export function itemsIn(f: Features, b: number, from: number, to: number, row: number, limit: number): number {
  let n = 0;
  for (let lane = 0; lane < f.lanes && n < limit; lane++) {
    const g = b * f.lanes + lane;
    segWindow(f.start, f.reach, f.seg[g], f.seg[g + 1], from, to, win);
    if (!f.hasRows) n += win[1] - win[0];
    else for (let k = win[0]; k < win[1] && n < limit; k++) if (f.on[k] < 0 || f.on[k] === row) n++;
  }
  return n < limit ? n : limit;
}
