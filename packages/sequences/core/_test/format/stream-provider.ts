/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The stream port's test double (seqspec:sources#the-loader): a scripted provider, its answers computed per request by
// a script, held until released when gated, every request and abort recorded. And the hg38-sized synthetic sources of
// ST1's bench: GRCh38's 24 chromosomes (3.1 G positions) as a letters stream's manifest, and a values stream over them
// whose bins are sampled on demand (a few samples a bin), so nothing is held and a request costs only its bins: never
// the procedural tile server, which caches every row it serves (047, RN3's review).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type {
  DataPatch,
  Manifest,
  SeqViewSpecDoc,
  StreamNode,
  StreamProvider,
  StreamRequest,
} from '../../src/format/index.js';
import { GeneticCode } from '../../src/model/index.js';

export interface Scripted {
  readonly provider: StreamProvider;
  /** Data requests in the order asked (the manifest's aside). */
  readonly requests: StreamRequest[];
  readonly aborted: StreamRequest[];
  /** Answers wait for `release` while true. */
  gated: boolean;
  /** Answers the first `n` held requests, in the order asked. */
  release(n?: number): void;
  readonly held: number;
}

const at = (r: StreamRequest) => ({
  source: r.source,
  system: r.system,
  start: r.start,
  stop: r.stop,
  row_start: r.row_start,
  row_stop: r.row_stop,
  level: r.level,
});

/** A provider over `manifest` answering each data request with `answer(r)` (an `Error`: that failure). */
export function scripted(
  manifest: Manifest,
  answer: (r: StreamRequest) => DataPatch | Error,
  o: { gated?: boolean } = {},
): Scripted {
  const waiting: (() => void)[] = [];
  const s: Scripted = {
    requests: [],
    aborted: [],
    gated: o.gated ?? false,
    release(n = Number.POSITIVE_INFINITY) {
      for (const go of waiting.splice(0, n)) go();
    },
    get held() {
      return waiting.length;
    },
    provider: {
      fetch(r, signal) {
        if (r.system === null) return Promise.resolve({ ...at(r), kind: 'manifest', manifest });
        s.requests.push(r);
        return new Promise<DataPatch>((ok, no) => {
          const go = () => {
            if (signal.aborted) return;
            const a = answer(r);
            if (a instanceof Error) no(a);
            else ok(a);
          };
          signal.addEventListener(
            'abort',
            () => {
              s.aborted.push(r);
              no(signal.reason);
            },
            { once: true },
          );
          if (s.gated) waiting.push(go);
          else queueMicrotask(go);
        });
      },
    },
  };
  return s;
}

/** GRCh38's chromosomes and their lengths: 3,088,269,832 positions. */
export const GRCH38: readonly (readonly [string, number])[] = [
  ['chr1', 248_956_422],
  ['chr2', 242_193_529],
  ['chr3', 198_295_559],
  ['chr4', 190_214_555],
  ['chr5', 181_538_259],
  ['chr6', 170_805_979],
  ['chr7', 159_345_973],
  ['chr8', 145_138_636],
  ['chr9', 138_394_717],
  ['chr10', 133_797_422],
  ['chr11', 135_086_622],
  ['chr12', 133_275_309],
  ['chr13', 114_364_328],
  ['chr14', 107_043_718],
  ['chr15', 101_991_189],
  ['chr16', 90_338_345],
  ['chr17', 83_257_441],
  ['chr18', 80_373_285],
  ['chr19', 58_617_616],
  ['chr20', 64_444_167],
  ['chr21', 46_709_983],
  ['chr22', 50_818_468],
  ['chrX', 156_040_895],
  ['chrY', 57_227_415],
];
/** BigWig-like zoom levels: positions per bin, finest first (a region is 256 bins: 8 kb to 268 Mb). */
export const LEVELS = [32, 256, 2048, 16_384, 131_072, 1_048_576];

/** The systems of a manifest over `chromosomes` (default GRCh38's). */
export const systemsOf = (chromosomes = GRCH38) =>
  chromosomes.map(([key, n]) => ({ key, extent: [0, n] as [number, number], numbering: null, label: key }));
export const lettersManifest = (chromosomes = GRCH38): Manifest => ({
  systems: systemsOf(chromosomes),
  levels: [],
  alphabet: 'dna',
});
export const valuesManifest = (chromosomes = GRCH38, levels = LEVELS): Manifest => ({
  systems: systemsOf(chromosomes),
  levels,
  domain: [0, 100],
});

/** A smooth synthetic score at position x of system s, in [5, 95]. */
export const score = (s: number, x: number) => 50 + 30 * Math.sin(x / 50_000 + s) + 15 * Math.sin(x / 700);

/**
 * The synthetic values answer of a request (`levels` the manifest's): bins with `min`, `max` and `mean` of four
 * samples each, or the data as a dense `value` column. Allocates only the patch.
 */
export function sampledValues(r: StreamRequest, keys: readonly string[], levels = LEVELS): DataPatch {
  const s = keys.indexOf(r.system as string);
  if (r.level < 0) {
    const value = new Float64Array(r.stop - r.start);
    for (let i = 0; i < value.length; i++) value[i] = score(s, r.start + i);
    return { ...at(r), kind: 'items', items: { value } };
  }
  const w = levels[r.level];
  const n = Math.ceil((r.stop - r.start) / w);
  const [min, max, mean] = [new Float64Array(n), new Float64Array(n), new Float64Array(n)];
  for (let j = 0; j < n; j++) {
    let lo = Number.POSITIVE_INFINITY;
    let hi = Number.NEGATIVE_INFINITY;
    let sum = 0;
    for (let q = 0; q < 4; q++) {
      const v = score(s, r.start + j * w + ((q + 0.5) * w) / 4);
      lo = Math.min(lo, v);
      hi = Math.max(hi, v);
      sum += v;
    }
    [min[j], max[j], mean[j]] = [lo, hi, sum / 4];
  }
  return { ...at(r), kind: 'bins', bins: { min, max, mean } };
}

/**
 * A 2bit-like synthetic genome's letter at position x of system s (ST3a): A, C, G, T by a hash of the place, runs of N
 * every 100 kb (1 kb each), soft-masked (lower-case) stretches of 300 every 2 kb.
 */
export function letterAt(s: number, x: number): string {
  if (x % 100_000 < 1000) return 'N';
  let h = Math.imul(x ^ (s * 0x9e3779b1), 0x85ebca6b);
  h = Math.imul(h ^ (h >>> 13), 0xc2b2ae35);
  const c = 'ACGT'[(h >>> 16) & 3];
  return x % 2000 < 300 ? c.toLowerCase() : c;
}
/** The synthetic letters answer of a data request: the region's letters, a string of `stop − start`. */
export function sampledLetters(r: StreamRequest, keys: readonly string[]): DataPatch {
  const s = keys.indexOf(r.system as string);
  let letters = '';
  for (let x = r.start; x < r.stop; x++) letters += letterAt(s, x);
  return { ...at(r), kind: 'letters', letters };
}

/** A hash of a place (fixtures): a non-negative 31-bit integer. */
const hash = (s: number, x: number) => {
  let h = Math.imul(x ^ (s * 0x9e3779b1), 0x85ebca6b);
  h = Math.imul(h ^ (h >>> 13), 0xc2b2ae35);
  return (h ^ (h >>> 16)) >>> 1;
};
/** An items manifest over `chromosomes` (ST3a-2): `count` bins at `levels`; its categorical fields' values (ST3b). */
export const itemsManifest = (
  chromosomes: readonly (readonly [string, number])[],
  levels: number[],
  vocabulary: Record<string, string[]> = {},
): Manifest => ({
  systems: systemsOf(chromosomes),
  levels,
  vocabulary,
});
/** A synthetic gene's `category` (ST3b), by a hash: the genes manifest's vocabulary, in its order. */
export const GENE_CATEGORIES = ['protein_coding', 'lncRNA', 'pseudogene'];
/**
 * Synthetic gene g of system s: about one each 37 kb, 8 to 30 kb long, 2 to 6 exons (rows of one `id`), the first 3 to
 * 6 kb (it carries the label), the others 0.4 to 1.3 kb.
 */
export function geneAt(s: number, g: number) {
  const start = g * 37_000 + (hash(s, g) % 9000);
  const length = 8000 + (hash(s, g + 7) % 22_000);
  const n = 2 + (hash(s, g + 13) % 5);
  const exons: [number, number][] = [];
  for (let k = 0; k < n; k++) {
    const a = start + Math.floor((k * length) / n);
    const w = k === 0 ? 3000 + (hash(s, g * 31) % 3000) : 400 + (hash(s, g * 31 + k) % 900);
    exons.push([a, Math.min(start + length, a + w)]);
  }
  const category = GENE_CATEGORIES[hash(s, g + 3) % 3];
  return { id: `gene-${s}-${g}`, label: `G${g}`, category, start, stop: start + length, exons };
}
/** A synthetic variant's alternative at x of system s (about one each kb, a second one at about a tenth), or "". */
export function variantAt(s: number, x: number, second = false): string {
  const h = hash(s + 99, x);
  return h % 997 === 0 && (!second || h % 10 === 0) ? 'ACGT'[(h >>> 10) & 3] : '';
}
/** A gene or variant source's answer (ST3a-2): its items over the region, or at a level a `count` per bin. */
export function sampledItems(
  r: StreamRequest,
  keys: readonly string[],
  kind: 'genes' | 'variants',
  levels: number[],
): DataPatch {
  const s = keys.indexOf(r.system as string);
  if (r.level >= 0) {
    const w = levels[r.level];
    const count = new Float64Array(Math.ceil((r.stop - r.start) / w));
    for (let j = 0; j < count.length; j++) {
      const a = r.start + j * w;
      // Genes starting in the bin; variants about one a kb, by a hash of the bin.
      count[j] =
        kind === 'genes'
          ? Math.max(0, Math.floor((a + w) / 37_000) - Math.floor(a / 37_000))
          : Math.round(w / 997) + (hash(s, a) % 3);
    }
    return { ...at(r), kind: 'bins', bins: { count } };
  }
  const start: number[] = [];
  const stop: number[] = [];
  const id: string[] = [];
  const label: string[] = [];
  const alt: string[] = [];
  const category: string[] = [];
  if (kind === 'genes')
    for (let g = Math.max(0, Math.floor((r.start - 40_000) / 37_000)); g * 37_000 < r.stop; g++) {
      const x = geneAt(s, g);
      if (x.stop <= r.start || x.start >= r.stop) continue;
      for (const [a, e] of x.exons) {
        start.push(a);
        stop.push(e);
        id.push(x.id);
        label.push(a === x.start ? x.label : '');
        category.push(x.category);
      }
    }
  else
    for (let x = r.start; x < r.stop; x++)
      for (const second of [false, true]) {
        const a = variantAt(s, x, second);
        if (a === '') continue;
        start.push(x);
        stop.push(x + 1);
        id.push(`var-${s}-${x}-${+second}`);
        label.push(`${a} at ${x + 1}`);
        alt.push(a);
      }
  return {
    ...at(r),
    kind: 'items',
    items: kind === 'genes' ? { start, stop, id, label, category } : { start, stop, id, label, alt },
  };
}

/**
 * The hg38-sized pair: `genome` (letters, its manifest only) and `phylop` (sampled values, or `answer`'s where it
 * gives one), each scripted.
 */
export function hg38Streams(
  o: {
    gated?: boolean;
    chromosomes?: readonly (readonly [string, number])[];
    answer?: (r: StreamRequest) => DataPatch | Error | null;
    /** The values manifest's dynamic nodes (ST1b-a). */
    nodes?: StreamNode[];
    /** `genome` answers its data with synthetic letters (ST3a); else every letters request fails. */
    letters?: boolean;
  } = {},
) {
  // ST3a-2: genes and variants (`streamedDoc({ items })`), gated with the rest.
  const keyed = (o.chromosomes ?? GRCH38).map((c) => c[0]);
  const items = (kind: 'genes' | 'variants', levels: number[]) =>
    scripted(
      itemsManifest(o.chromosomes ?? GRCH38, levels, VOCABULARY[kind]),
      (r) => sampledItems(r, keyed, kind, levels),
      {
        gated: o.gated,
      },
    );
  const genes = items('genes', GENE_LEVELS);
  const variants = items('variants', VARIANT_LEVELS);
  const chromosomes = o.chromosomes ?? GRCH38;
  const keys = chromosomes.map((c) => c[0]);
  const lettered = (r: StreamRequest) =>
    o.letters ? sampledLetters(r, keys) : new Error('the bench asks no letters without a sequence track');
  const genome = scripted(lettersManifest(chromosomes), lettered, { gated: o.gated });
  const values = (r: StreamRequest) => o.answer?.(r) ?? sampledValues(r, keys);
  const phylop = scripted({ ...valuesManifest(chromosomes), nodes: o.nodes }, values, { gated: o.gated });
  const by: Record<string, Scripted> = { genome, phylop, genes, variants };
  const stream = (src: { stream: string | null }) => (src.stream ? (by[src.stream]?.provider ?? null) : null);
  return { genome, phylop, genes, variants, stream };
}
/** The genes' and variants' categorical fields (ST3b). */
export const VOCABULARY = { genes: { category: GENE_CATEGORIES }, variants: { alt: ['A', 'C', 'G', 'T'] } };
/** The genes' and variants' `count` levels (ST3a-2): their data where a device px holds fewer positions than the first. */
export const GENE_LEVELS = [2048, 65_536, 1_048_576];
export const VARIANT_LEVELS = [1024, 32_768, 1_048_576];

/**
 * A genome browser's document over the pair: the complex from `genome`, a ruler, with `sequence` its letters (ST3a),
 * with `items` the genes and variants (ST3a-2), and `phylop`'s values unless `values` is false (`browser.main`).
 */
export function streamedDoc(
  o: { sequence?: boolean; items?: boolean; values?: boolean; variantsAs?: 'variants' | 'points' } = {},
): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'hg38' });
  const genome = b.data.source({ ref: 'genome', stream: 'genome', data_kind: 'letters' });
  const hg38 = b.data.complex({ ref: 'hg38', source_ref: genome, label: 'GRCh38' });
  b.data.source({ ref: 'phylop_src', stream: 'phylop', data_kind: 'values', coordinate_system: 'hg38' });
  const phylop = b.data.annotation({ ref: 'phylop', geometry: 'values', source_ref: 'phylop_src' });
  const v = b.view({ ref: 'browser', label: 'hg38', style: 'bare' });
  v.block({ complex: hg38 });
  v.track({ label: 'Ruler' }).ruler();
  if (o.sequence) v.track({ label: 'Sequence' }).letters({ color: b.colors.scheme('nucleotide') });
  if (o.items) {
    for (const k of ['genes', 'variants']) {
      b.data.source({ ref: `${k}_src`, stream: k, data_kind: 'items', coordinate_system: 'hg38' });
      b.data.annotation({ ref: k, geometry: 'ranges', source_ref: `${k}_src` });
    }
    v.track({ label: 'Genes', height: 4 }).ranges({ data: 'genes', labels: true });
    if (o.variantsAs === 'points') v.track({ label: 'Variants' }).points({ data: 'variants' });
    else v.track({ label: 'Variants' }).variants({ data: 'variants' });
  }
  if (o.values !== false) v.track({ label: 'phyloP' }).values({ data: phylop });
  return b.getDocument();
}

// ---- ST7: a protein on its genome ----
type Section = ReturnType<ReturnType<ReturnType<typeof createSeqViewSpec>['view']>['section']>;
/**
 * A coding gene on chr2 of the synthetic genome: three exons over 220 kb, each in its own region of 65,536; the
 * introns split codons 500 and 900.
 */
export const CODING_CDS: [number, number][] = [
  [100_010_000, 100_011_501],
  [100_090_000, 100_091_200],
  [100_230_000, 100_231_301],
];
/** Its protein: the synthetic letters of the CDS translated (a stop or an N read as X), 1,334 residues. */
export function codingProtein(): string {
  let dna = '';
  for (const [a, e] of CODING_CDS) for (let x = a; x < e; x++) dna += letterAt(1, x).toUpperCase();
  let protein = '';
  for (let k = 0; k + 3 <= dna.length; k += 3) {
    const r = GeneticCode.translate(1, dna.charCodeAt(k), dna.charCodeAt(k + 1), dna.charCodeAt(k + 2));
    protein += r < 0 || r === 42 ? 'X' : String.fromCharCode(r);
  }
  return protein;
}
/**
 * The genome's genes, variants and phyloP on the protein's axis, through its codon map (ST7): a section of the protein
 * alone, and with `genome` a section of chr2 at the gene beside it, the same tracks shown directly; `strand` "-" maps
 * the protein on the reverse strand.
 */
export function codingDoc(o: { genome?: boolean; strand?: '+' | '-' } = {}): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'A protein on its genome' });
  const genome = b.data.source({ ref: 'genome', stream: 'genome', data_kind: 'letters' });
  const hg38 = b.data.complex({ ref: 'hg38', source_ref: genome, label: 'GRCh38' });
  for (const k of ['genes', 'variants']) {
    b.data.source({ ref: `${k}_src`, stream: k, data_kind: 'items', coordinate_system: 'hg38' });
    b.data.annotation({ ref: k, geometry: 'ranges', source_ref: `${k}_src` });
  }
  b.data.source({ ref: 'phylop_src', stream: 'phylop', data_kind: 'values', coordinate_system: 'hg38' });
  b.data.annotation({ ref: 'phylop', geometry: 'values', source_ref: 'phylop_src' });
  const protein = b.data.sequence({ ref: 'p1', label: 'P1', alphabet: 'protein', value: codingProtein() });
  // On the minus strand (a test's) the same exons read downward, the highest first: the protein's residues stay (its
  // letters then code nothing exactly).
  const minus = o.strand === '-';
  const cds = minus ? CODING_CDS.slice().reverse() : CODING_CDS;
  b.data.mapping({ ref: 'p1.codon', from: protein, to: 'hg38/chr2', method: 'codon', cds, strand: o.strand ?? '+' });
  const v = b.view({ ref: 'coding', label: 'A protein on its genome' });
  const tracks = (s: Section) => {
    s.track({ label: '' }).ruler();
    s.track({ label: 'Genes', height: 3 }).ranges({ data: 'genes', labels: true });
    s.track({ label: 'Variants' }).variants({ data: 'variants' });
    s.track({ label: 'phyloP' }).values({ data: 'phylop' });
  };
  const p = v.section({ ref: 'coding.protein', label: "P1, with chr2's tracks through its codon map" });
  p.block({ sequence: protein });
  tracks(p);
  p.track({ label: 'P1' }).letters();
  if (o.genome) {
    const g = v.section({ ref: 'coding.genome', label: "chr2 at P1's gene, the same tracks directly" });
    g.block({ complex: hg38 });
    g.focus({
      ref: 'coding.focus',
      target: 'hg38/chr2',
      start: CODING_CDS[0][0] - 5000,
      stop: CODING_CDS[2][1] + 5000,
    });
    tracks(g);
  }
  return b.getDocument();
}

// ---- ST7: item coordinate systems ----
/** P1's transcript on chr2: its exons (the CDS) as rows of one item `T1` that carry its `map` (a codon map's params). */
export const TRANSCRIPT = 'T1';
/**
 * A genes answer with `T1` added, whole (every exon a row, each carrying `map`), in every region it meets, as a
 * transcripts file answers; `map` in place of its codon map over `CODING_CDS`, `cds` in place of its exons.
 */
export function transcriptItems(
  r: StreamRequest,
  keys: readonly string[],
  map: Record<string, unknown> | null = null,
  cds: readonly [number, number][] = CODING_CDS,
  id = TRANSCRIPT,
) {
  const p = sampledItems(r, keys, 'genes', GENE_LEVELS);
  if (p.kind !== 'items' || r.system !== 'chr2') return p;
  const t = p.items as Record<string, unknown[]>;
  const n = (t.start ?? []).length;
  t.map = Array.from({ length: n }, () => null);
  if (cds[cds.length - 1][1] <= r.start || cds[0][0] >= r.stop) return p;
  for (const [a, e] of cds) {
    const row = { start: a, stop: e, id, label: id, category: 'protein_coding' };
    for (const [k, v] of Object.entries(row)) t[k] = [...(t[k] ?? []), v];
    t.map.push(map ?? { method: 'codon', cds });
  }
  return p;
}
/**
 * P1's domains on the streamed transcript that codes it (ST7): P1 is mapped base for base onto the item coordinate
 * system `genes_src/T1` (the transcript's residues, made from its `map` when its region arrives); a genome section
 * shows the genes and P1's domains along T1's exons, and with `item` a section of T1's own residues shows the domains
 * and the genome's variants and phyloP through its codon map.
 */
export function transcriptDoc(
  o: { item?: boolean; named?: boolean; to?: string; target?: boolean; templates?: boolean; rows?: boolean } = {},
): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'A protein on its streamed transcript' });
  const genome = b.data.source({ ref: 'genome', stream: 'genome', data_kind: 'letters' });
  const hg38 = b.data.complex({ ref: 'hg38', source_ref: genome, label: 'GRCh38' });
  for (const k of ['genes', 'variants']) {
    b.data.source({ ref: `${k}_src`, stream: k, data_kind: 'items', coordinate_system: 'hg38' });
    b.data.annotation({ ref: k, geometry: 'ranges', source_ref: `${k}_src` });
  }
  b.data.source({ ref: 'phylop_src', stream: 'phylop', data_kind: 'values', coordinate_system: 'hg38' });
  b.data.annotation({ ref: 'phylop', geometry: 'values', source_ref: 'phylop_src' });
  const protein = b.data.sequence({ ref: 'p1', label: 'P1', alphabet: 'protein', value: codingProtein() });
  const t1 = `genes_src/${TRANSCRIPT}`;
  // `named` false: nothing names T1 (no item coordinate system is made); `to`: the mapping's end in its place (its
  // unit, or another sequence: a retarget); `target`: the domains target T1 itself, no mapping.
  const named = o.named !== false;
  b.data.sequence({ ref: 'p1copy', label: 'P1 copy', alphabet: 'protein', value: codingProtein() });
  if (named && !o.target)
    b.data.mapping({ ref: 'p1.t1', from: protein, to: o.to ?? t1, method: 'segments', segments: [[0, 0, null]] });
  const value = { start: [20, 480, 880], stop: [200, 620, 1000], label: ['Alpha', 'Beta', 'Gamma'] };
  b.data.annotation({ ref: 'domains', geometry: 'ranges', target: named && o.target ? t1 : protein, value });
  // `rows`: an annotation whose rows target T1 one by one (not supported yet: said).
  if (o.rows) b.data.annotation({ ref: 'by_row', geometry: 'ranges', value: { target: [t1], start: [0], stop: [5] } });
  const v = b.view({ ref: 'tx', label: 'A protein on its streamed transcript' });
  const g = v.section({ ref: 'tx.genome', label: 'chr2 at T1' });
  g.block({ complex: hg38 });
  g.focus({ ref: 'tx.focus', target: 'hg38/chr2', start: CODING_CDS[0][0] - 5000, stop: CODING_CDS[2][1] + 5000 });
  g.track({ label: '' }).ruler();
  g.track({ label: 'Genes', height: 3 }).ranges({ data: 'genes', labels: true });
  g.track({ label: 'P1 domains' }).ranges({ data: 'domains', labels: true });
  if (o.templates) {
    g.track({ label: 'Items', for_each: { kind: 'item_system', source: 'genes_src' } }).letters({ of: '$node' });
    g.track({ label: 'Sequences', for_each: { kind: 'sequence', source: 'genes_src' } }).letters({ of: '$node' });
  }
  if (o.item) {
    const s = v.section({ ref: 'tx.item', label: "T1's own residues" });
    s.block({ sequence: t1 });
    s.track({ label: '' }).ruler();
    s.track({ label: 'P1 domains' }).ranges({ data: 'domains', labels: true });
    s.track({ label: 'Variants' }).variants({ data: 'variants' });
    s.track({ label: 'phyloP' }).values({ data: 'phylop' });
  }
  return b.getDocument();
}
