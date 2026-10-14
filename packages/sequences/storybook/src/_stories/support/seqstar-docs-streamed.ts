/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The documents and providers of the streaming stories (`Seq*/Dynamic nodes`, `Seq*/Streamed data`, and the atlas of
// `Seq*/Performance`), kept out of the story files so the smoke test (packages/sequences/storybook/_test/smoke.test.ts)
// can open each on a headless viewer. The providers are the committed fixtures of the format's stream-provider.ts.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import {
  type DataPatch,
  type Loaders,
  type Manifest,
  type SeqViewSpecDoc,
  type StreamNode,
  type StreamRequest,
  TILE,
} from '@molstar/seq-core/format';
import { type TileTransport, tilesStream } from '@molstar/seq-core/stream/tiles';
import { Hit, type SeqViewerContext } from '@molstar/seq-viewer';
import {
  codingDoc,
  GENE_LEVELS,
  GRCH38,
  hg38Streams,
  itemsManifest,
  LEVELS,
  lettersManifest,
  sampledItems,
  sampledLetters,
  sampledValues,
  scripted,
  transcriptDoc,
  transcriptItems,
  VARIANT_LEVELS,
  VOCABULARY,
  valuesManifest,
} from '../../../../core/_test/format/stream-provider.js';

// ---- Seq*/Dynamic nodes ----
const DYNAMIC_CHR: [string, number][] = [['chrT', 2_000_000]];
export const PEAKS = { A: 300_000, B: 1_000_000, C: 1_700_000 } as const;
const SPAN = 4096;
const genes: StreamNode = {
  kind: 'annotation',
  key: 'genes',
  params: { geometry: 'ranges', target: 'chrT/chrT', value: { start: [200_000, 950_000], stop: [420_000, 1_100_000] } },
};
/** A region of data's peak: A, B or C where the region holds one, else one named by the region. */
function peak(r: StreamRequest): StreamNode {
  const named = Object.entries(PEAKS).find(([, at]) => at >= r.start && at < r.stop)?.[0];
  const at = named ? PEAKS[named as keyof typeof PEAKS] : r.start + 1000;
  const value = { start: [at - 200], stop: [at + 200], label: [named ? `peak ${named}` : 'peak'] };
  return {
    kind: 'annotation',
    key: named ?? `peak-${r.start / SPAN}`,
    params: { geometry: 'ranges', target: 'chrT/chrT', value },
  };
}

/** The provider, its answers delayed by the knob's latency (read when asked, so the loaders stay one object). */
export const dynamicKnobs = { latency: 300 };
export function dynamicLoaders(): Loaders {
  const genome = scripted(lettersManifest(DYNAMIC_CHR), () => new Error('no letters here'));
  const answer = (r: StreamRequest): DataPatch => {
    const p = sampledValues(r, ['chrT']);
    return r.level < 0 ? { ...p, nodes: [peak(r)] } : p;
  };
  const signal = scripted({ ...valuesManifest(DYNAMIC_CHR, [32, 256, 2048, 16_384]), nodes: [genes] }, answer);
  const late = (inner: typeof signal.provider): typeof signal.provider => ({
    fetch: (r, sig) =>
      new Promise((ok) => setTimeout(ok, r.system === null ? 0 : dynamicKnobs.latency)).then(() => inner.fetch(r, sig)),
  });
  const by: Record<string, typeof signal.provider> = { genome: genome.provider, signal: late(signal.provider) };
  return {
    load: () => Promise.reject(new Error('no files')),
    stream: (src) => (src.stream ? (by[src.stream] ?? null) : null),
  };
}

export function dynamicDoc() {
  const b = createSeqViewSpec({ title: 'Dynamic nodes' });
  const genome = b.data.source({ ref: 'genome', stream: 'genome', data_kind: 'letters' });
  const chrT = b.data.complex({ ref: 'chrT', source_ref: genome, label: 'chrT' });
  b.data.source({ ref: 'sig', stream: 'signal', data_kind: 'values', coordinate_system: 'chrT' });
  const signal = b.data.annotation({ ref: 'signal', geometry: 'values', source_ref: 'sig' });
  const v = b.view({ ref: 'dyn', label: 'Dynamic nodes' });
  v.block({ complex: chrT });
  v.track({ label: '' }).ruler();
  v.track({ label: 'Signal' }).values({ data: signal });
  v.track({ label: "Genes (the manifest's)" }).ranges({ data: 'sig/genes' });
  for (const k of Object.keys(PEAKS))
    v.track({ label: `Peak ${k} (a patch's)` }).ranges({ data: `sig/${k}`, labels: true });
  return b.getDocument();
}

/**
 * Genes keyed G1 to G5 come with the `genes` source's manifest; each names a nested source of its transcripts, and each
 * transcript one of its domains (two deep). Three templates (`for_each` over `genes` at depth 0, 1 and 2) give a track
 * per gene, per transcript and per domain; a nested source is asked for only because a template reaches its nodes, and
 * not past `max_depth`.
 */
const GENES = [200_000, 520_000, 900_000, 1_300_000, 1_650_000];
export const items = (key: string, label: string, at: number, len: number): StreamNode => ({
  kind: 'annotation',
  key,
  params: { geometry: 'ranges', target: 'chrT/chrT', label, value: { start: [at], stop: [at + len] } },
});
const source = (key: string, stream: string): StreamNode => ({
  kind: 'source',
  key,
  params: { stream, data_kind: 'items' },
});
/** Every manifest by stream id: the genes', each gene's transcripts', each transcript's domains'. */
function keyed(): Record<string, StreamNode[]> {
  const out: Record<string, StreamNode[]> = { genes: [] };
  GENES.forEach((at, g) => {
    const G = `G${g + 1}`;
    out.genes.push(items(G, `Gene ${G}`, at, 120_000), source(`${G}/tx`, `tx-${G}`));
    out[`tx-${G}`] = [1, 2].flatMap((k) => {
      const T = `${G}.${k}`;
      out[`dom-${T}`] = [items('d1', `${T} domain`, at + 20_000 * k, 15_000)];
      return [items(T, `Transcript ${T}`, at + 10_000 * k, 90_000), source(`${T}/dom`, `dom-${T}`)];
    });
  });
  return out;
}
/** The keyed sources over the same genome; `asked`: the manifests asked for, in order (nested ones after a delay). */
export function keyedLoaders(asked: string[]): Loaders {
  const manifests = keyed();
  const genome = scripted(lettersManifest(DYNAMIC_CHR), () => new Error('no letters here'));
  const of = (id: string): typeof genome.provider => ({
    fetch: (r) =>
      new Promise((ok) => setTimeout(ok, id === 'genes' ? 0 : dynamicKnobs.latency)).then(() => {
        asked.push(id);
        const m: Manifest = {
          systems: [{ key: 'chrT', extent: [0, DYNAMIC_CHR[0][1]], numbering: null, label: null }],
          levels: [],
        };
        return { ...emptyAt(r), kind: 'manifest', manifest: { ...m, nodes: manifests[id] } } as DataPatch;
      }),
  });
  return {
    load: () => Promise.reject(new Error('no files')),
    stream: (src) =>
      src.stream === 'genome' ? genome.provider : src.stream && src.stream in manifests ? of(src.stream) : null,
  };
}
const emptyAt = (r: StreamRequest) => ({
  source: r.source,
  system: null,
  start: 0,
  stop: 0,
  row_start: 0,
  row_stop: 0,
  level: -1,
});

export function keyedDoc(maxDepth: number) {
  const b = createSeqViewSpec({ title: 'Templates' });
  const genome = b.data.source({ ref: 'genome', stream: 'genome', data_kind: 'letters' });
  const chrT = b.data.complex({ ref: 'chrT', source_ref: genome, label: 'chrT' });
  b.data.source({ ref: 'genes', stream: 'genes', data_kind: 'items', coordinate_system: 'chrT', max_depth: maxDepth });
  const v = b.view({ ref: 'tpl', label: 'Templates' });
  v.block({ complex: chrT });
  v.track({ label: '' }).ruler();
  // No track label: each instance is labeled by its node.
  for (const depth of [0, 1, 2])
    v.track({ for_each: { kind: 'annotation', source: 'genes', depth } }).ranges({ data: '$node', labels: true });
  return b.getDocument();
}

// ---- Seq*/Streamed data ----
export const STREAMED_SECTION = 'streamed.main';
/** Clear of the fixture's N runs (each 100 kb's first 1 kb). */
export const STREAMED_AT = 12_345_100;
const STREAMED_CHR: [string, number][] = [['chrS', 20_000_000]];

/**
 * The provider, its answers delayed by the knob's latency and failing at the knob's share (both read when asked, so
 * the loaders stay one object).
 */
export const streamedKnobs = { latency: 250, failure: 0 };
export const late = <T>(answer: () => Promise<T> | T): Promise<T> =>
  new Promise((ok) => setTimeout(ok, streamedKnobs.latency)).then(() =>
    Math.random() < streamedKnobs.failure ? Promise.reject(new Error("503 from the story's provider")) : answer(),
  );
export function streamedLoaders(asked?: (source: string, r: StreamRequest) => void): Loaders {
  const keys = ['chrS'];
  const sources = {
    genome: scripted(lettersManifest(STREAMED_CHR), (r: StreamRequest): DataPatch => sampledLetters(r, keys)),
    genes: scripted(itemsManifest(STREAMED_CHR, GENE_LEVELS, VOCABULARY.genes), (r) =>
      sampledItems(r, keys, 'genes', GENE_LEVELS),
    ),
    variants: scripted(itemsManifest(STREAMED_CHR, VARIANT_LEVELS, VOCABULARY.variants), (r) =>
      sampledItems(r, keys, 'variants', VARIANT_LEVELS),
    ),
    signal: scripted(valuesManifest(STREAMED_CHR, LEVELS), (r) => sampledValues(r, keys)),
  };
  const slow = (name: string, inner: (typeof sources)['genome']['provider']): typeof inner => ({
    fetch: (r, sig) => {
      // The request as made (`asked`: the prefetch story's log), before the latency.
      if (r.system !== null) asked?.(name, r);
      return r.system === null ? inner.fetch(r, sig) : late(() => inner.fetch(r, sig));
    },
  });
  const by = Object.fromEntries(Object.entries(sources).map(([k, x]) => [k, slow(k, x.provider)]));
  return {
    load: () => Promise.reject(new Error('no files')),
    stream: (src) => (src.stream ? (by[src.stream] ?? null) : null),
  };
}

/**
 * The streamed chromosome's document. With `snapshots` (title, start, stop), a focus node at the first and a snapshot
 * for each setting it (the Prefetch story: entering one is drawn from the cache when its neighbors were hinted).
 */
export function streamedDoc(o: { snapshots?: readonly (readonly [string, number, number])[] } = {}) {
  const b = createSeqViewSpec({ title: 'Streamed data' });
  const genome = b.data.source({ ref: 'genome', stream: 'genome', data_kind: 'letters' });
  const chr = b.data.complex({ ref: 'chrS', source_ref: genome, label: 'chrS' });
  const v = b.view({ ref: 'streamed', label: 'Streamed data' });
  v.block({ complex: chr });
  const [first] = o.snapshots ?? [];
  if (first) v.focus({ ref: 'streamed.focus', target: 'chrS/chrS', start: first[1], stop: first[2] });
  for (const [title, start, stop] of o.snapshots ?? [])
    b.snapshot({
      ref: title,
      title,
      duration_ms: 2000,
      transition_ms: 0,
      ops: [{ op: 'set', ref: 'streamed.focus', params: { start, stop } }],
    });
  v.track({ label: '' }).ruler();
  v.track({ label: 'Sequence' }).letters({ color: b.colors.scheme('nucleotide') });
  for (const k of ['genes', 'variants']) {
    b.data.source({ ref: `${k}_src`, stream: k, data_kind: 'items', coordinate_system: 'chrS' });
    b.data.annotation({ ref: k, geometry: 'ranges', source_ref: `${k}_src` });
  }
  b.data.source({ ref: 'signal_src', stream: 'signal', data_kind: 'values', coordinate_system: 'chrS' });
  const signal = b.data.annotation({ ref: 'signal', geometry: 'values', source_ref: 'signal_src' });
  v.track({ label: 'Genes', height: 4 }).ranges({
    data: 'genes',
    labels: true,
    color: { theme: 'field_palette', params: { field: 'category' } },
  });
  v.track({ label: 'Variants' }).variants({
    data: 'variants',
    color: { theme: 'letter', params: { colors: 'Nucleotide' } },
  });
  v.track({ label: 'Signal' }).values({ data: signal });
  return b.getDocument();
}

// ---- Seq*/Streamed data: a protein on its genome (ST7) ----
export const CODING_PROTEIN = 'coding.protein';
/** The protein's section beside chr2's at its gene (`codingDoc`: core/_test/format/stream-provider.ts). */
export const codingStoryDoc = () => codingDoc({ genome: true });
/**
 * Residue i of the protein hovered (the first row whose hit lands on it, as the pointer would) and its label computed:
 * its codon line (which asks for the genome's letters until they are held); null: not hovered.
 */
export function hoverResidue(ctx: SeqViewerContext, i: number): string | null {
  const l = ctx.layout.get(CODING_PROTEIN)?.current;
  if (!l) return null;
  const h = Hit.create();
  let y = 1;
  for (; y < l.height; y += 2) {
    ctx.layout.hitTest(CODING_PROTEIN, l.dataLeft + 1, y, h);
    if (h.block === 0 && h.index >= 0) break;
  }
  for (let x = l.dataLeft; x < l.width; x++) {
    ctx.layout.hitTest(CODING_PROTEIN, x, y, h);
    if (h.block !== 0 || h.index !== i) continue;
    ctx.hover.setFromHit(h, 'pointer');
    return ctx.labels.compute(ctx.hover.live.current)?.lines.find((x) => x.key === 'codon')?.value ?? null;
  }
  return null;
}
/** The synthetic genome's first three chromosomes: genes, variants and phyloP, answered after the knobs' latency. */
export function codingLoaders(): Loaders {
  // Letters too: a residue's codon in its label reads them.
  const s = hg38Streams({ chromosomes: GRCH38.slice(0, 3), letters: true });
  return {
    load: () => Promise.reject(new Error('no files')),
    stream: (src) => {
      const p = s.stream(src);
      return p && { fetch: (r, sig) => (r.system === null ? p.fetch(r, sig) : late(() => p.fetch(r, sig))) };
    },
  };
}

// ---- Seq*/Streamed data: a protein on its streamed transcript (ST7) ----
export const TRANSCRIPT_GENOME = 'tx.genome';
export const TRANSCRIPT_SECTION = 'tx.item';
/** The genome section and the transcript's own (`transcriptDoc`: core/_test/format/stream-provider.ts). */
export const transcriptStoryDoc = () => transcriptDoc({ item: true });
/** As `codingLoaders`, the genes with the transcript `T1` and its `map` (a transcripts file's rows). */
export function transcriptLoaders(): Loaders {
  const three = GRCH38.slice(0, 3);
  const s = hg38Streams({ chromosomes: three, letters: true });
  const keys = three.map((c) => c[0]);
  const genes = scripted(itemsManifest(three, GENE_LEVELS, VOCABULARY.genes), (r) => transcriptItems(r, keys));
  return {
    load: () => Promise.reject(new Error('no files')),
    stream: (src) => {
      const p = src.stream === 'genes' ? genes.provider : s.stream(src);
      return p && { fetch: (r, sig) => (r.system === null ? p.fetch(r, sig) : late(() => p.fetch(r, sig))) };
    },
  };
}

// ---- Seq*/Performance: the glyph atlas's lazy alignment ----
const ATLAS_URI = 'https://example.org/atlas-lazy';
export const ATLAS_ROWS = 100_000;
export const ATLAS_COLUMNS = 5_000;
const ATLAS_LETTERS = 'ACDEFGHIKLMNPQRSTVWY';
const encoder = new TextEncoder();

/** Tiles of letters (and their share) from a hash of the cell; rows and labels by index. Nothing is kept. */
const atlasTransport: TileTransport = async (path) => {
  const [kind, ...rest] = path.split('/');
  const a = rest.map(Number);
  if (kind === 'tiles') {
    const out = new Uint8Array(2 * TILE * TILE);
    for (let i = 0; i < TILE * TILE; i++) {
      const h = Math.imul(i * 2654435761 + a[2] * 97 + a[3] * 31 + a[0], 40503) >>> 0;
      out[i] = ATLAS_LETTERS.charCodeAt(h % 20);
      out[TILE * TILE + i] = 120 + (h % 120);
    }
    return out;
  }
  if (kind === 'labels')
    return encoder.encode(JSON.stringify(Array.from({ length: a[1] - a[0] }, (_, i) => `SEQ${a[0] + i}`)));
  if (kind === 'rows') {
    let row = '';
    for (let c = 0; c < ATLAS_COLUMNS; c++) {
      const h = Math.imul((a[0] + 1) * 73856093 + c * 19349663, 2654435761) >>> 0;
      row += h % 6 === 0 ? '-' : ATLAS_LETTERS[(h >>> 8) % 20];
    }
    return encoder.encode(JSON.stringify({ label: `SEQ${a[0]}`, row }));
  }
  throw new Error(`404 ${path}`);
};
export const atlasLoaders: Loaders = {
  load: () => Promise.reject(new Error('the atlas story loads nothing')),
  stream: tilesStream((u) => (u === ATLAS_URI ? atlasTransport : null)),
};

export function atlasDoc(): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'A lazy alignment of 100,000 × 5,000' });
  const aln = b.data.alignment({
    ref: 'lazy',
    row_count: ATLAS_ROWS,
    column_count: ATLAS_COLUMNS,
    source: { uri: ATLAS_URI, protocol: 'seqstar-tiles/1' },
  } as never);
  const v = b.view({ ref: 'g', label: 'Lazy', height: 'box' });
  v.block({ alignment: aln });
  v.track({ label: '' }).ruler();
  v.track({ label: 'Rows' }).rows();
  return b.getDocument();
}
