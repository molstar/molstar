/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { GffRow } from './gff-parse.js';

// GFF3's and GTF's hierarchy assembled by the line's type (_spec/README.md "Data", GFF3 and GTF; EC-227): a pure
// function over parsed rows, so a tile's rows read as the whole file's. A line is a transcript when its type is one of
// `transcriptTypes`, a piece when its type is a piece type, a gene when a transcript names it as its parent. A piece
// whose parent is of another type (a prokaryote's gene → CDS) makes that parent a transcript (counted by type). A child
// naming several parents is a piece of each. A `Parent` chain deeper than `MaxDepth` or a cycle leaves its rows as items
// of their own; every other feature is an item of its own by `ID` (rows sharing an `ID` its pieces).

/**
 * The Sequence Ontology's transcript types written in Ensembl, GENCODE, NCBI and Prokka files; any type ending in
 * `_transcript` is one too.
 */
export const DefaultTranscriptTypes: readonly string[] = [
  'mRNA',
  'transcript',
  'primary_transcript',
  'lnc_RNA',
  'ncRNA',
  'rRNA',
  'tRNA',
  'snRNA',
  'snoRNA',
  'miRNA',
  'scRNA',
  'misc_RNA',
  'pseudogenic_transcript',
  'unconfirmed_transcript',
  'V_gene_segment',
  'D_gene_segment',
  'J_gene_segment',
  'C_gene_segment',
];

/** A transcript's pieces' types (GTF's `start_codon` and `stop_codon` widen the CDS hull only). */
export const PieceTypes: readonly string[] = [
  'exon',
  'CDS',
  'five_prime_UTR',
  'three_prime_UTR',
  'UTR',
  'start_codon',
  'stop_codon',
];

/** The deepest `Parent` chain read as a hierarchy. */
export const MaxDepth = 16;

/** An item: a transcript, a gene (its transcripts merged) or a feature of its own. Positions 0-based, end exclusive. */
export interface GffItem {
  readonly id: string;
  readonly seqid: string;
  readonly type: string;
  readonly role: 'transcript' | 'gene' | 'feature';
  /** Its extent: its own lines' (a synthesized transcript's or gene's: the hull of what it holds). */
  readonly start: number;
  readonly end: number;
  readonly strand: -1 | 0 | 1;
  /** Its pieces, `[start, end, start, end, …]`, sorted and apart. */
  readonly pieces: readonly number[];
  /** The hull of its CDS (and GTF's codons); absent without one. */
  readonly thickStart?: number;
  readonly thickEnd?: number;
  /** The phase of its first CDS piece in its direction. */
  readonly cdsPhase?: number;
  /** Its gene's `ID`, and the gene's `Name` (GTF's `gene_name`). */
  readonly gene?: string;
  readonly geneName?: string;
  /** Its line's attributes (a synthesized one's: its first child's). */
  readonly attributes: Readonly<Record<string, string>>;
  /** Its line's multi-valued attributes (GFF3's `Alias`, …), each value decoded after the split. */
  readonly lists?: Readonly<Record<string, readonly string[]>>;
  /** Its first line. */
  readonly line: number;
  /** Made from its children: no line of its own (GTF 2.2's transcripts, a missing parent). */
  readonly synthetic?: boolean;
}

export interface AssembleOptions {
  /** The types read as transcripts (`DefaultTranscriptTypes`; `*_transcript` always). */
  readonly transcriptTypes?: readonly string[];
}

/** What assembling found besides the items: what a runtime says once. */
export interface Assembly {
  readonly transcripts: GffItem[];
  /** The genes of transcripts, each its transcripts' pieces merged. */
  readonly genes: GffItem[];
  /** Every other feature, an item of its own. */
  readonly features: GffItem[];
  /** Types outside `transcriptTypes` read as transcripts (a piece's parent), with their counts. */
  readonly promoted: ReadonlyMap<string, number>;
  /** Lines left as items of their own: in a cycle, deeper than `MaxDepth`, or naming a parent no line has. */
  readonly cycles: readonly number[];
  readonly deep: readonly number[];
  readonly orphans: readonly number[];
  /** Transcripts made from their pieces (no line of their own). */
  readonly synthetic: number;
  /** Children outside their parent's extent (lines). */
  readonly outside: readonly number[];
}

/** Whether `type` is a transcript's by the list (and `*_transcript`). */
export function transcriptTest(types: readonly string[] = DefaultTranscriptTypes): (type: string) => boolean {
  const set = new Set(types);
  return (type) => set.has(type) || type.endsWith('_transcript');
}

/** Intervals `[s, e, …]` merged where they overlap or touch, sorted. */
export function mergeIntervals(list: number[]): number[] {
  const n = list.length / 2;
  if (n <= 1) return list;
  const order = Array.from({ length: n }, (_, i) => i).sort((a, b) => list[2 * a] - list[2 * b]);
  const out: number[] = [];
  for (const i of order) {
    const s = list[2 * i];
    const e = list[2 * i + 1];
    if (out.length && s <= out[out.length - 1]) out[out.length - 1] = Math.max(out[out.length - 1], e);
    else out.push(s, e);
  }
  return out;
}

const isPiece = new Set(PieceTypes);
const isCds = (t: string) => t === 'CDS' || t === 'start_codon' || t === 'stop_codon';

interface Feature {
  readonly id: string;
  readonly rows: number[];
  start: number;
  end: number;
}

/** A transcript's pieces and thick part from its children's rows. */
function piecesOf(rows: readonly GffRow[], children: readonly number[], strand: number) {
  const exons: number[] = [];
  const parts: number[] = [];
  let thickStart = Number.POSITIVE_INFINITY;
  let thickEnd = Number.NEGATIVE_INFINITY;
  let phaseAt = -1;
  for (const c of children) {
    const r = rows[c];
    if (r.type === 'exon') exons.push(r.start, r.end);
    else if (r.type !== 'start_codon' && r.type !== 'stop_codon') parts.push(r.start, r.end);
    if (isCds(r.type)) {
      thickStart = Math.min(thickStart, r.start);
      thickEnd = Math.max(thickEnd, r.end);
    }
    // the first CDS in the transcript's direction: the lowest start forward, the highest end reversed
    if (r.type === 'CDS' && (phaseAt < 0 || (strand < 0 ? r.end > rows[phaseAt].end : r.start < rows[phaseAt].start)))
      phaseAt = c;
  }
  const pieces = mergeIntervals(exons.length ? exons : parts);
  return {
    pieces,
    ...(Number.isFinite(thickStart) ? { thickStart, thickEnd } : {}),
    ...(phaseAt >= 0 && rows[phaseAt].phase >= 0 ? { cdsPhase: rows[phaseAt].phase } : {}),
  };
}

/**
 * The hierarchy of `rows` (one sequence's or several: ids are the file's): transcripts with their pieces, genes, the
 * other features, and what was read past (see `Assembly`). Never throws for a file's content.
 */
export function assembleTranscripts(rows: readonly GffRow[], options: AssembleOptions = {}): Assembly {
  const isTranscript = transcriptTest(options.transcriptTypes);
  // features by id: rows sharing an `ID` are one (discontinuous) feature
  const byId = new Map<string, Feature>();
  for (let i = 0; i < rows.length; i++) {
    const id = rows[i].id;
    if (id === undefined) continue;
    let f = byId.get(id);
    if (!f) byId.set(id, (f = { id, rows: [], start: rows[i].start, end: rows[i].end }));
    f.rows.push(i);
    f.start = Math.min(f.start, rows[i].start);
    f.end = Math.max(f.end, rows[i].end);
  }
  const parentsOf = (id: string) => rows[byId.get(id)!.rows[0]].parents;

  // the chains: a feature in a cycle or deeper than MaxDepth is flat (its rows items of their own)
  const depth = new Map<string, number>();
  const flat = new Set<string>();
  const cycles: number[] = [];
  const deep: number[] = [];
  // a feature's depth by an explicit stack (a chain of 20,000 never overflows the call stack), saturated past
  // `MaxDepth` (deeper is flat whatever its depth); −1 in or under a cycle
  const onPath = new Set<string>();
  const stack: { id: string; next: number; d: number }[] = [];
  const visit = (root: string): number => {
    const known = depth.get(root);
    if (known !== undefined) return known;
    stack.push({ id: root, next: 0, d: 0 });
    onPath.add(root);
    while (stack.length) {
      const f = stack[stack.length - 1];
      const ps = parentsOf(f.id);
      if (f.d >= 0 && f.next < ps.length) {
        const p = ps[f.next++];
        if (!byId.has(p)) continue;
        const pd = depth.get(p);
        if (pd !== undefined) f.d = pd < 0 ? -1 : Math.max(f.d, Math.min(pd + 1, MaxDepth + 1));
        else if (onPath.has(p)) f.d = -1;
        else {
          stack.push({ id: p, next: 0, d: 0 });
          onPath.add(p);
        }
        continue;
      }
      stack.pop();
      onPath.delete(f.id);
      depth.set(f.id, f.d);
      const child = stack[stack.length - 1];
      if (child) child.d = f.d < 0 || child.d < 0 ? -1 : Math.max(child.d, Math.min(f.d + 1, MaxDepth + 1));
    }
    return depth.get(root)!;
  };
  for (const id of byId.keys()) {
    const d = visit(id);
    if (d < 0) {
      flat.add(id);
      for (const r of byId.get(id)!.rows) cycles.push(rows[r].line);
    } else if (d > MaxDepth) {
      flat.add(id);
      for (const r of byId.get(id)!.rows) deep.push(rows[r].line);
    }
  }
  // a parent read as one: it exists and is not flat
  const live = (p: string) => byId.has(p) && !flat.has(p);

  // pieces by their transcript (a child naming several parents a piece of each); parents missing are made from them
  const children = new Map<string, number[]>();
  const missing = new Map<string, number[]>();
  const orphans: number[] = [];
  const promoted = new Map<string, number>();
  const transcriptIds = new Set<string>();
  const used = new Set<number>();
  // the live parents of transcripts are genes: a piece naming one directly never makes it a transcript (its mRNA
  // keeps its gene); the piece is a feature of its own
  const genesOf = new Set<string>();
  for (const [id, f] of byId) {
    if (flat.has(id) || !isTranscript(rows[f.rows[0]].type)) continue;
    transcriptIds.add(id);
    for (const p of rows[f.rows[0]].parents) if (live(p)) genesOf.add(p);
  }
  for (let i = 0; i < rows.length; i++) {
    const r = rows[i];
    if (!isPiece.has(r.type) || (r.id !== undefined && flat.has(r.id))) continue;
    let any = false;
    for (const p of r.parents) {
      if (live(p)) {
        if (genesOf.has(p) && !transcriptIds.has(p)) continue;
        if (!transcriptIds.has(p)) {
          // a prokaryote's gene → CDS: the parent a transcript-level item
          const t = rows[byId.get(p)!.rows[0]].type;
          if (isTranscript(t)) continue;
          transcriptIds.add(p);
          promoted.set(t, (promoted.get(t) ?? 0) + 1);
        }
        let c = children.get(p);
        if (!c) children.set(p, (c = []));
        c.push(i);
        any = true;
      } else if (!byId.has(p)) {
        orphans.push(r.line);
        let c = missing.get(p);
        if (!c) missing.set(p, (c = []));
        c.push(i);
        any = true;
      }
    }
    if (any) used.add(i);
  }
  // a transcript with no line (GTF 2.2, a missing parent): made from its pieces, its extent their hull
  let synthetic = 0;
  const made = new Map<string, GffItem>();
  for (const [id, list] of missing) {
    const first = rows[list[0]];
    let start = Number.POSITIVE_INFINITY;
    let end = Number.NEGATIVE_INFINITY;
    for (const c of list) {
      start = Math.min(start, rows[c].start);
      end = Math.max(end, rows[c].end);
    }
    synthetic++;
    // GTF's piece lines carry the gene (`gene_id`); a GFF3 piece names only its missing parent
    const gene = first.attributes.gene_id;
    made.set(id, {
      id,
      seqid: first.seqid,
      type: 'transcript',
      role: 'transcript',
      start,
      end,
      strand: first.strand,
      ...piecesOf(rows, list, first.strand),
      ...(gene ? { gene } : {}),
      attributes: first.attributes,
      line: first.line,
      synthetic: true,
    });
  }

  const outside: number[] = [];
  const transcripts: GffItem[] = [];
  for (const id of transcriptIds) {
    const f = byId.get(id)!;
    const head = rows[f.rows[0]];
    const kids = children.get(id) ?? [];
    for (const c of kids) if (rows[c].start < f.start || rows[c].end > f.end) outside.push(rows[c].line);
    const got = kids.length ? piecesOf(rows, kids, head.strand) : { pieces: [f.start, f.end] };
    const gene = head.parents.find((p) => live(p) && !transcriptIds.has(p)) ?? head.parents.find((p) => !byId.has(p));
    if (gene !== undefined && !byId.has(gene)) orphans.push(head.line);
    transcripts.push({
      id,
      seqid: head.seqid,
      type: head.type,
      role: 'transcript',
      start: f.start,
      end: f.end,
      strand: head.strand,
      ...got,
      ...(gene ? { gene } : {}),
      attributes: head.attributes,
      ...(head.lists ? { lists: head.lists } : {}),
      line: head.line,
    });
    for (const r of f.rows) used.add(r);
  }
  for (const t of made.values()) transcripts.push(t);
  transcripts.sort((a, b) => a.line - b.line);

  // genes: the parents of transcripts; a gene no line has (GTF without gene lines) made from its transcripts
  const byGene = new Map<string, GffItem[]>();
  for (const t of transcripts) {
    if (t.gene === undefined) continue;
    let g = byGene.get(t.gene);
    if (!g) byGene.set(t.gene, (g = []));
    g.push(t);
  }
  const genes: GffItem[] = [];
  const named = new Map<string, string>();
  for (const [id, ts] of byGene) {
    const f = byId.get(id);
    const head = f ? rows[f.rows[0]] : undefined;
    if (f) for (const r of f.rows) used.add(r);
    const all: number[] = [];
    let ts0 = Number.POSITIVE_INFINITY;
    let te0 = Number.NEGATIVE_INFINITY;
    let start = Number.POSITIVE_INFINITY;
    let end = Number.NEGATIVE_INFINITY;
    for (const t of ts) {
      // no spread: a transcript of 200,000 exons passes no 400,000 arguments
      for (let k = 0; k < t.pieces.length; k++) all.push(t.pieces[k]);
      if (t.thickStart !== undefined) {
        ts0 = Math.min(ts0, t.thickStart);
        te0 = Math.max(te0, t.thickEnd!);
      }
      start = Math.min(start, t.start);
      end = Math.max(end, t.end);
    }
    const attributes = head?.attributes ?? ts[0].attributes;
    const name = attributes.Name ?? attributes.gene_name;
    if (name !== undefined) named.set(id, name);
    genes.push({
      id,
      seqid: head?.seqid ?? ts[0].seqid,
      type: head?.type ?? 'gene',
      role: 'gene',
      start: f ? f.start : start,
      end: f ? f.end : end,
      strand: head?.strand ?? ts[0].strand,
      pieces: mergeIntervals(all),
      ...(Number.isFinite(ts0) ? { thickStart: ts0, thickEnd: te0 } : {}),
      gene: id,
      ...(name !== undefined ? { geneName: name } : {}),
      attributes,
      ...(head?.lists ? { lists: head.lists } : {}),
      line: head?.line ?? ts[0].line,
      ...(f ? {} : { synthetic: true }),
    });
  }
  genes.sort((a, b) => a.line - b.line);
  // each transcript its gene's name (GTF writes it on the transcript's own line too)
  for (let k = 0; k < transcripts.length; k++) {
    const t = transcripts[k];
    const name = (t.gene !== undefined ? named.get(t.gene) : undefined) ?? t.attributes.gene_name;
    if (name !== undefined) transcripts[k] = { ...t, geneName: name };
  }

  // every other feature an item of its own: by `ID` (rows sharing it its pieces), else the line alone
  const features: GffItem[] = [];
  const seenKeys = new Map<string, number>();
  const done = new Set<string>();
  for (let i = 0; i < rows.length; i++) {
    if (used.has(i)) continue;
    const r = rows[i];
    if (r.id !== undefined && !flat.has(r.id) && byId.has(r.id)) {
      if (done.has(r.id)) continue;
      done.add(r.id);
      const f = byId.get(r.id)!;
      const pieces: number[] = [];
      for (const k of f.rows) pieces.push(rows[k].start, rows[k].end);
      for (const p of r.parents) if (!byId.has(p)) orphans.push(r.line);
      features.push(featureItem(r, r.id, f.start, f.end, mergeIntervals(pieces)));
      continue;
    }
    // a line without an id (or flat): its identity `type:start-end`, `#2` on a later identical one
    const base = r.id ?? `${r.type}:${r.start + 1}-${r.end}`;
    const k = (seenKeys.get(base) ?? 0) + 1;
    seenKeys.set(base, k);
    features.push(featureItem(r, k === 1 ? base : `${base}#${k}`, r.start, r.end, [r.start, r.end]));
  }
  return { transcripts, genes, features, promoted, cycles, deep, orphans, synthetic, outside };
}

function featureItem(r: GffRow, id: string, start: number, end: number, pieces: number[]): GffItem {
  return {
    id,
    seqid: r.seqid,
    type: r.type,
    role: 'feature',
    start,
    end,
    strand: r.strand,
    pieces,
    attributes: r.attributes,
    ...(r.lists ? { lists: r.lists } : {}),
    line: r.line,
  };
}
