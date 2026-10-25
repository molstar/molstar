/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Mappings and paths (seqspec:model#coordinate-maps). A mapping pairs elements of one coordinate system with another's:
// built once into a coordinate map from its method's own input (coordinates.ts), with `forward`, `status` and a codon's
// `spans` per from element and nothing per to position (G1-1). A row and its alignment's columns are joined by an
// implicit edge; paths are found breadth-first under the copy rule and composed into one coordinate map when resolved,
// never on hover.

import {
  CustomData,
  Issue,
  type MappingId,
  type NameId,
  Provenance,
  type ProvenanceInput,
  type RowIndex,
  type SeqIndex,
  type UnitId,
} from './common.js';
import { type CoordinateMap, type CoordinateSystem, CoordinateSystemKind, Coords } from './coordinates.js';
import { GeneticCode } from './genetic-code.js';
import { NO_LABEL, Numbering, type PositionLabel } from './numbering.js';
import { KeyKind, type KeyTarget, type SeqSet } from './set.js';

export type MappingMethod = 'identity' | 'segments' | 'pairwise' | 'by_numbering' | 'codon';
export const MappingStatus = { exact: 0, conflict: 1, none: 2 } as const;

export interface Mapping {
  readonly id: MappingId;
  readonly key: NameId;
  readonly method: MappingMethod;
  readonly from: CoordinateSystem;
  readonly to: CoordinateSystem;
  /** Built once from the method's input; to → from reads it (`Coords.forEachFromAt`, `forEachFrom`). */
  readonly map: CoordinateMap;
  /** From index → to index, -1: none; a span mapping: the first of the span. One array read per lookup. */
  readonly forward: Int32Array;
  /** Per from position: MappingStatus. */
  readonly status: Uint8Array;
  /** To positions per from position: 3 for `codon`, else 1. */
  readonly span: number;
  /** span > 1: per from position, its `span` to positions in reading order, -1: none. */
  readonly spans: Int32Array | null;
  /** `codon` only (seqspec:model#codon-maps). */
  readonly codon: CodonInfo | null;
  /** Never null: from the method when the input gave none. */
  readonly provenance: Provenance;
  readonly custom: CustomData;
}

/** A nucleotide's phase is its sub-position modulo the span, read through the map. */
export interface CodonInfo {
  /** 1, or -1: the CDS reads the reverse complement. */
  readonly strand: number;
  /** NCBI's table id. */
  readonly geneticCode: number;
}

/** Endpoints by key; a row is `"<alignment>/<row>"`. `align` is the compiler's: it passes `pairwise`. */
export type MappingInput = { key: string; from: string; to: string; provenance?: ProvenanceInput; custom?: object } & (
  | { method: 'identity'; fromRange?: [number, number] }
  | { method: 'segments'; segments: [number, number, number | null][] }
  | { method: 'pairwise'; cigar: string }
  | { method: 'by_numbering'; numbering: string; toNumbering?: string }
  | CodonInput
);
/** A protein (`from`) and the DNA or RNA coding it (`to`): `cds` ranges in transcript order, or one from `start`. */
export type CodonInput = {
  method: 'codon';
  cds?: [number, number][];
  start?: number;
  fromRange?: [number, number];
  strand?: '+' | '-';
  geneticCode?: number;
};

/**
 * One step of a path: a mapping (`reversed` when walked to → from); else, with `mapping` -1, a row edge
 * (`alignment`, `row`: row → columns, reversed columns → row) or a unit step (`unit`: unit → its molecule
 * within its range, reversed molecule → unit). `unit` is an addition to seqspec:model#indices-and-labels: it names the
 * unit step.
 */
export interface PathStep {
  mapping: number;
  reversed: boolean;
  alignment: number;
  row: number;
  unit: number;
}
export interface ComposedPath {
  readonly from: CoordinateSystem;
  readonly to: CoordinateSystem;
  /**
   * The steps composed and clipped to both ends' domains. Its `from` is the path's, except through a codon step taken
   * from the nucleotides (the protein side is `from`) and where the path's from side is a source used twice (the copies'
   * side is `from`): then it is the path's `to`.
   */
  readonly map: CoordinateMap;
  /** Elements of the path's `from` system. */
  readonly length: number;
  /**
   * Sparse (ST7, seqspec:model#paths): what the path pairs on its `from` side, as sorted disjoint runs [start, stop)
   * (pairs), each with its first slot in `heldTo`. Nothing is held outside them, so a protein's path through a
   * chromosome holds its CDS, never the chromosome. Read through `Paths`.
   */
  readonly runs: Int32Array;
  readonly slots: Int32Array;
  /** Per held from element (by slot): its first to position (a residue: its codon's first), -1: none. */
  readonly heldTo: Int32Array;
  /** 3 through a codon step taken from the protein side, with `heldSpans`; else 1, null. */
  readonly span: number;
  /** span > 1: per held from element (by slot), its `span` to positions in reading order, -1: none. */
  readonly heldSpans: Int32Array | null;
  /** The codon step's MappingId; -1: none. */
  readonly codon: number;
}

/**
 * The row → bound coordinate systems table of one alignment (seqspec:model#locations-and-loci): per bound row, the
 * units its mappings reach (a molecule-level binding reaches every unit over the molecule), each with its mapping's
 * map, the row on its `from` side or its `to` side (callers check the unit's range). Built when the set is assembled.
 */
export interface RowBindings {
  /** Bound rows, increasing. */
  readonly rows: Int32Array;
  /** rows.length + 1: each row's entries. */
  readonly start: Int32Array;
  /** Per bound row: how many mappings bind it; more than one makes loci conversions ambiguous. */
  readonly bindings: Int32Array;
  /** Per entry. */
  readonly unit: Int32Array;
  readonly mapping: Int32Array;
  readonly map: readonly CoordinateMap[];
  /** 1: the row is the map's `from`, 0: its `to`. */
  readonly rowIsFrom: Uint8Array;
  /** Per entry: its row's index in `rows`. */
  readonly entryRow: Int32Array;
  /** The unit index: the units reached, increasing; unit k's entries are byUnit[unitStart[k]..unitStart[k + 1]). */
  readonly units: Int32Array;
  readonly unitStart: Int32Array;
  readonly byUnit: Int32Array;
}

/** What the set keeps per mapping to rebuild it when an endpoint changes. */
interface Built extends Mapping {
  readonly input: MappingInput;
  readonly deps: readonly unknown[];
}
/** The set's internal fields this module reads. */
interface Assembled {
  readonly bindings: ReadonlyMap<number, RowBindings>;
}

/** Counted steps; unit steps are free. */
const MAX_STEPS = 4;
/** Shortest paths compared for ambiguity, at most. */
const MAX_PATHS = 16;
const DEFAULT_METHOD = {
  identity: 'identity',
  segments: 'curated',
  pairwise: 'alignment',
  by_numbering: 'numbering',
  codon: 'translation',
};
/** By KeyKind (molecule, unit, alignment, row; set.ts's, written out: the modules import each other). */
const SPACE_OF_KEY = [
  CoordinateSystemKind.molecule,
  CoordinateSystemKind.unit,
  CoordinateSystemKind.columns,
  CoordinateSystemKind.row,
];
const target: KeyTarget = { kind: -1, id: -1, row: -1 };
const label: PositionLabel = { num: 0, suffix: -1 as NameId, name: -1 as NameId };

/** Elements in a coordinate system; -1: no such system (a removed id, a row not held). */
function systemLength(set: SeqSet, s: CoordinateSystem): number {
  if (s.kind === CoordinateSystemKind.molecule) return set.molecules[s.id]?.length ?? -1;
  if (s.kind === CoordinateSystemKind.unit)
    return set.units.alive[s.id] ? (set.molecules[set.units.molecule[s.id]]?.length ?? -1) : -1;
  const a = set.alignments[s.id];
  if (s.kind === CoordinateSystemKind.row) return a ? a.rows.lengthOf(s.row as RowIndex) : -1;
  return s.kind === CoordinateSystemKind.columns && a ? a.columnCount : -1;
}

/** A coordinate system's domain: a unit's range, else [0, length). */
function domainOf(set: SeqSet, s: CoordinateSystem): [number, number] {
  if (s.kind === CoordinateSystemKind.unit) return [set.units.rangeStart[s.id], set.units.rangeStop[s.id]];
  return [0, Math.max(0, systemLength(set, s))];
}

/** The coordinate system a key names: a molecule, a unit, a row (`"R/3"`) or an alignment's columns; false when none. */
function systemOf(set: SeqSet, key: string, out: CoordinateSystem): boolean {
  if (!set.byName(key, target) || SPACE_OF_KEY[target.kind] === undefined) return false;
  out.kind = SPACE_OF_KEY[target.kind];
  out.id = target.id;
  out.row = target.kind === KeyKind.row ? target.row : -1;
  return systemLength(set, out) >= 0;
}

const keyOf = (s: CoordinateSystem) => `${s.kind}:${s.id}:${s.row}`;

function moleculeOf(set: SeqSet, s: CoordinateSystem) {
  return set.molecules[s.kind === CoordinateSystemKind.unit ? set.units.molecule[s.id] : s.id] ?? null;
}
/**
 * A molecule's or unit's coordinate system over a cyclic molecule (seqspec:model#cyclic-systems); a row or columns
 * never is.
 */
const isCyclic = (set: SeqSet, s: CoordinateSystem) =>
  (s.kind === CoordinateSystemKind.molecule || s.kind === CoordinateSystemKind.unit) && !!moleculeOf(set, s)?.cyclic;

function lettersOf(set: SeqSet, s: CoordinateSystem): { code: Uint8Array; element: Uint8Array | null } {
  if (s.kind !== CoordinateSystemKind.row) {
    const m = moleculeOf(set, s);
    return { code: m?.code ?? new Uint8Array(), element: m?.element ?? null };
  }
  const rows = (set.alignments[s.id] as NonNullable<SeqSet['alignments'][number]>).rows;
  const slot = rows.slotOf(s.row as RowIndex);
  return { code: rows.residues.subarray(rows.residueStart[slot], rows.residueStart[slot + 1]), element: null };
}

/** A system's numbering by name: a unit's own wins over its molecule's; a row's are `ordinal` and `row`. */
function numberingOf(set: SeqSet, s: CoordinateSystem, name: NameId) {
  if (s.kind === CoordinateSystemKind.row)
    return set.alignments[s.id]?.numberingsOf(s.row as RowIndex)?.byName(name) ?? null;
  const own = s.kind === CoordinateSystemKind.unit ? set.units.numberings[s.id]?.byName(name) : null;
  return own ?? moleculeOf(set, s)?.numberings.byName(name) ?? null;
}

/** What a mapping's arrays depend on beside its endpoints: the molecule, a unit's numberings, the alignment. */
function depsOf(set: SeqSet, s: CoordinateSystem): unknown[] {
  if (s.kind === CoordinateSystemKind.row) return [set.alignments[s.id], null];
  return [moleculeOf(set, s), s.kind === CoordinateSystemKind.unit ? set.units.numberings[s.id] : null];
}

/** Both endpoints, or null when one is gone or not a sequence, ligand, unit or held row. */
function endsOf(set: SeqSet, input: MappingInput): [CoordinateSystem, CoordinateSystem] | null {
  const from: CoordinateSystem = { kind: -1, id: -1, row: -1 };
  const to: CoordinateSystem = { kind: -1, id: -1, row: -1 };
  const ok = (key: string, s: CoordinateSystem) => systemOf(set, key, s) && s.kind !== CoordinateSystemKind.columns;
  return ok(input.from, from) && ok(input.to, to) ? [from, to] : null;
}

const isT = (c: number) => c === 84 || c === 85;

/**
 * Per path `from` element its first to position (`own`: the map's `from` is the path's, a residue's codon's first, and
 * a span map's `spans`); else the first from element holding it (the map read to → from). Edges.
 */
function arraysOf(m: CoordinateMap, n: number, own: boolean): { forward: Int32Array; spans: Int32Array | null } {
  const forward = new Int32Array(Math.max(0, n)).fill(-1);
  const spans = own && m.span > 1 ? new Int32Array(m.span * n).fill(-1) : null;
  for (let s = 0; s < m.count; s++) {
    for (let k = 0; k < m.length[s]; k++) {
      const q = m.at[s] + k;
      const t = m.dir[s] > 0 ? m.toLo[s] + k : m.toLo[s] + m.length[s] - 1 - k;
      if (!own) {
        if (forward[t] < 0) forward[t] = (q / m.span) | 0;
      } else if (spans) {
        spans[q] = t;
        if (q % m.span === 0) forward[q / m.span] = t;
      } else if (forward[q] < 0) forward[q] = t;
    }
  }
  return { forward, spans };
}

/**
 * A composed map's arrays held over its `from` side's support only (ST7): its runs, merged where they touch, and per
 * held element its first to position (`own`: as `arraysOf`'s; else the first from element holding it). Resolve-time.
 */
function heldOf(m: CoordinateMap, n: number, own: boolean) {
  const ends: [number, number][] = [];
  for (let s = 0; s < m.count; s++) {
    const a = own ? Math.floor(m.at[s] / m.span) : m.toLo[s];
    const e = own ? Math.floor((m.at[s] + m.length[s] - 1) / m.span) + 1 : m.toLo[s] + m.length[s];
    if (Math.min(e, n) > Math.max(a, 0)) ends.push([Math.max(a, 0), Math.min(e, n)]);
  }
  ends.sort((x, y) => x[0] - y[0]);
  const list: number[] = [];
  for (const [a, e] of ends)
    if (list.length && a <= list[list.length - 1]) list[list.length - 1] = Math.max(list[list.length - 1], e);
    else list.push(a, e);
  const runs = Int32Array.from(list);
  const slots = new Int32Array(runs.length >> 1);
  let total = 0;
  for (let k = 0; k < slots.length; k++) {
    slots[k] = total;
    total += runs[2 * k + 1] - runs[2 * k];
  }
  const held = { runs, slots, heldTo: new Int32Array(total).fill(-1), heldSpans: null as Int32Array | null };
  const spans = own && m.span > 1 ? new Int32Array(m.span * total).fill(-1) : null;
  held.heldSpans = spans;
  for (let s = 0; s < m.count; s++) {
    for (let k = 0; k < m.length[s]; k++) {
      const q = m.at[s] + k;
      const t = m.dir[s] > 0 ? m.toLo[s] + k : m.toLo[s] + m.length[s] - 1 - k;
      const i = own ? (q / m.span) | 0 : t;
      const slot = i < n ? slotIn(runs, slots, i) : -1;
      if (slot < 0) continue;
      if (!own) {
        if (held.heldTo[slot] < 0) held.heldTo[slot] = (q / m.span) | 0;
      } else if (spans) {
        spans[m.span * slot + (q % m.span)] = t;
        if (q % m.span === 0) held.heldTo[slot] = t;
      } else if (held.heldTo[slot] < 0) held.heldTo[slot] = t;
    }
  }
  return held;
}

/** Element i's slot among `runs` (sorted, disjoint pairs), -1: not held. Allocation-free, a binary search. */
function slotIn(runs: Int32Array, slots: Int32Array, i: number): number {
  let lo = 0;
  let hi = slots.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (runs[2 * mid + 1] <= i) lo = mid + 1;
    else hi = mid;
  }
  return lo < slots.length && runs[2 * lo] <= i ? slots[lo] + i - runs[2 * lo] : -1;
}

/**
 * A codon mapping's segments: the CDS ranges read in transcript order (on `-` each from its stop down, complemented),
 * codon k paired with residue r0 + k, one segment per range; the statuses by the genetic code once the map gives the
 * spans. Null when `fail` reported why.
 */
function codons(
  set: SeqSet,
  input: CodonInput & { from: string; to: string },
  from: CoordinateSystem,
  to: CoordinateSystem,
  quads: number[],
  fail: (code: string, at: string, m: string) => null,
) {
  const [na, nb, int] = [systemLength(set, from), systemLength(set, to), Number.isInteger];
  const alphabetOf = (s: CoordinateSystem) =>
    s.kind === CoordinateSystemKind.row ? set.alignments[s.id]?.alphabet : moleculeOf(set, s)?.alphabet;
  const [fa, ta] = [alphabetOf(from), alphabetOf(to)];
  if (fa !== 'protein' || (ta !== 'dna' && ta !== 'rna')) {
    return fail(
      'model.mapping.codon',
      '',
      `"${input.from}" must be a protein and "${input.to}" DNA or RNA, not ${fa} and ${ta}.`,
    );
  }
  const [strand, gc] = [input.strand === '-' ? -1 : 1, input.geneticCode ?? 1];
  if (!GeneticCode.has(gc)) return fail('model.mapping.codon', '/geneticCode', `NCBI defines no genetic code ${gc}.`);
  const [r0, r1] = input.fromRange ?? [0, na];
  if (!(int(r0) && int(r1) && r0 >= 0 && r0 < r1 && r1 <= na)) {
    return fail('model.mapping.range', '/fromRange', `[${r0}, ${r1}) is not within "${input.from}" (${na} residues).`);
  }
  const [start, cyclic] = [input.start ?? 0, isCyclic(set, to)];
  // On a cyclic `to` a CDS crosses the origin once at most: as two ranges, the second from the origin.
  const last = start + 3 * (r1 - r0);
  let laps = 0;
  const across: [number, number][] = [
    [start, nb],
    [0, last - nb],
  ];
  const cds = input.cds ?? (cyclic && last > nb ? (strand > 0 ? across : across.reverse()) : [[start, last]]);
  let total = 0;
  for (let k = 0; k < cds.length; k++) {
    const [s, t] = cds[k];
    const at = input.cds ? `/cds/${k}` : '/start';
    if (!(int(s) && int(t) && s >= 0 && s < t && t <= nb)) {
      return fail('model.mapping.range', at, `[${s}, ${t}) is not within "${input.to}" (${nb} positions).`);
    }
    const back = k > 0 && (strand > 0 ? s < cds[k - 1][1] : t > cds[k - 1][0]);
    if (back && !(cyclic && laps++ === 0 && (strand > 0 ? t <= cds[0][0] : s >= cds[0][1]))) {
      const order = strand > 0 ? 'ascending' : 'descending (strand -)';
      return fail('model.mapping.codon', at, `range ${k} overlaps range ${k - 1} or is out of ${order} order.`);
    }
    total += t - s;
  }
  if (total % 3) return fail('model.mapping.codon', '/cds', `the CDS has ${total} nucleotides, not a multiple of 3.`);
  const paired = Math.min(total / 3, r1 - r0);
  // Codons past the residues (a stop codon) pair with nothing: the support ends at the last paired residue.
  let q = 3 * r0;
  const end = 3 * (r0 + paired);
  for (const [s, t] of cds) {
    const n = Math.min(t - s, end - q);
    if (n > 0) quads.push(q, n, strand > 0 ? s : t - n, strand);
    q += t - s;
  }
  const status = (spans: Int32Array) => {
    const [aa, nt] = [lettersOf(set, from).code, lettersOf(set, to).code];
    const base = (j: number) => (strand > 0 ? nt[j] : GeneticCode.complement(nt[j]));
    const out = new Uint8Array(na).fill(MappingStatus.none);
    for (let r = r0; r < r0 + paired; r++) {
      const [x, y, z, res] = [base(spans[3 * r]), base(spans[3 * r + 1]), base(spans[3 * r + 2]), aa[r]];
      const exact =
        GeneticCode.translate(gc, x, y, z) === res ||
        (r === r0 && res === 77 && GeneticCode.isStart(gc, x, y, z)) ||
        (res === 85 && isT(x) && y === 71 && z === 65) || // Sec from TGA
        (res === 79 && isT(x) && y === 65 && z === 71); // Pyl from TAG
      out[r] = exact ? MappingStatus.exact : MappingStatus.conflict;
    }
    return out;
  };
  return { status, info: { strand, geneticCode: gc } as CodonInfo };
}

/** A mapping from its input over the set's current endpoints, or null with errors (seqspec:format#mappings). */
function build(set: SeqSet, input: MappingInput, id: number, issues: Issue[], path: string): Built | null {
  const err = (code: string, at: string, m: string) => {
    issues.push(Issue.error(code, `${path}${at}`, `Mapping "${input.key}": ${m}`));
    return null;
  };
  const ends = endsOf(set, input);
  if (!ends) {
    const what = 'a sequence, a ligand, a unit or a held alignment row';
    return err('model.mapping.endpoint', '', `"${input.from}" and "${input.to}" must each name ${what}.`);
  }
  const [from, to] = ends;
  const na = systemLength(set, from);
  const nb = systemLength(set, to);
  /** Segments as (at, length, toLo, dir): the method's own runs, never a per-position array. */
  const quads: number[] = [];
  const int = Number.isInteger;
  const sides = `both sides (${na} and ${nb} positions)`;
  if (input.method === 'identity') {
    const [s, t] = input.fromRange ?? [0, Math.min(na, nb)];
    if (!(int(s) && int(t) && s >= 0 && s < t && t <= Math.min(na, nb))) {
      return err('model.mapping.range', '/fromRange', `[${s}, ${t}) is not within ${sides}.`);
    }
    quads.push(s, t - s, s, 1);
  } else if (input.method === 'segments') {
    const [fc, tc] = [isCyclic(set, from), isCyclic(set, to)];
    /** Per quad, its input segment: one across a cyclic origin is split there. */
    const src: number[] = [];
    for (let k = 0; k < input.segments.length; k++) {
      const [fs, ts, len] = input.segments[k];
      const n = len ?? Math.min(na - fs, nb - ts);
      const fits = (x: number, c: boolean, max: number) => x >= 0 && x < max && n <= (c ? max : max - x);
      if (!(int(fs) && int(ts) && int(n) && n > 0 && fits(fs, fc, na) && fits(ts, tc, nb))) {
        return err('model.mapping.range', `/segments/${k}`, `segment ${k} is not within ${sides}.`);
      }
      for (let [a, b, r] = [fs, ts, n]; r > 0; ) {
        const l = Math.min(r, na - a, nb - b);
        quads.push(a, l, b, 1);
        src.push(k);
        [a, b, r] = [(a + l) % na, (b + l) % nb, r - l];
      }
    }
    // From ranges never overlap; to ranges may (a source used twice).
    const order = Array.from(src.keys()).sort((a, b) => quads[4 * a] - quads[4 * b]);
    for (let i = 1; i < order.length; i++) {
      const [p, k] = [order[i - 1], order[i]];
      if (quads[4 * k] < quads[4 * p] + quads[4 * p + 1]) {
        const later = Math.max(src[p], src[k]);
        return err('model.mapping.range', `/segments/${later}`, `segment ${later} overlaps an earlier one.`);
      }
    }
  } else if (input.method === 'pairwise') {
    if (!/^([0-9]*[MID=X])*$/.test(input.cigar)) {
      return err('model.mapping.cigar', '/cigar', 'a CIGAR is counted operations M, I, D, = and X (`3M1I2M`).');
    }
    let i = 0;
    let j = 0;
    for (const [, count, op] of input.cigar.matchAll(/([0-9]*)([MID=X])/g)) {
      const n = count ? Number(count) : 1;
      if (op !== 'I' && op !== 'D') quads.push(i, Math.min(n, na - i, nb - j), j, 1);
      if (op !== 'D') i += n;
      if (op !== 'I') j += n;
    }
    if (i !== na || j !== nb) {
      return err('model.mapping.cigar', '/cigar', `it covers ${i} and ${j} positions, not ${na} and ${nb}.`);
    }
  } else if (input.method === 'by_numbering') {
    const toName = input.toNumbering ?? input.numbering;
    const nf = numberingOf(set, from, set.strings.find(input.numbering));
    const nt = numberingOf(set, to, set.strings.find(toName));
    if (!nf || !nt) {
      const m = `"${input.from}" has no numbering "${input.numbering}" or "${input.to}" none "${toName}".`;
      return err('model.mapping.numbering', '/numbering', m);
    }
    // Pairs of equal labels, run-length encoded; a to position takes the first from position that labels it.
    const used = new Uint8Array(nb);
    for (let i = 0; i < na; i++) {
      if (!Numbering.labelOf(nf, i as SeqIndex, label)) continue;
      let j = -1;
      if (nt.method === 'names') j = label.name >= 0 ? Numbering.indexOfName(nt, label.name) : -1;
      else if (label.num !== NO_LABEL) j = Numbering.indexOf(nt, label.num, label.suffix);
      if (j < 0 || used[j]) continue;
      used[j] = 1;
      const k = quads.length - 4;
      if (k >= 0 && quads[k] + quads[k + 1] === i && quads[k + 2] + quads[k + 1] === j) quads[k + 1]++;
      else quads.push(i, 1, j, 1);
    }
  } else if (input.method !== 'codon') {
    return err('model.mapping.method', '/method', `unknown method "${(input as { method: string }).method}".`);
  }

  // Statuses from both sides' letters, or elements when both are ligands; a codon's by translation.
  const coded = input.method === 'codon' ? codons(set, input, from, to, quads, err) : undefined;
  if (coded === null) return null;
  const map = Coords.create(from, to, coded ? 3 : 1, quads);
  const { forward, spans } = arraysOf(map, na, true);
  const a = lettersOf(set, from);
  const b = lettersOf(set, to);
  // Letters not held on either side (a lazy molecule, an item coordinate system): no status, none (unknown).
  const lazy = (x: CoordinateSystem) =>
    (x.kind === CoordinateSystemKind.molecule || x.kind === CoordinateSystemKind.unit) && !!moleculeOf(set, x)?.lazy;
  const unknown = lazy(from) || lazy(to);
  const status = coded && spans && !unknown ? coded.status(spans) : new Uint8Array(na).fill(MappingStatus.none);
  for (let i = 0; !coded && !unknown && i < na; i++) {
    const j = forward[i];
    if (j < 0) continue;
    const same = a.element && b.element ? a.element[i] === b.element[j] : a.code[i] === b.code[j];
    status[i] = same ? MappingStatus.exact : MappingStatus.conflict;
  }
  if (!map.count) {
    issues.push(Issue.warning('model.mapping.empty', path, `Mapping "${input.key}" pairs nothing.`));
  }

  const given = input.provenance ? Provenance.create(input.provenance, set.strings) : null;
  const method = set.strings.intern(DEFAULT_METHOD[input.method]);
  const none = -1 as NameId;
  const provenance: Provenance =
    given && given.method >= 0
      ? given
      : given
        ? { ...given, method }
        : { source: none, method, note: none, confidence: Number.NaN, references: new Int32Array(), parameters: null };
  return {
    id: id as MappingId,
    key: set.strings.intern(input.key),
    method: input.method,
    from,
    to,
    map,
    forward,
    status,
    span: coded ? 3 : 1,
    spans,
    codon: coded?.info ?? null,
    provenance,
    custom: CustomData.create(input.custom),
    input,
    deps: [...depsOf(set, from), ...depsOf(set, to)],
  };
}

/** A kept mapping against the set's current endpoints (set.ts `assemble`): 0 unchanged, 1 to rebuild, 2 orphaned. */
function check(set: SeqSet, mapping: Mapping): number {
  const m = mapping as Built;
  const ends = endsOf(set, m.input);
  if (!ends) return 2;
  const deps = [...depsOf(set, ends[0]), ...depsOf(set, ends[1])];
  const same = Coords.sameSystem(ends[0], m.from) && Coords.sameSystem(ends[1], m.to);
  return same && deps.every((d, i) => d === m.deps[i]) ? 0 : 1;
}

/** The row → bound systems tables, by AlignmentId. */
function bind(set: SeqSet): Map<number, RowBindings> {
  const unitsOf = new Map<number, number[]>();
  const { units } = set;
  for (let u = 0; u < units.count; u++) {
    if (!units.alive[u]) continue;
    const list = unitsOf.get(units.molecule[u]);
    if (list) list.push(u);
    else unitsOf.set(units.molecule[u], [u]);
  }
  type Entry = { unit: number; mapping: number; map: CoordinateMap; rowIsFrom: number };
  const perRow = new Map<number, Map<number, { bindings: number; entries: Entry[] }>>();
  for (const m of set.mappings) {
    if (!m) continue;
    for (const side of [0, 1]) {
      const row = side ? m.to : m.from;
      const other = side ? m.from : m.to;
      if (row.kind !== CoordinateSystemKind.row || other.kind === CoordinateSystemKind.row) continue;
      const rows = perRow.get(row.id) ?? perRow.set(row.id, new Map()).get(row.id);
      const at = rows?.get(row.row) ?? rows?.set(row.row, { bindings: 0, entries: [] }).get(row.row);
      if (!at) continue;
      at.bindings++;
      const reached = other.kind === CoordinateSystemKind.unit ? [other.id] : (unitsOf.get(other.id) ?? []);
      for (const unit of reached) at.entries.push({ unit, mapping: m.id, map: m.map, rowIsFrom: side ? 0 : 1 });
    }
  }
  const out = new Map<number, RowBindings>();
  for (const [alignment, rows] of perRow) {
    const order = [...rows.keys()].sort((x, y) => x - y);
    const start = new Int32Array(order.length + 1);
    const entries: Entry[] = [];
    order.forEach((r, k) => {
      entries.push(...(rows.get(r)?.entries ?? []));
      start[k + 1] = entries.length;
    });
    const entryRow = new Int32Array(entries.length);
    for (let k = 0; k < order.length; k++) entryRow.fill(k, start[k], start[k + 1]);
    const byUnit = Int32Array.from(entries.keys()).sort((a, b) => entries[a].unit - entries[b].unit || a - b);
    const units = [...new Set(Array.from(byUnit, (e) => entries[e].unit))];
    const unitStart = new Int32Array(units.length + 1);
    for (let i = 0, k = 0; i < byUnit.length; i++) if (entries[byUnit[i]].unit !== units[k]) unitStart[++k] = i;
    unitStart[units.length] = byUnit.length;
    out.set(alignment, {
      rows: Int32Array.from(order),
      start,
      bindings: Int32Array.from(order, (r) => rows.get(r)?.bindings ?? 0),
      unit: Int32Array.from(entries, (e) => e.unit),
      mapping: Int32Array.from(entries, (e) => e.mapping),
      map: entries.map((e) => e.map),
      rowIsFrom: Uint8Array.from(entries, (e) => e.rowIsFrom),
      entryRow,
      units: Int32Array.from(units),
      unitStart,
      byUnit,
    });
  }
  return out;
}

function bindingsOf(set: SeqSet, alignment: number): RowBindings | null {
  return (set as SeqSet & Assembled).bindings?.get(alignment) ?? null;
}

/** The index of `x` in increasing `xs`, -1: absent. Allocation-free. */
function find(xs: Int32Array, x: number): number {
  let lo = 0;
  let hi = xs.length - 1;
  while (lo <= hi) {
    const mid = (lo + hi) >> 1;
    if (xs[mid] === x) return mid;
    if (xs[mid] < x) lo = mid + 1;
    else hi = mid - 1;
  }
  return -1;
}
/** A row's index in its alignment's table, -1: not bound. */
const boundRow = (b: RowBindings, row: number) => find(b.rows, row);
/** A unit's index in the unit index, -1: no row reaches it. */
const boundUnit = (b: RowBindings, unit: number) => find(b.units, unit);

/** `codon`: the step is a codon mapping's, at most one per path; `unit`: a free step. */
type Edge = { step: PathStep; to: CoordinateSystem; codon: boolean };
const graphs = new WeakMap<SeqSet, Map<string, Edge[]>>();
const step = (mapping: number, reversed: boolean, alignment: number, row: number, unit: number): PathStep => ({
  mapping,
  reversed,
  alignment,
  row,
  unit,
});
const rowEdge = (alignment: number, row: number, reversed: boolean): Edge => ({
  step: step(-1, reversed, alignment, row, -1),
  to: reversed
    ? { kind: CoordinateSystemKind.row, id: alignment, row }
    : { kind: CoordinateSystemKind.columns, id: alignment, row: -1 },
  codon: false,
});

/**
 * Every step from every coordinate system, once per set: mappings both ways, unit steps, and row edges of the rows a
 * mapping names (the columns never fan out to every held row; a search adds its own end rows).
 */
function graph(set: SeqSet): Map<string, Edge[]> {
  const cached = graphs.get(set);
  if (cached) return cached;
  const g = new Map<string, Edge[]>();
  const add = (a: CoordinateSystem, e: Edge) => {
    const list = g.get(keyOf(a));
    if (list) list.push(e);
    else g.set(keyOf(a), [e]);
  };
  const rows = new Set<string>();
  for (const m of set.mappings) {
    if (!m) continue;
    add(m.from, { step: step(m.id, false, -1, -1, -1), to: m.to, codon: m.span > 1 });
    add(m.to, { step: step(m.id, true, -1, -1, -1), to: m.from, codon: m.span > 1 });
    for (const s of [m.from, m.to]) {
      if (s.kind !== CoordinateSystemKind.row || rows.has(keyOf(s))) continue;
      rows.add(keyOf(s));
      add(s, rowEdge(s.id, s.row, false));
      add({ kind: CoordinateSystemKind.columns, id: s.id, row: -1 }, rowEdge(s.id, s.row, true));
    }
  }
  for (let u = 0; u < set.units.count; u++) {
    if (!set.units.alive[u]) continue;
    const unit = { kind: CoordinateSystemKind.unit, id: u, row: -1 };
    const molecule = { kind: CoordinateSystemKind.molecule, id: set.units.molecule[u], row: -1 };
    add(unit, { step: step(-1, false, -1, -1, u), to: molecule, codon: false });
    add(molecule, { step: step(-1, true, -1, -1, u), to: unit, codon: false });
  }
  graphs.set(set, g);
  return g;
}

const unitEnds = new WeakMap<SeqSet, boolean>();
/**
 * Whether a mapping ends on a unit: the only way one copy can reach another (seqspec:model#copies). Without one, a
 * unit's search to a sibling unit is a dead end, so a caller over a merged block skips it (3J3Q: 1,356 copies, 1.8M
 * searches).
 */
function unitsMapped(set: SeqSet): boolean {
  let out = unitEnds.get(set);
  if (out === undefined) {
    out = set.mappings.some(
      (m) => m !== null && (m.from.kind === CoordinateSystemKind.unit || m.to.kind === CoordinateSystemKind.unit),
    );
    unitEnds.set(set, out);
  }
  return out;
}

/**
 * A search state: a system, whether a codon step was taken, the molecule an up step reached (the copy rule: no down
 * step from it) and whether the last step was a down step (no up step right after: it would only clip).
 */
type State = { key: string; system: CoordinateSystem; codon: boolean; up: number; down: boolean };
type Pred = { prev: State; step: PathStep };
const searches = new WeakMap<SeqSet, Map<string, { n: number; steps: PathStep[] }>>();
/** States the last search expanded (a counter). */
let expanded = 0;

/**
 * The shortest path from one coordinate system to another (breadth-first, at most 4 counted steps, unit steps free, at
 * most one codon step, the copy rule in the search) into `out`: its step count, -1 none, -2 ambiguous (several
 * shortest candidates with different composed maps; `out` holds the first). A path that pairs nothing is not a
 * candidate. Cached per set. Edges.
 */
function path(set: SeqSet, from: CoordinateSystem, to: CoordinateSystem, out: PathStep[]): number {
  out.length = 0;
  if (systemLength(set, from) < 0 || systemLength(set, to) < 0) return -1;
  if (keyOf(from) === keyOf(to)) return 0;
  const cache = searches.get(set) ?? searches.set(set, new Map()).get(set);
  const key = `${keyOf(from)}>${keyOf(to)}`;
  const found = cache?.get(key) ?? search(set, from, to);
  cache?.set(key, found);
  out.push(...found.steps);
  return found.n;
}

function search(set: SeqSet, from: CoordinateSystem, to: CoordinateSystem): { n: number; steps: PathStep[] } {
  expanded = 0;
  const g = graph(set);
  const goal = keyOf(to);
  const state = (system: CoordinateSystem, codon: boolean, up: number, down: boolean): State => ({
    key: `${keyOf(system)}|${+codon}|${up}|${+down}`,
    system,
    codon,
    up,
    down,
  });
  /** The end rows' own edges: the start row onto its columns, the goal row from its columns. */
  const edgesOf = (s: CoordinateSystem): readonly Edge[] => {
    const list = g.get(keyOf(s)) ?? [];
    const extra: Edge[] = [];
    if (s === from && s.kind === CoordinateSystemKind.row) extra.push(rowEdge(s.id, s.row, false));
    if (s.kind === CoordinateSystemKind.columns && to.kind === CoordinateSystemKind.row && to.id === s.id)
      extra.push(rowEdge(to.id, to.row, true));
    return extra.length ? [...list, ...extra] : list;
  };
  const start = state(from, false, -1, false);
  const level = new Map<string, number>([[start.key, 0]]);
  const preds = new Map<string, Pred[]>();
  const reach = (s: State, e: Edge, into: State[], at: number) => {
    if (e.codon && s.codon) return;
    const unit = e.step.unit >= 0;
    if (unit && (e.step.reversed ? s.up === s.system.id : s.down)) return;
    const n = state(e.to, s.codon || e.codon, unit && !e.step.reversed ? e.to.id : s.up, unit && e.step.reversed);
    const seen = level.get(n.key);
    if (seen !== undefined && seen !== at) return;
    const pred = { prev: s, step: e.step };
    if (seen === at) preds.get(n.key)?.push(pred);
    else {
      level.set(n.key, at);
      preds.set(n.key, [pred]);
      into.push(n);
    }
  };
  let frontier = [start];
  for (let d = 0; d <= MAX_STEPS && frontier.length; d++) {
    // Free steps within the level, then the goals it reached; else the counted steps make the next level.
    const all = frontier.slice();
    for (const s of frontier) {
      if (keyOf(s.system) !== goal) for (const e of edgesOf(s.system)) if (e.step.unit >= 0) reach(s, e, all, d);
    }
    const goals = all.filter((s) => keyOf(s.system) === goal);
    const found = goals.length ? candidates(set, from, goals, preds) : null;
    if (found) return found;
    // A goal reached only by paths that pair nothing may be reached again, one level on.
    for (const s of goals) level.delete(s.key);
    const next: State[] = [];
    for (const s of all) {
      if (keyOf(s.system) === goal || d === MAX_STEPS) continue;
      expanded++;
      for (const e of edgesOf(s.system)) if (e.step.unit < 0) reach(s, e, next, d + 1);
    }
    frontier = next;
  }
  return { n: -1, steps: [] };
}

/** The candidates' steps (at most 16), those that pair nothing dropped; -2 when two differ. Null when none pairs. */
function candidates(set: SeqSet, from: CoordinateSystem, goals: State[], preds: Map<string, Pred[]>) {
  const paths: PathStep[][] = [];
  const walk = (key: string, tail: PathStep[]) => {
    if (paths.length >= MAX_PATHS) return;
    const list = preds.get(key);
    if (!list) paths.push(tail);
    else for (const p of list) walk(p.prev.key, [p.step, ...tail]);
  };
  for (const s of goals) walk(s.key, []);
  const live = paths.filter((p) => p.length && compose(set, from, p).map.count > 0);
  if (!live.length) return null;
  const first = compose(set, from, live[0]);
  const same = live.every((p) => sameComposed(first, compose(set, from, p)));
  return { n: same ? live[0].length : -2, steps: live[0] };
}

/** Two composed paths that place everything alike: their maps, so `Paths.at`, `Paths.spanAt` and to → from alike. */
function sameComposed(a: ComposedPath, b: ComposedPath): boolean {
  return Coords.areEqual(a.map, b.map);
}

/** A step's map as stored, whether the path walks it reversed, and the system it reaches. */
function stepMap(set: SeqSet, s: PathStep): { map: CoordinateMap; reversed: boolean; to: CoordinateSystem } {
  if (s.mapping >= 0) {
    const m = set.mappings[s.mapping];
    if (!m) throw new Error(`Mappings.compose: no mapping ${s.mapping}.`);
    return { map: m.map, reversed: s.reversed, to: s.reversed ? m.from : m.to };
  }
  if (s.unit >= 0) {
    const map = Coords.ofUnit(set, s.unit as UnitId);
    return { map, reversed: s.reversed, to: s.reversed ? map.from : map.to };
  }
  const a = set.alignments[s.alignment];
  const map = a ? Coords.ofRow(a, s.row as RowIndex) : null;
  if (!map) throw new Error(`Mappings.compose: row ${s.row} of alignment ${s.alignment} is not held.`);
  return { map, reversed: s.reversed, to: s.reversed ? map.from : map.to };
}

const composed = new WeakMap<SeqSet, Map<string, ComposedPath>>();

/**
 * One coordinate map along `steps` from `from`, clipped to every system's domain on the way, with its first to position
 * per held from element and a codon step's spans (sparse: over what it pairs); cached per set and path. A codon step taken from the nucleotides is
 * composed from the other end, so the protein side stays `from`. Edges.
 */
function compose(set: SeqSet, from: CoordinateSystem, steps: readonly PathStep[]): ComposedPath {
  const key = `${keyOf(from)}|${steps.map((s) => `${s.mapping},${+s.reversed},${s.alignment},${s.row},${s.unit}`)}`;
  const cache = composed.get(set) ?? composed.set(set, new Map()).get(set);
  const hit = cache?.get(key);
  if (hit) return hit;
  const parts = steps.map((s) => stepMap(set, s));
  const systems = [from, ...parts.map((p) => p.to)];
  const n = steps.length;
  const to = systems[n];
  const c = parts.findIndex((p) => p.map.span > 1);
  const flip = c >= 0 && parts[c].reversed;
  let m = Coords.identity(flip ? to : from, domainOf(set, flip ? to : from));
  for (let i = 0; i < n; i++) {
    const k = flip ? n - 1 - i : i;
    const next = systems[flip ? k : k + 1];
    const map = parts[k].reversed === flip ? parts[k].map : Coords.invert(parts[k].map);
    m = Coords.compose(m, map, next, domainOf(set, next));
  }
  // A source used twice on the path's from side: the copies' side becomes `from`, so from ranges never overlap.
  if (!flip && m.order === 2) {
    const inverse = Coords.invert(m);
    if (inverse.order !== 2) m = inverse;
  }
  const own = Coords.sameSystem(m.from, from);
  const length = Math.max(0, systemLength(set, from));
  const held = heldOf(m, length, own);
  const span = held.heldSpans ? m.span : 1;
  const result = { from, to, map: m, length, ...held, span, codon: c >= 0 ? steps[c].mapping : -1 };
  cache?.set(key, result);
  return result;
}

export const Mappings = {
  build,
  /** What set.ts needs to keep, rebuild or orphan a mapping. */
  dependents: {
    build,
    check,
    inputOf: (m: Mapping) => (m as Built).input,
    keysOf: (m: Mapping) => [(m as Built).input.from, (m as Built).input.to],
  },
  bind,
  bindingsOf,
  boundRow,
  boundUnit,
  systemOf,
  systemLength,
  domainOf,
  path,
  unitsMapped,
  compose,
  sameComposed,
  /** States the last uncached search expanded. */
  expanded: () => expanded,
};

/** `Paths.fromAt`'s scratch: a hit, and the first element visited. */
const hit = { element: -1, phase: -1 };
const first = { element: -1 };
function firstOf(h: { element: number }) {
  if (first.element < 0 || h.element < first.element) first.element = h.element;
}
/** A visitor of runs, with its context (no closure per call). */
type RunVisit<C> = (start: number, stop: number, ctx: C) => void;
/** The path's map is read from its `from` side (else from its `to`: a codon step taken from the nucleotides). */
const ownSide = (p: ComposedPath) => Coords.sameSystem(p.map.from, p.from);

/**
 * Reading a composed path (seqspec:model#paths): per element of its `from` system, and by ranges either way (a streamed
 * request through a path, ST7). Allocation-free: binary searches over the held runs and the map's segments.
 */
export const Paths = {
  /** Element i's first to position (a residue: its codon's first); -1: none. */
  at(p: ComposedPath, i: number): number {
    const k = i >= 0 ? slotIn(p.runs, p.slots, i) : -1;
    return k < 0 ? -1 : p.heldTo[k];
  },
  /** The j-th of element i's `span` to positions in reading order (span 1: `at`, j 0); -1: none. */
  spanAt(p: ComposedPath, i: number, j: number): number {
    if (p.heldSpans === null) return j === 0 ? Paths.at(p, i) : -1;
    const k = i >= 0 ? slotIn(p.runs, p.slots, i) : -1;
    return k < 0 ? -1 : p.heldSpans[p.span * k + j];
  },
  /** From elements [start, stop) as runs of the to positions they reach, merged where they touch. */
  forEachTo<C>(p: ComposedPath, start: number, stop: number, visit: RunVisit<C>, ctx: C): void {
    if (stop <= start) return;
    if (ownSide(p)) Coords.forEachTo(p.map, start, stop, visit, ctx);
    else Coords.forEachFrom(p.map, start, stop, visit, ctx);
  },
  /** To positions [start, stop) as runs of the from elements reaching them (the inverse image, all copies). */
  forEachFrom<C>(p: ComposedPath, start: number, stop: number, visit: RunVisit<C>, ctx: C): void {
    if (stop <= start) return;
    if (ownSide(p)) Coords.forEachFrom(p.map, start, stop, visit, ctx);
    else Coords.forEachTo(p.map, start, stop, visit, ctx);
  },
  /**
   * The k-th from position reaching to element i in its reading order: through a codon step taken from the nucleotides
   * (the to side the protein's), a residue's k-th base, both parts of a codon an intron splits; else (k 0) the first
   * from element reaching it. -1: none.
   */
  fromAt(p: ComposedPath, i: number, k: number): number {
    if (i < 0) return -1;
    if (!ownSide(p)) return k < p.map.span ? Coords.toAt(p.map, i, k) : -1;
    if (k > 0) return -1;
    first.element = -1;
    Coords.forEachFromAt(p.map, i, hit, firstOf, null);
    return first.element;
  },
  /** From positions a to position covers (3 when a residue is the to side's of a codon step; 1/3 the other way). */
  ratio(p: ComposedPath): number {
    return ownSide(p) ? 1 / p.map.span : p.map.span;
  },
};
