/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The kernel's providers (seqspec:model#providers): annotations computed from other data, as `AnnotationInput`s the set
// builds like any other. They run at the edges (once per input, cached by the compiler) and allocate freely.

import { type Alignment, GAP } from './alignment.js';
import { type AnnotationInput, Pairs } from './annotation.js';
import { Issue, type NameId, type RowIndex, type UnitId } from './common.js';
import { CoordinateSystemKind } from './coordinates.js';
import { Mappings, type PathStep, Paths } from './mapping.js';
import type { Alphabet } from './molecule.js';
import { KeyKind, type KeyTarget, type SeqSet } from './set.js';

export type ColumnMeasure = 'conservation' | 'consensus' | 'identity' | 'gap_fraction';
export interface CoverageOptions {
  key: string;
  /** Unit keys; null: every live unit. */
  of?: string[] | null;
  /** A sequence key: positions are the reference's, through a mapping path to each unit. */
  reference?: string | null;
}
export interface ColumnStatsOptions {
  key: string;
  alignment: string;
  measure: ColumnMeasure;
  /** Rows counted; null: every row. */
  rowSubset?: number[] | null;
  gaps?: 'ignore' | 'count';
  /** null: the alignment's. */
  referenceRow?: number | null;
}
export type ProfileOptions = Omit<ColumnStatsOptions, 'measure' | 'referenceRow'>;
export interface PairedOptions {
  key: string;
  /** Pair annotations by key. */
  pairs: string[];
  /** Unit keys; null: every endpoint (an interface). */
  partner?: string[] | null;
  categories?: string[] | null;
}
export interface ProviderResult {
  annotation: AnnotationInput | null;
  issues: Issue[];
}

/** `coverage`'s default palette (a format `Palette`). */
export const COVERAGE_PALETTE = {
  kind: 'categorical',
  colors: { observed: '#4e79a7', unobserved: '#bab0ac', outside: '#eeeeee', conflict: '#e15759' },
} as const;
const COVERAGE = ['observed', 'unobserved', 'outside', 'conflict'] as const;

/** The k standard letters counted per alphabet, in code order; `T` and `U` are one symbol. */
const STANDARD: Record<Exclude<Alphabet, 'generic'>, string> = {
  protein: 'ACDEFGHIKLMNPQRSTVWY',
  dna: 'ACGT',
  rna: 'ACGU',
};

/** Symbol index per code (-1: left out); the gap is symbol k. */
function symbolsOf(alphabet: Alphabet) {
  if (alphabet === 'generic') return null;
  const letters = STANDARD[alphabet];
  const index = new Int8Array(256).fill(-1);
  for (let i = 0; i < letters.length; i++) index[letters.charCodeAt(i)] = i;
  if (alphabet !== 'protein') index[84] = index[85] = 3;
  index[GAP] = letters.length;
  return { k: letters.length, letters, index };
}

const keyOf = (set: SeqSet, kind: number, id: number) => {
  const t: KeyTarget = { kind: -1, id: -1, row: -1 };
  return (k: string) => (set.byName(k, t) && t.kind === kind ? t.id : id);
};

/** Ranges of equal codes (-1 none) of `codes` as rows: target, start, stop, category. */
function runs(
  codes: Int8Array,
  names: readonly string[],
  target: string,
  out: Record<string, unknown[]>,
  extra?: object,
) {
  for (let i = 0; i < codes.length; ) {
    let j = i + 1;
    while (j < codes.length && codes[j] === codes[i]) j++;
    if (codes[i] >= 0) {
      out.target.push(target);
      out.start.push(i);
      out.stop.push(j);
      out.category?.push(names[codes[i]]);
      for (const [k, v] of Object.entries(extra ?? {})) out[k].push(v);
    }
    i = j;
  }
}

/** Ranges by category: observed, unobserved, outside, conflict (against `reference` through a mapping path). */
function coverage(set: SeqSet, o: CoverageOptions): ProviderResult {
  const issues: Issue[] = [];
  const { units, strings } = set;
  const unitOf = keyOf(set, KeyKind.unit, -1);
  const list = o.of ? o.of.map(unitOf) : Array.from({ length: units.count }, (_, u) => (units.alive[u] ? u : -1));
  (o.of ?? []).forEach((k, j) => {
    if (list[j] < 0) issues.push(Issue.error('model.provider.input', `/of/${j}`, `coverage: no unit "${k}".`));
  });
  const ref = o.reference == null ? -1 : keyOf(set, KeyKind.molecule, -1)(o.reference);
  const refMol = set.molecules[ref];
  if (o.reference != null && !refMol) {
    issues.push(Issue.error('model.provider.input', '/reference', `coverage: no sequence "${o.reference}".`));
    return { annotation: null, issues };
  }
  const cols: Record<string, unknown[]> = { target: [], start: [], stop: [], category: [], id: [] };
  const steps: PathStep[] = [];
  for (const u of list) {
    const mol = u < 0 ? null : set.molecules[units.molecule[u]];
    if (!mol) continue;
    const key = strings.get(units.key[u]);
    const seen = new Uint8Array(mol.length);
    const obs = units.observed[u];
    for (let i = units.rangeStart[u]; i < units.rangeStop[u]; i++) seen[i] = obs ? 1 : 2;
    for (let k = 0; obs && k < obs.length; k += 2) seen.fill(2, obs[k], obs[k + 1]);
    if (!refMol) {
      runs(
        Int8Array.from(seen, (s) => (s === 2 ? 0 : s === 1 ? 1 : 2)),
        COVERAGE,
        key,
        cols,
        { id: key },
      );
      continue;
    }
    const from = { kind: CoordinateSystemKind.molecule, id: ref, row: -1 };
    const found = Mappings.path(set, from, { kind: CoordinateSystemKind.unit, id: u, row: -1 }, steps);
    if (found < 0) {
      // -2: more than one shortest path, and they disagree.
      const [code, how] =
        found === -2 ? ['model.coverage.ambiguous', 'ambiguously'] : ['model.provider.unmapped', 'not'];
      issues.push(Issue.warning(code, '/reference', `coverage: unit "${key}" is ${how} mapped to "${o.reference}".`));
      continue;
    }
    const path = Mappings.compose(set, from, steps);
    const codes = Int8Array.from({ length: path.length }, (_, i) => {
      const j = Paths.at(path, i);
      return j < 0 || !seen[j] ? 2 : seen[j] === 1 ? 1 : mol.code[j] !== refMol.code[i] ? 3 : 0;
    });
    runs(codes, COVERAGE, o.reference as string, cols, { id: key });
  }
  return { annotation: { key: o.key, geometry: 'ranges', target: o.reference ?? undefined, columns: cols }, issues };
}

/** The alignment and its counted rows, or null with an issue. */
function alignmentRows(set: SeqSet, o: ProfileOptions, what: string, issues: Issue[]) {
  const a: Alignment | null = set.alignments[keyOf(set, KeyKind.alignment, -1)(o.alignment)] ?? null;
  const fail = (code: string, path: string, m: string) => {
    issues.push(Issue.error(code, path, `${what}: ${m}`));
    return null;
  };
  if (!a) return fail('model.provider.input', '/alignment', `no alignment "${o.alignment}".`);
  if (a.lazy)
    return fail('model.provider.lazy', '/alignment', "a progressive alignment's statistics come from its source.");
  const symbols = symbolsOf(a.alphabet);
  if (!symbols) return fail('model.provider.alphabet', '/alignment', 'a generic alphabet has no standard letters.');
  const rows = o.rowSubset ?? Array.from({ length: a.rowCount }, (_, r) => r);
  const bad = rows.findIndex((r) => !(Number.isInteger(r) && r >= 0 && r < a.rowCount));
  if (bad >= 0) return fail('model.provider.input', `/row_subset/${bad}`, `no row ${rows[bad]}.`);
  return { a, rows, symbols, kk: symbols.k + (o.gaps === 'ignore' ? 0 : 1) };
}

/** Per column: counts of the k letters and the gap (symbol k) into `counts`; returns n, the counted symbols. */
function countColumn(
  a: Alignment,
  rows: readonly number[],
  index: Int8Array,
  c: number,
  kk: number,
  counts: Float64Array,
) {
  counts.fill(0);
  const C = a.columnCount;
  for (const r of rows) {
    const s = index[a.rows.codes[a.rows.slotOf(r as RowIndex) * C + c]];
    if (s >= 0) counts[s]++;
  }
  let n = 0;
  for (let j = 0; j < kk; j++) n += counts[j];
  return n;
}

const entropy = (counts: Float64Array, kk: number, n: number) => {
  let h = 0;
  for (let j = 0; j < kk; j++) if (counts[j] > 0) h -= (counts[j] / n) * Math.log2(counts[j] / n);
  return h;
};

/** Values over columns, or for `consensus` one range per column (`label` the letter, `value` its share). */
function columnStats(set: SeqSet, o: ColumnStatsOptions): ProviderResult {
  const issues: Issue[] = [];
  const ctx = alignmentRows(set, o, 'column_stats', issues);
  if (!ctx) return { annotation: null, issues };
  const { a, rows, symbols, kk } = ctx;
  const { k, letters, index } = symbols;
  const C = a.columnCount;
  const ref = o.referenceRow ?? a.referenceRow;
  if (o.measure === 'identity' && a.rows.slotOf(ref as RowIndex) < 0) {
    issues.push(Issue.error('model.provider.input', '/reference_row', `column_stats: no reference row ${ref}.`));
    return { annotation: null, issues };
  }
  const counts = new Float64Array(k + 1);
  const value = new Float64Array(C);
  const label: (string | null)[] = [];
  for (let c = 0; c < C; c++) {
    const n = countColumn(a, rows, index, c, kk, counts);
    if (o.measure === 'conservation') value[c] = n ? 1 - entropy(counts, kk, n) / Math.log2(kk) : Number.NaN;
    else if (o.measure === 'gap_fraction') value[c] = rows.length ? counts[k] / rows.length : Number.NaN;
    else if (o.measure === 'identity') {
      const s = index[a.rows.codes[a.rows.slotOf(ref as RowIndex) * C + c]];
      value[c] = s >= 0 && s < k && n ? counts[s] / n : Number.NaN;
    } else {
      let best = -1;
      for (let j = 0; j < k; j++) if (counts[j] > 0 && (best < 0 || counts[j] > counts[best])) best = j;
      label.push(best < 0 ? null : letters[best]);
      value[c] = best < 0 ? Number.NaN : counts[best] / n;
    }
  }
  const base = { key: o.key, target: o.alignment, label: o.measure };
  if (o.measure !== 'consensus') return { annotation: { ...base, geometry: 'values', columns: { value } }, issues };
  const start = Array.from({ length: C }, (_, c) => c);
  return { annotation: { ...base, geometry: 'ranges', columns: { start, label, value } }, issues };
}

/** A matrix over columns × the alphabet's standard letters: frequencies over the counted symbols. */
function profile(set: SeqSet, o: ProfileOptions): ProviderResult {
  const issues: Issue[] = [];
  const ctx = alignmentRows(set, o, 'profile', issues);
  if (!ctx) return { annotation: null, issues };
  const { a, rows, symbols, kk } = ctx;
  const counts = new Float64Array(symbols.k + 1);
  const value: number[][] = [];
  for (let c = 0; c < a.columnCount; c++) {
    const n = countColumn(a, rows, symbols.index, c, kk, counts);
    value.push(Array.from(counts.subarray(0, symbols.k), (x) => (n ? x / n : Number.NaN)));
  }
  const index = Array.from({ length: a.columnCount }, (_, c) => c);
  const annotation: AnnotationInput = {
    key: o.key,
    geometry: 'matrix',
    target: o.alignment,
    states: [...symbols.letters],
    columns: { index, value },
  };
  return { annotation, issues };
}

/** Ranges: positions with a pair to a partner (a pocket, on each partner's counterpart units), or every endpoint. */
function paired(set: SeqSet, o: PairedOptions): ProviderResult {
  const issues: Issue[] = [];
  const { strings, units } = set;
  const marks = new Map<string, Uint8Array>();
  const mark = (key: string, length: number, s: number, t: number) => {
    const m = marks.get(key) ?? marks.set(key, new Uint8Array(length)).get(key);
    m?.fill(1, Math.max(0, s), Math.min(length, t));
  };
  const partners = (o.partner ?? []).map(keyOf(set, KeyKind.unit, -1));
  partners.forEach((u, j) => {
    if (u < 0)
      issues.push(Issue.error('model.provider.input', `/partner/${j}`, `paired: no unit "${o.partner?.[j]}".`));
  });
  o.pairs.forEach((key, j) => {
    const ann = set.annotations[keyOf(set, KeyKind.annotation, -1)(key)];
    if (!ann?.pairs) {
      issues.push(Issue.error('model.provider.input', `/pairs/${j}`, `paired: "${key}" is not a pairs annotation.`));
      return;
    }
    const pc = ann.pairs;
    const ok = (row: number) => !o.categories || o.categories.includes(strings.get(ann.category?.[row] ?? -1));
    if (o.partner) {
      for (const p of partners) {
        if (p < 0) continue;
        Pairs.partnersOf(
          ann,
          p as UnitId,
          (u, i, row) => {
            const mol = set.molecules[units.molecule[u]];
            if (!ok(row) || !mol) return;
            const [s, t] = i < 0 ? [units.rangeStart[u], units.rangeStop[u]] : [i, i + 1];
            mark(strings.get(units.key[u]), mol.length, s, t);
          },
          null,
        );
      }
      return;
    }
    for (let row = 0; row < ann.count; row++) {
      if (!ok(row)) continue;
      for (const side of [0, 1]) {
        const t = side ? pc.bTarget : ann.target;
        const s = side ? pc.bStart[row] : ann.start[row];
        const e = side ? pc.bStop[row] : ann.stop[row];
        const { key: tk, length } = targetKey(set, t, row);
        if (tk) mark(tk, length, s < 0 ? 0 : s, s < 0 ? length : Math.max(e, s + 1));
      }
    }
  });
  const cols: Record<string, unknown[]> = { target: [], start: [], stop: [] };
  for (const [key, m] of marks)
    runs(
      Int8Array.from(m, (x) => x - 1),
      [],
      key,
      cols,
    );
  return { annotation: { key: o.key, geometry: 'ranges', columns: cols }, issues };
}

/** A target's key and length (a whole-unit endpoint over a molecule-level target: the molecule's length). */
function targetKey(
  set: SeqSet,
  t: { molecule: Int32Array; unit: Int32Array; alignment: Int32Array; row: Int32Array },
  i: number,
) {
  const { strings, units } = set;
  if (t.unit[i] >= 0 || t.molecule[i] >= 0) {
    const mol = set.molecules[t.molecule[i]];
    const key = t.unit[i] >= 0 ? strings.get(units.key[t.unit[i]]) : strings.get(mol?.key ?? (-1 as NameId));
    return { key, length: mol?.length ?? 0 };
  }
  const a = set.alignments[t.alignment[i]];
  if (!a) return { key: '', length: 0 };
  const row = t.row[i];
  const key = row >= 0 ? `${strings.get(a.key)}/${row}` : strings.get(a.key);
  return { key, length: row >= 0 ? a.rows.lengthOf(row as RowIndex) : a.columnCount };
}

export const Providers = { coverage, columnStats, profile, paired, symbolsOf };
