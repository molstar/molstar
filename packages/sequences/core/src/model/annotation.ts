/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Annotations (seqspec:model#annotations): tables of items over targets, with a geometry. Built once from columns at
// load:
// targets resolved and labels converted through each target's numbering once, rows kept in input order, a
// per-slot interval index, dense values per target, fields, pairs with a per-unit index, matrices.
// Queries (`itemsOver`, `valueAt`, `Pairs.touching`, `partnersOf`, `endUnits`) take callbacks and allocate nothing.

import {
  type AnnotationId,
  CustomData,
  columnNumbers,
  Issue,
  type NameId,
  Provenance,
  type ProvenanceInput,
  type SeqIndex,
  type UnitId,
} from './common.js';
import { type CoordinateSystem, CoordinateSystemKind, Coords } from './coordinates.js';
import { LABEL, NO_LABEL, Numbering, type PositionLabel } from './numbering.js';
import { KeyKind, type KeyTarget, type SeqSet, SeqSets } from './set.js';

export type Geometry = 'ranges' | 'values' | 'pairs' | 'matrix' | 'grid' | 'rows';
const GEOMETRIES: readonly string[] = ['ranges', 'values', 'pairs', 'matrix', 'grid', 'rows'];

/** Per row: what it is over. A unit target also names its molecule. */
export interface Targets {
  /** MoleculeId, or -1. */
  readonly molecule: Int32Array;
  /** UnitId, or -1: every unit over the molecule. */
  readonly unit: Int32Array;
  /** AlignmentId, or -1. */
  readonly alignment: Int32Array;
  /** RowIndex: the row's own residues (for `rows`, the row); -1: the columns. */
  readonly row: Int32Array;
}

/** Side a is the main targets and start/stop. */
export interface PairColumns {
  readonly bTarget: Targets;
  /** bStart -1: the whole unit. */
  readonly bStart: Int32Array;
  readonly bStop: Int32Array;
  /** 1: a → b, -1: b → a, 0: none. */
  readonly direction: Int8Array | null;
}

/**
 * count × states, position-major. `max` (an addition to seqspec:model#indices-and-labels): per row, the largest value,
 * for drawing below a pixel.
 */
export interface MatrixColumns {
  readonly states: Int32Array;
  readonly values: Float64Array;
  readonly max: Float64Array;
}

/**
 * A grid's values (seqspec:model#annotations): per row its block of (stop − start) × (bStop − bStart) values, a-major,
 * NaN missing: the row's own typed array where it gave one (a Float32Array kept as given, no copy), else its numbers as
 * doubles. Its b side is `pairs` (`bStart` never -1: a whole b target is its domain).
 */
export interface GridColumns {
  readonly rows: readonly (Float32Array | Float64Array)[];
}

export type Field =
  | { readonly name: NameId; readonly kind: 'number'; readonly data: Float64Array }
  | { readonly name: NameId; readonly kind: 'string'; readonly data: Int32Array }
  | { readonly name: NameId; readonly kind: 'json'; readonly data: readonly unknown[] };

export interface Annotation {
  readonly id: AnnotationId;
  readonly key: NameId;
  readonly geometry: Geometry;
  readonly semanticType: NameId;
  readonly label: NameId;
  /** Rows (0 for a streamed annotation), in input order. */
  readonly count: number;
  /** Rows served by a streamed source by region, never held (seqspec:sources#streamed-sources). */
  readonly streamed: boolean;
  readonly itemCount: number;
  /** ItemIndex per row. */
  readonly item: Int32Array;
  /** NameId per item: its `id`, or -1. */
  readonly itemKey: Int32Array;
  readonly target: Targets;
  /** The distinct targets, in first-row order. */
  readonly slots: Targets;
  /** Per row: its slot. */
  readonly slot: Int32Array;
  /** SeqIndex or ColumnIndex per row; start === stop: a boundary; -1 for `rows`. */
  readonly start: Int32Array;
  readonly stop: Int32Array;
  /** NameId per row, -1: missing. */
  readonly category: Int32Array | null;
  readonly rowLabel: Int32Array | null;
  readonly tooltip: Int32Array | null;
  readonly alt: Int32Array | null;
  /** Values given per position: one array per slot, NaN missing; no rows. */
  readonly dense: readonly Float64Array[] | null;
  /** Values per row, NaN missing; null for dense values. */
  readonly sparse: Float64Array | null;
  readonly fields: readonly Field[];
  readonly links: readonly { label: NameId; template: NameId }[];
  readonly pairs: PairColumns | null;
  readonly matrix: MatrixColumns | null;
  readonly grid: GridColumns | null;
  readonly provenance: Provenance | null;
  readonly custom: CustomData;
}

/** A column: values, a dictionary-coded string column, or arbitrary JSON values (seqspec:format#positions, `Table`). */
export type AnnotationColumn =
  | ArrayLike<unknown>
  | { dict: readonly (string | null)[]; codes: ArrayLike<number | null> }
  | { json: readonly unknown[] };

/**
 * Columns as seqspec:format#annotation-tables names them; positions as indices, or labels with `numbering`. `target`,
 * `numbering` and `bNumbering` are the node's params (additions to the spec): every row's without a column.
 */
export interface AnnotationInput {
  key: string;
  geometry: Geometry;
  semanticType?: string;
  label?: string;
  target?: string;
  numbering?: string;
  bNumbering?: string;
  columns: Record<string, AnnotationColumn>;
  states?: string[];
  /** Streamed: its coordinate systems' refs (one per manifest system), no rows. */
  streamed?: readonly string[];
  links?: { label: string; url: string }[];
  provenance?: ProvenanceInput;
  custom?: object;
  /** Runtime properties (seqspec:model#custom-data), attached to the built annotation's `custom.props`. */
  props?: Record<string, unknown>;
}

interface PairIndex {
  /** Distinct units touched, increasing, and their entries (sorted by start, whole-unit first). */
  readonly units: Int32Array;
  readonly unitStart: Int32Array;
  readonly eRow: Int32Array;
  readonly eSide: Int8Array;
  readonly eStart: Int32Array;
  /** Effective end: at least start + 1. */
  readonly eEnd: Int32Array;
  /** Per entry: the largest `eEnd` of its unit's entries so far, so `touching` searches. */
  readonly eRunMax: Int32Array;
  /** Per (row * 2 + side): the units that endpoint touches. */
  readonly sideStart: Int32Array;
  readonly sideUnits: Int32Array;
}

/** What the set keeps beside the spec's fields: the indices, and what to rebuild from. */
interface Built extends Annotation {
  /** Rows by (slot, start); per slot [slotStart[s], slotStart[s + 1]) of it. */
  readonly order: Int32Array;
  readonly slotStart: Int32Array;
  readonly sStart: Int32Array;
  /** Effective end in `order` (a boundary covers its position), and its running maximum per slot. */
  readonly sEnd: Int32Array;
  readonly runMax: Int32Array;
  readonly pairIndex: PairIndex | null;
  /** Per item, the epoch it was last visited in (`beginItems`). */
  readonly marks: Int32Array;
  epoch: number;
  readonly input: AnnotationInput;
  readonly targetKeys: readonly string[];
  readonly deps: readonly unknown[];
  /** Built on first use: item → rows (`itemRows`) and key → item (`itemOf`). */
  byItem?: { readonly start: Int32Array; readonly rows: Int32Array };
  itemByKey?: Map<number, number>;
}

export const ColumnRole = { target: 0, position: 1, descriptive: 2 } as const;
/** Keyed by column name with no prototype, so a column named `constructor` is a field, not a role (CF-02). */
const roles = (...byRole: string[]): Record<string, number> =>
  Object.assign(
    Object.create(null),
    Object.fromEntries(byRole.flatMap((names, role) => (names ? names.split(' ').map((n) => [n, role]) : []))),
  );
/**
 * Per geometry, the columns the kernel reads and their role; every other column is a field. A target or positional
 * column of another geometry is `model.columns.unknown`. Read by the validator too.
 */
export const ColumnRoles: Readonly<Record<Geometry, Readonly<Record<string, number>>>> = {
  ranges: roles('target', 'start stop beg end ins_code boundary', 'value id category label tooltip alt'),
  values: roles('target', 'start stop beg end index pos ins_code', 'value category label tooltip'),
  pairs: roles(
    'target b_target',
    'start stop beg end ins_code b_start b_stop b_beg b_end b_ins_code',
    'id category label tooltip numbering b_numbering direction',
  ),
  matrix: roles('target', 'index pos ins_code', 'value'),
  grid: roles(
    'target b_target',
    'start stop beg end index pos ins_code b_start b_stop b_beg b_end b_ins_code',
    'value id label tooltip numbering b_numbering',
  ),
  rows: roles('row state node'),
};
/** Columns that place rows, in some geometry. */
const PLACING = new Set(
  Object.values(ColumnRoles).flatMap((r) => Object.keys(r).filter((k) => r[k] !== ColumnRole.descriptive)),
);
/** A key table's member columns: read as string fields too. */
const KEY_COLUMNS = new Set(['state', 'node']);
const scratch = new Int32Array(2);
const MAX_ISSUES = 10;
const isMissing = (v: unknown) => v == null || (typeof v === 'number' && Number.isNaN(v));

function reader(col: AnnotationColumn | undefined): { n: number; at: (i: number) => unknown } | null {
  if (!col) return null;
  if (Array.isArray(col) || ArrayBuffer.isView(col) || !('json' in col || 'codes' in col)) {
    const a = col instanceof Float32Array ? columnNumbers(col, Number.NaN) : (col as ArrayLike<unknown>);
    return { n: a.length, at: (i) => a[i] };
  }
  if ('json' in col) return { n: col.json.length, at: (i) => col.json[i] };
  return {
    n: col.codes.length,
    at: (i) => {
      const c = col.codes[i];
      return c == null || c < 0 ? null : (col.dict[c] ?? null);
    },
  };
}

const emptyTargets = (n: number): Targets => ({
  molecule: new Int32Array(n).fill(-1),
  unit: new Int32Array(n).fill(-1),
  alignment: new Int32Array(n).fill(-1),
  row: new Int32Array(n).fill(-1),
});
const pick = (t: Targets, rows: readonly number[]): Targets => ({
  molecule: Int32Array.from(rows, (r) => t.molecule[r]),
  unit: Int32Array.from(rows, (r) => t.unit[r]),
  alignment: Int32Array.from(rows, (r) => t.alignment[r]),
  row: Int32Array.from(rows, (r) => t.row[r]),
});

/** A resolved target: its space (molecule, unit, row or columns), its length and its numberings by name. */
interface Resolved {
  molecule: number;
  unit: number;
  alignment: number;
  row: number;
  length: number;
  space: number;
}

function resolveKey(set: SeqSet, key: string, rowOf: number): Resolved | null {
  const t: KeyTarget = { kind: -1, id: -1, row: -1 };
  if (!set.byName(key, t)) return null;
  const r: Resolved = { molecule: -1, unit: -1, alignment: -1, row: -1, length: -1, space: -1 };
  if (t.kind === KeyKind.molecule || t.kind === KeyKind.unit) {
    r.unit = t.kind === KeyKind.unit ? t.id : -1;
    r.molecule = r.unit >= 0 ? set.units.molecule[r.unit] : t.id;
    r.length = set.molecules[r.molecule]?.length ?? -1;
    r.space = r.unit >= 0 ? CoordinateSystemKind.unit : CoordinateSystemKind.molecule;
  } else if (t.kind === KeyKind.alignment || t.kind === KeyKind.row) {
    const a = set.alignments[t.id];
    r.alignment = t.id;
    r.row = t.kind === KeyKind.row ? t.row : rowOf;
    // A `rows` record needs no residues: a lazy alignment's rows count though not held.
    r.length = !a ? -1 : r.row >= 0 && rowOf < 0 ? a.rows.lengthOf(r.row as never) : a.columnCount;
    r.space = r.row >= 0 ? CoordinateSystemKind.row : CoordinateSystemKind.columns;
    if (rowOf >= 0 && (t.kind === KeyKind.row || !a || rowOf >= a.rowCount)) return null;
  } else return null;
  return r.length < 0 ? null : r;
}

function numberingFor(set: SeqSet, r: Resolved, name: NameId): Numbering | null {
  if (r.space === CoordinateSystemKind.unit) return SeqSets.numberingOf(set, r.unit as UnitId, name);
  if (r.space === CoordinateSystemKind.molecule) return set.molecules[r.molecule]?.numberings.byName(name) ?? null;
  if (r.space === CoordinateSystemKind.row)
    return set.alignments[r.alignment]?.numberingsOf(r.row as never)?.byName(name) ?? null;
  return null;
}

/** What an annotation's arrays depend on: per target key its resolution and the objects it reads. */
function depsOf(set: SeqSet, keys: readonly string[], unitTable: boolean): unknown[] {
  const deps: unknown[] = unitTable ? [set.units] : [];
  const t: KeyTarget = { kind: -1, id: -1, row: -1 };
  for (const k of keys) {
    if (!set.byName(k, t)) deps.push(null, null, null);
    else if (t.kind === KeyKind.unit)
      deps.push(t.id, set.molecules[set.units.molecule[t.id]], set.units.numberings[t.id]);
    else if (t.kind === KeyKind.molecule) deps.push(t.id, set.molecules[t.id], null);
    else deps.push(t.id, set.alignments[t.id], t.row);
  }
  return deps;
}

/** An annotation from its columns over the set's current targets, or null with errors. */
function build(set: SeqSet, input: AnnotationInput, id: number, issues: Issue[], path: string): Built | null {
  const { strings } = set;
  let reported = 0;
  /** A row left out: the first few reported, then a count. */
  const report = (code: string, col: string, i: number, m: string, warning = false) => {
    const message = `Annotation "${input.key}", row ${i}: ${m}`;
    const issue = (warning ? Issue.warning : Issue.error)(code, `${path}/columns/${col}/${i}`, message);
    if (reported++ < MAX_ISSUES) issues.push(issue);
  };
  const fail = (code: string, at: string, m: string) => {
    issues.push(Issue.error(code, `${path}${at}`, `Annotation "${input.key}": ${m}`));
    return null;
  };
  const g = input.geometry;
  if (!GEOMETRIES.includes(g)) return fail('model.annotation.geometry', '/geometry', `unknown geometry "${g}".`);
  // Every column, and the ones this geometry reads (`ColumnRoles`); the others are fields.
  const cols = new Map<string, NonNullable<ReturnType<typeof reader>>>();
  const read = new Map<string, NonNullable<ReturnType<typeof reader>>>();
  let n = -1;
  for (const [name, col] of Object.entries(input.columns ?? {})) {
    const r = reader(col);
    if (!r) continue;
    if (n >= 0 && r.n !== n) return fail('model.columns.length', `/columns/${name}`, `columns differ in length.`);
    n = r.n;
    if (ColumnRoles[g][name] !== undefined) read.set(name, r);
    else if (PLACING.has(name)) {
      const m = `Annotation "${input.key}": "${name}" is not a ${g} column.`;
      issues.push(Issue.error('model.columns.unknown', `${path}/columns/${name}`, m));
      continue;
    }
    cols.set(name, r);
  }
  const lazy = input.streamed != null;
  if (lazy && !input.streamed?.length)
    return fail('model.annotation.target', '', 'a streamed annotation needs a coordinate system.');
  const count = lazy ? 0 : Math.max(0, n);
  const has = (c: string) => read.has(c);
  const formA = [
    has('start') || has('stop') || has('boundary') || has('index'),
    has('beg') || has('end') || has('pos'),
  ];
  const formB = [has('b_start') || has('b_stop'), has('b_beg') || has('b_end')];
  if ((formA[0] && formA[1]) || (formB[0] && formB[1])) {
    return fail('model.position.form', '/columns', 'one side mixes indices and labels.');
  }
  // An insertion code completes an author number: it needs that side's labels.
  for (const [b, labels] of [
    ['', formA[1]],
    ['b_', formB[1]],
  ] as const) {
    if (has(`${b}ins_code`) && !labels) {
      const m = `${b}ins_code completes an author number in ${b ? 'b_beg' : 'pos or beg'} ("52" and "A" name 52A); this table has no ${b ? 'b_beg' : 'pos or beg'}: give positions as labels in an author numbering, or leave ${b}ins_code out.`;
      return fail('model.position.form', `/columns/${b}ins_code`, m);
    }
  }
  const states = input.states ?? [];
  if (g === 'matrix' && !states.length) return fail('model.annotation.states', '/states', 'a matrix needs states.');
  // Values with no position column are dense: one per position in index order.
  const dense = g === 'values' && !formA[0] && !formA[1];

  // Targets: each distinct key resolved once, each (target, numbering) looked up once.
  const slotOfKey = new Map<string, number>();
  const slotRes: Resolved[] = [];
  const numberings = new Map<string, Numbering | null>();
  const targetKeys = new Set<string>();
  const targetOf = (key: unknown, rowOf: number): number => {
    if (typeof key !== 'string') return -1;
    const k = `${key}#${rowOf}`;
    let s = slotOfKey.get(k);
    if (s !== undefined) return s;
    targetKeys.add(key);
    const r = resolveKey(set, key, rowOf);
    s = r ? slotRes.push(r) - 1 : -1;
    slotOfKey.set(k, s);
    return s;
  };
  /** The one slot of a key table keyed by `state` or `node`: a target of nothing. */
  let keyedSlot = -1;
  const memberSlot = () => {
    if (keyedSlot < 0)
      keyedSlot = slotRes.push({ molecule: -1, unit: -1, alignment: -1, row: -1, length: 0, space: -1 }) - 1;
    return keyedSlot;
  };
  const numberingOf = (slot: number, name: unknown): Numbering | null => {
    if (typeof name !== 'string') return null;
    const k = `${slot}:${name}`;
    if (!numberings.has(k)) numberings.set(k, numberingFor(set, slotRes[slot], strings.intern(name)));
    return numberings.get(k) ?? null;
  };

  const lab1: PositionLabel = { num: NO_LABEL, suffix: -1 as NameId, name: -1 as NameId };
  const lab2: PositionLabel = { num: NO_LABEL, suffix: -1 as NameId, name: -1 as NameId };
  const scratchLabel: PositionLabel = { num: NO_LABEL, suffix: -1 as NameId, name: -1 as NameId };
  const toLabel = (v: unknown, num: Numbering, out: PositionLabel) => {
    out.num = NO_LABEL;
    out.suffix = -1 as NameId;
    out.name = -1 as NameId;
    if (num.method === 'names') {
      out.name = strings.find(String(v));
      return out.name >= 0;
    }
    if (typeof v === 'number') {
      out.num = v;
      return Number.isInteger(v);
    }
    const m = typeof v === 'string' ? LABEL.exec(v) : null;
    if (!m) return false;
    out.num = Number(m[1]);
    out.suffix = (m[2] ? strings.intern(m[2]) : -1) as NameId;
    return true;
  };
  const range = new Float64Array(2);
  const num = (x: unknown) => (typeof x === 'number' ? x : Number.NaN);
  const v = (c: string, i: number) => read.get(c)?.at(i);
  /**
   * Row i's label on a side joined with its insertion code (`ins_code`: none when missing, "", "?" or "."), as
   * [beg, end]; null (reported) leaves the row out: a code that is not letters, one on a range, or on a label that has one.
   */
  const coded = (i: number, b: string, l0: unknown, e0: unknown): [unknown, unknown] | null => {
    const c = v(`${b}ins_code`, i);
    if (isMissing(c) || c === '' || c === '?' || c === '.') return [l0, e0];
    const col = `${b}ins_code`;
    if (typeof c !== 'string' || !/^[A-Za-z]+$/.test(c)) {
      report(
        'model.position.syntax',
        col,
        i,
        `insertion code ${JSON.stringify(c)} is not letters ("A", "AA"); none is "", "?" or ".".`,
      );
      return null;
    }
    const pos = b ? 'b_beg' : has('pos') ? 'pos' : 'beg';
    if (!isMissing(e0) && String(e0) !== String(l0)) {
      const m = `${col} "${c}" names one residue, but ${b}beg ${JSON.stringify(l0)} to ${b}end ${JSON.stringify(e0)} is a range: write a range's ends as whole labels (${b}beg "${l0}${c}"), or give ${pos} alone.`;
      report('model.position.form', col, i, m);
      return null;
    }
    if (typeof l0 === 'string' && !/^-?[0-9]+$/.test(l0)) {
      const m = `${pos} ${JSON.stringify(l0)} is not an author number to complete with ${col} "${c}": give the number alone (${pos} "52", ${col} "A"), or the whole label without ${col}.`;
      report('model.position.form', col, i, m);
      return null;
    }
    const label = `${l0}${c}`;
    return [label, label];
  };
  /** What a numbering's labels are, for a label it lacks: their first and last, and the labels on that number. */
  const around = (n: Numbering, l: PositionLabel): string => {
    let [first, last] = ['', ''];
    const same: string[] = [];
    for (let k = 0; k < n.length; k++) {
      if (!Numbering.labelOf(n, k as SeqIndex, scratchLabel)) continue;
      const f = Numbering.format(n, k as SeqIndex, strings);
      if (!first) first = f;
      last = f;
      if (l.num !== NO_LABEL && scratchLabel.num === l.num && same.length < 8) same.push(f);
    }
    const at = same.length ? `; on ${l.num}: ${same.join(', ')}` : '';
    return first ? ` (its labels run ${first}-${last}${at})` : '';
  };
  /** Side a (b: "b_") of row i into `range`; false (with an issue) leaves the row out. */
  const position = (i: number, slot: number, b: string, numbering: unknown): boolean => {
    const r = slotRes[slot];
    // On a cyclic molecule a range may wrap its origin (start > stop, seqspec:model#cyclic-systems): side a of ranges,
    // values, matrices.
    const cyclic = !b && g !== 'pairs' && g !== 'grid' && r.molecule >= 0 && set.molecules[r.molecule]?.cyclic === true;
    const s0 = b ? v('b_start', i) : (v('start', i) ?? v('index', i));
    let l0 = b ? v('b_beg', i) : (v('beg', i) ?? v('pos', i));
    let e0 = v(b ? 'b_end' : 'end', i);
    if (!isMissing(l0)) {
      const c = coded(i, b, l0, e0);
      if (!c) return false;
      [l0, e0] = c;
    }
    const boundary = !b && !isMissing(v('boundary', i));
    if (boundary) {
      range[0] = range[1] = num(v('boundary', i));
    } else if (isMissing(s0) && isMissing(l0)) {
      if (b) {
        range[0] = range[1] = -1;
        return true;
      }
      // A grid's side a without a position is its target's whole domain.
      if (g === 'grid') {
        range[0] = 0;
        range[1] = r.length;
        return true;
      }
      report('model.position.missing', has('beg') ? 'beg' : 'start', i, 'no position.');
      return false;
    } else if (!isMissing(s0)) {
      range[0] = num(s0);
      const t = v(b ? 'b_stop' : 'stop', i);
      range[1] = isMissing(t) ? range[0] + 1 : num(t);
    } else {
      const num = r.space === CoordinateSystemKind.columns ? null : numberingOf(slot, numbering);
      if (!num) {
        report('model.position.form', `${b}beg`, i, 'labels need a numbering of the target.');
        return false;
      }
      const end = isMissing(e0) ? l0 : e0;
      const okBeg = toLabel(l0, num, lab1);
      const ok = okBeg && toLabel(end, num, lab2);
      if (!ok || !Numbering.rangeOf(num, lab1, lab2, scratch, cyclic)) {
        const named = typeof numbering === 'string' ? ` in numbering "${numbering}"` : '';
        // Written backwards (`beg` after `end`, either placed by order): an error that says what to write (CF-01).
        // On a cyclic sequence only a whole lap is refused (LA-20).
        if (ok && Numbering.backwards(num, lab1, lab2)) {
          const m = cyclic
            ? `beg "${l0}" comes right after end "${end}"${named}: a range across the origin covers less than a whole lap (for all of the sequence, write its first and last labels).`
            : `beg "${l0}" comes after end "${end}"${named}: write beg ≤ end (a range wraps only on a cyclic sequence).`;
          report('model.position.range', `${b}beg`, i, m);
          return false;
        }
        // Else the label that is not there, by name: beg's when it places nowhere, else end's.
        if (ok) Numbering.placeRange(num, lab1, lab2, scratch);
        const missingBeg = !okBeg || (ok && scratch[0] < 0);
        const more = reported < MAX_ISSUES && ok ? around(num, missingBeg ? lab1 : lab2) : '';
        const which = missingBeg ? `${b}beg` : `${b}end`;
        report('model.position.unknown', which, i, `no position "${missingBeg ? l0 : end}"${named}${more}.`, true);
        return false;
      }
      range[0] = scratch[0];
      range[1] = scratch[1];
    }
    const [s, t] = range;
    const wraps = cyclic && !boundary && s > t && t >= 0 && s < r.length;
    if (
      !(Number.isInteger(s) && Number.isInteger(t) && s >= 0 && (boundary ? s === t : s < t || wraps) && t <= r.length)
    ) {
      const why =
        s > t && !boundary
          ? g === 'grid' || g === 'pairs'
            ? `; a ${g} row's sides never wrap a cyclic origin: write the run across the origin as two rows, [${s}, ${r.length}) and [0, ${t})`
            : "; start > stop wraps only a cyclic sequence's origin"
          : '';
      report('model.position.range', `${b}start`, i, `[${s}, ${t}) is not within 0..${r.length}${why}.`);
      return false;
    }
    return true;
  };

  /**
   * A grid row's values, checked against its block's size (`a` × `b` positions): false (reported) leaves the row out. A
   * Float32Array or Float64Array is kept as given; a list of numbers becomes doubles (null: NaN).
   */
  const gridRows: (Float32Array | Float64Array)[] = [];
  const gridRow = (i: number, a: number, b: number): boolean => {
    const list = v('value', i);
    const typed = list instanceof Float32Array || list instanceof Float64Array;
    const n = typed || Array.isArray(list) ? (list as ArrayLike<unknown>).length : -1;
    if (n !== a * b) {
      const got = n < 0 ? 'no list' : `${n}`;
      const m = `${a} × ${b} = ${a * b} values are needed (side a's positions × side b's), a-major; got ${got}.`;
      report('model.annotation.grid', 'value', i, m);
      return false;
    }
    gridRows.push(typed ? list : Float64Array.from(list as unknown[], num));
    return true;
  };
  const kept: number[] = [];
  const wrapRuns = new Int32Array(4);
  const slotA: number[] = [];
  const slotB: number[] = [];
  const startA: number[] = [];
  const stopA: number[] = [];
  const startB: number[] = [];
  const stopB: number[] = [];
  const denseFill: number[] = [];
  const denseBySlot: Float64Array[] = [];
  const denseOrder: number[] = [];
  for (let i = 0; i < count; i++) {
    const rowOf = g === 'rows' ? num(v('row', i)) : -1;
    // A key table's member by `state` or `node` has no target: it is a field, and `target.row` is -1.
    const keyed = g === 'rows' && !has('row') && (has('state') || has('node'));
    if (g === 'rows' && !keyed && !(Number.isInteger(rowOf) && rowOf >= 0)) {
      report('model.annotation.row', 'row', i, 'an alignment row index is needed.');
      continue;
    }
    const key = v('target', i) ?? input.target;
    const slot = keyed ? memberSlot() : targetOf(key, rowOf);
    if (slot < 0) {
      report('model.annotation.target', 'target', i, `"${key}" names no target here.`);
      continue;
    }
    if (g === 'rows') {
      range[0] = range[1] = -1;
    } else if (dense) {
      // Dense values are no rows: each written straight into its target's array, index order.
      const at = denseFill[slot] ?? 0;
      if (at >= slotRes[slot].length) {
        report('model.position.range', 'value', i, 'more values than positions.');
        continue;
      }
      denseFill[slot] = at + 1;
      if (!denseBySlot[slot]) {
        denseBySlot[slot] = new Float64Array(slotRes[slot].length).fill(Number.NaN);
        denseOrder.push(slot);
      }
      denseBySlot[slot][at] = num(v('value', i));
      continue;
    } else if (!position(i, slot, '', v('numbering', i) ?? input.numbering)) continue;
    const s = range[0];
    const t = range[1];
    if (g === 'pairs' || g === 'grid') {
      const bs = targetOf(v('b_target', i) ?? key, -1);
      if (bs < 0) {
        report('model.annotation.target', 'b_target', i, 'side b names no target.');
        continue;
      }
      const bn = v('b_numbering', i) ?? input.bNumbering ?? v('numbering', i) ?? input.numbering;
      if (!position(i, bs, 'b_', bn)) continue;
      // A grid's b side without a range is its b target's whole domain.
      if (g === 'grid' && range[0] < 0) {
        range[0] = 0;
        range[1] = slotRes[bs].length;
      }
      if (g === 'grid' && !gridRow(i, t - s, range[1] - range[0])) continue;
      slotB.push(bs);
      startB.push(range[0]);
      stopB.push(range[1]);
    }
    // A range across a cyclic origin is stored as its two runs, rows of one item.
    const runs = s > t ? Coords.cyclicRuns(s, t, slotRes[slot].length, wrapRuns) : 1;
    for (let k = 0; k < runs; k++) {
      kept.push(i);
      slotA.push(slot);
      startA.push(s > t ? wrapRuns[2 * k] : s);
      stopA.push(s > t ? wrapRuns[2 * k + 1] : t);
    }
  }
  if (reported > MAX_ISSUES) {
    const m = `Annotation "${input.key}": ${reported - MAX_ISSUES} more rows left out.`;
    issues.push(Issue.warning('model.annotation.rows', path, m));
  }
  for (const t of input.streamed ?? []) {
    slotA.push(targetOf(t, -1));
    if (slotA[slotA.length - 1] < 0) return fail('model.annotation.target', '/target', `"${t}" names no target.`);
  }

  // Slots: side a's targets of the rows kept, in first-row order.
  const K = kept.length;
  const all: Targets = emptyTargets(slotRes.length);
  slotRes.forEach((r, s) => {
    all.molecule[s] = r.molecule;
    all.unit[s] = r.unit;
    all.alignment[s] = r.alignment;
    all.row[s] = r.row;
  });
  const used = dense && !lazy ? denseOrder : [...new Set(slotA)];
  const slotTargets = pick(all, used);
  // A target's slot through an index: one record per row of a `rows` table is a slot each.
  const slotIndex = new Int32Array(slotRes.length);
  used.forEach((s, k) => {
    slotIndex[s] = k;
  });
  const slot = Int32Array.from(slotA.slice(0, K), (s) => slotIndex[s]);
  const start = Int32Array.from(startA);
  const stop = Int32Array.from(stopA);

  // Items: rows sharing an `id` are one item; others one each.
  const item = new Int32Array(K);
  const itemKeys: number[] = [];
  const byId = new Map<string, number>();
  kept.forEach((i, k) => {
    const raw = v('id', i);
    if (k > 0 && kept[k - 1] === i) item[k] = item[k - 1];
    else if (isMissing(raw)) item[k] = itemKeys.push(-1) - 1;
    else {
      const idKey = String(raw);
      let it = byId.get(idKey);
      if (it === undefined) {
        it = itemKeys.push(strings.intern(idKey)) - 1;
        byId.set(idKey, it);
      }
      item[k] = it;
    }
  });

  const names = (c: string) =>
    has(c) ? Int32Array.from(kept, (i) => (isMissing(v(c, i)) ? -1 : strings.intern(String(v(c, i))))) : null;
  const hasValue = has('value') && g !== 'matrix' && g !== 'grid';
  const sparse = hasValue ? Float64Array.from(kept, (i) => num(v('value', i))) : null;
  const denseArrays = dense && !lazy ? used.map((s) => denseBySlot[s]) : null;

  const fields: Field[] = [];
  for (const [name, col] of Object.entries(input.columns ?? {})) {
    if (!cols.has(name) || (read.has(name) && !KEY_COLUMNS.has(name))) continue;
    const at = (cols.get(name) as NonNullable<ReturnType<typeof reader>>).at;
    const values = kept.map(at);
    const given = values.filter((x) => !isMissing(x));
    const fname = strings.intern(name);
    const isJson = !Array.isArray(col) && !ArrayBuffer.isView(col) && 'json' in col;
    if (!isJson && given.every((x) => typeof x === 'number')) {
      fields.push({ name: fname, kind: 'number', data: Float64Array.from(values, num) });
    } else if (!isJson && given.every((x) => typeof x === 'string')) {
      const data = Int32Array.from(values, (x) => (isMissing(x) ? -1 : strings.intern(x as string)));
      fields.push({ name: fname, kind: 'string', data });
    } else fields.push({ name: fname, kind: 'json', data: values.map((x) => (x === undefined ? null : x)) });
  }

  const links: { label: NameId; template: NameId }[] = [];
  (input.links ?? []).forEach((l, k) => {
    if (typeof l.url === 'string' && l.url.startsWith('https://')) {
      links.push({ label: strings.intern(l.label), template: strings.intern(l.url) });
    } else
      issues.push(Issue.error('model.annotation.link', `${path}/links/${k}`, 'a link must be an https:// template.'));
  });

  let pairs: PairColumns | null = null;
  if (g === 'pairs' || g === 'grid') {
    const dir = has('direction')
      ? Int8Array.from(kept, (i) => (v('direction', i) === 'ab' ? 1 : v('direction', i) === 'ba' ? -1 : 0))
      : null;
    pairs = {
      bTarget: pick(all, slotB),
      bStart: Int32Array.from(startB),
      bStop: Int32Array.from(stopB),
      direction: dir,
    };
  }

  let matrix: MatrixColumns | null = null;
  if (g === 'matrix') {
    const S = states.length;
    const values = new Float64Array(K * S).fill(Number.NaN);
    const max = new Float64Array(K).fill(Number.NaN);
    kept.forEach((i, k) => {
      const row = v('value', i);
      if (!Array.isArray(row) || row.length !== S) {
        if (row != null) report('model.annotation.matrix', 'value', i, `${S} values are needed.`);
        return;
      }
      for (let j = 0; j < S; j++) {
        const x = num(row[j]);
        values[k * S + j] = x;
        if (x > max[k] || Number.isNaN(max[k])) max[k] = x;
      }
    });
    matrix = { states: Int32Array.from(states, (s) => strings.intern(s)), values, max };
  }

  const grid = g === 'grid' ? { rows: gridRows } : null;
  const custom = CustomData.create(input.custom, input.props);
  const unitTable = g === 'pairs' || g === 'grid';
  const ann: Built = {
    id: id as AnnotationId,
    key: strings.intern(input.key),
    geometry: g,
    semanticType: (input.semanticType == null ? -1 : strings.intern(input.semanticType)) as NameId,
    label: strings.intern(input.label ?? input.key),
    count: K,
    streamed: input.streamed != null,
    itemCount: itemKeys.length,
    item,
    itemKey: Int32Array.from(itemKeys),
    target: pick(all, slotA.slice(0, K)),
    slots: slotTargets,
    slot,
    start,
    stop,
    category: names('category'),
    rowLabel: names('label'),
    tooltip: names('tooltip'),
    alt: names('alt'),
    dense: denseArrays,
    sparse,
    fields,
    links,
    pairs,
    matrix,
    grid,
    provenance: input.provenance ? Provenance.create(input.provenance, strings) : null,
    custom,
    ...intervalIndex(slot, start, stop, used.length),
    pairIndex: pairs ? pairIndex(set, slotA, start, stop, all, pairs, slotB) : null,
    marks: new Int32Array(itemKeys.length),
    epoch: 0,
    input,
    targetKeys: [...targetKeys],
    deps: depsOf(set, [...targetKeys], unitTable),
  };
  return ann;
}

/** Rows sorted by (slot, start), stable; effective ends and their running maximum per slot. */
function intervalIndex(slot: Int32Array, start: Int32Array, stop: Int32Array, slots: number) {
  const order = Int32Array.from(slot.keys()).sort((a, b) => slot[a] - slot[b] || start[a] - start[b] || a - b);
  const slotStart = new Int32Array(slots + 1);
  for (let k = 0; k < slot.length; k++) slotStart[slot[k] + 1]++;
  for (let s = 0; s < slots; s++) slotStart[s + 1] += slotStart[s];
  const sStart = Int32Array.from(order, (r) => start[r]);
  const sEnd = Int32Array.from(order, (r) => Math.max(stop[r], start[r] + 1));
  const runMax = new Int32Array(order.length);
  for (let i = 0; i < order.length; i++) {
    const first = i === 0 || slot[order[i - 1]] !== slot[order[i]];
    runMax[i] = first ? sEnd[i] : Math.max(runMax[i - 1], sEnd[i]);
  }
  return { order, slotStart, sStart, sEnd, runMax };
}

/** Per unit, the endpoints touching it; a molecule-level side touches every live unit over the molecule. */
function pairIndex(
  set: SeqSet,
  slot: readonly number[],
  start: Int32Array,
  stop: Int32Array,
  slots: Targets,
  pairs: PairColumns,
  slotB: readonly number[],
): PairIndex {
  const unitsOf = new Map<number, number[]>();
  const reach = (s: number): readonly number[] => {
    if (slots.unit[s] >= 0) return [slots.unit[s]];
    if (slots.molecule[s] < 0) return [];
    let list = unitsOf.get(slots.molecule[s]);
    if (!list) {
      list = [];
      for (let u = 0; u < set.units.count; u++)
        if (set.units.alive[u] && set.units.molecule[u] === slots.molecule[s]) list.push(u);
      unitsOf.set(slots.molecule[s], list);
    }
    return list;
  };
  const sideStart = new Int32Array(slot.length * 2 + 1);
  const sideUnits: number[] = [];
  const entries: { unit: number; row: number; side: number; start: number; end: number }[] = [];
  for (let r = 0; r < slot.length; r++) {
    for (let side = 0; side < 2; side++) {
      const s = side ? pairs.bStart[r] : start[r];
      const end = side ? pairs.bStop[r] : stop[r];
      for (const u of reach(side ? slotB[r] : slot[r])) {
        sideUnits.push(u);
        entries.push({ unit: u, row: r, side, start: s, end: Math.max(end, s + 1) });
      }
      sideStart[r * 2 + side + 1] = sideUnits.length;
    }
  }
  entries.sort((a, b) => a.unit - b.unit || a.start - b.start || a.row - b.row);
  const units = [...new Set(entries.map((e) => e.unit))];
  const unitStart = new Int32Array(units.length + 1);
  let k = 0;
  entries.forEach((e, i) => {
    if (e.unit !== units[k]) unitStart[++k] = i;
  });
  unitStart[units.length] = entries.length;
  const eEnd = Int32Array.from(entries, (e) => e.end);
  const eRunMax = new Int32Array(entries.length);
  for (let i = 0; i < entries.length; i++) {
    eRunMax[i] = i > 0 && entries[i - 1].unit === entries[i].unit ? Math.max(eRunMax[i - 1], eEnd[i]) : eEnd[i];
  }
  return {
    units: Int32Array.from(units),
    unitStart,
    eRow: Int32Array.from(entries, (e) => e.row),
    eSide: Int8Array.from(entries, (e) => e.side),
    eStart: Int32Array.from(entries, (e) => e.start),
    eEnd,
    eRunMax,
    sideStart,
    sideUnits: Int32Array.from(sideUnits),
  };
}

/** A kept annotation against the set's current targets (set.ts `assemble`): 0 unchanged, 1 to rebuild, 2 orphaned. */
function check(set: SeqSet, ann: Annotation): number {
  const a = ann as Built;
  const deps = depsOf(set, a.targetKeys, a.geometry === 'pairs');
  if (deps.length === a.deps.length && deps.every((d, i) => d === a.deps[i])) return 0;
  const t: KeyTarget = { kind: -1, id: -1, row: -1 };
  return a.targetKeys.every((k) => !set.byName(k, t)) ? 2 : 1;
}

/** The slot over a space; -1: no rows there. Allocation-free. */
function slotOf(ann: Annotation, space: CoordinateSystem): number {
  const t = ann.slots;
  for (let s = 0; s < t.unit.length; s++) {
    const ok =
      space.kind === CoordinateSystemKind.unit
        ? t.unit[s] === space.id
        : space.kind === CoordinateSystemKind.molecule
          ? t.unit[s] < 0 && t.molecule[s] === space.id
          : t.alignment[s] === space.id && t.row[s] === (space.kind === CoordinateSystemKind.row ? space.row : -1);
    if (ok) return s;
  }
  return -1;
}

/** In `order` within the slot: [first row that may cover `start`, first row starting at or after `stop`). */
const bounds = { lo: 0, hi: 0 };
function search(x: Built, slot: number, start: number, stop: number): boolean {
  if (!(slot >= 0 && slot < x.slotStart.length - 1)) return false;
  let a = x.slotStart[slot];
  let b = x.slotStart[slot + 1];
  while (a < b) {
    const m = (a + b) >> 1;
    if (x.sStart[m] < stop) a = m + 1;
    else b = m;
  }
  bounds.hi = a;
  a = x.slotStart[slot];
  while (a < b) {
    const m = (a + b) >> 1;
    if (x.runMax[m] > start) b = m;
    else a = m + 1;
  }
  bounds.lo = a;
  return true;
}

/** Rows of the slot over [start, stop), by start (a boundary at k counts at position k). Allocation-free. */
function itemsOver<C>(
  ann: Annotation,
  slot: number,
  start: number,
  stop: number,
  visit: (row: number, ctx: C) => void,
  ctx: C,
): void {
  const x = ann as Built;
  if (!search(x, slot, start, stop)) return;
  for (let i = bounds.lo, hi = bounds.hi; i < hi; i++) if (x.sEnd[i] > start) visit(x.order[i], ctx);
}

/** Dense: one read; sparse: the last covering row's value (row order); NaN: missing. Allocation-free. */
function valueAt(ann: Annotation, slot: number, index: number): number {
  if (ann.dense) {
    const d = ann.dense[slot];
    return d && index >= 0 && index < d.length ? d[index] : Number.NaN;
  }
  const x = ann as Built;
  if (!ann.sparse || !search(x, slot, index, index + 1)) return Number.NaN;
  let last = -1;
  for (let i = bounds.lo; i < bounds.hi; i++) if (x.sEnd[i] > index && x.order[i] > last) last = x.order[i];
  return last < 0 ? Number.NaN : ann.sparse[last];
}

/** Starts a pass in which `visitItem` reports each item once (location-level item queries). */
function beginItems(ann: Annotation): void {
  (ann as Built).epoch++;
}
/** True the first time an item is seen since `beginItems`. */
function firstVisit(ann: Annotation, item: number): boolean {
  const x = ann as Built;
  if (x.marks[item] === x.epoch) return false;
  x.marks[item] = x.epoch;
  return true;
}

const stringAt = (set: SeqSet, a: Int32Array | null, row: number) =>
  a && a[row] >= 0 ? set.strings.get(a[row]) : undefined;
const n = (x: number) => (Number.isNaN(x) ? undefined : x);

/**
 * A row's column or field by name, as given (tooltips, links, a host's panel; the format's filters, split rows and
 * palettes read it per row, so it makes no closure); undefined: missing.
 */
function rowValue(ann: Annotation, set: SeqSet, row: number, name: string): unknown {
  // A column this geometry does not read is a field (`ColumnRoles`).
  if (ColumnRoles[ann.geometry][name] !== undefined && !KEY_COLUMNS.has(name))
    switch (name) {
      case 'category':
        return stringAt(set, ann.category, row);
      case 'label':
        return stringAt(set, ann.rowLabel, row);
      case 'tooltip':
        return stringAt(set, ann.tooltip, row);
      case 'alt':
        return stringAt(set, ann.alt, row);
      case 'id':
        return ann.itemKey[ann.item[row]] >= 0 ? set.strings.get(ann.itemKey[ann.item[row]]) : undefined;
      case 'value':
        return ann.sparse ? n(ann.sparse[row]) : undefined;
      case 'start':
      case 'stop':
        return ann[name][row];
    }
  const id = set.strings.find(name);
  for (const f of ann.fields) {
    if (f.name !== id) continue;
    if (f.kind === 'number') return n(f.data[row]);
    if (f.kind === 'string') return f.data[row] >= 0 ? set.strings.get(f.data[row]) : undefined;
    return f.data[row] ?? undefined;
  }
  return undefined;
}

/** A row's link: `{field}` replaced by the row's value, URL-encoded; null when the row misses a field. Edges. */
function linkOf(ann: Annotation, set: SeqSet, link: number, row: number): string | null {
  let missing = false;
  const url = set.strings.get(ann.links[link]?.template ?? -1).replace(/\{([^}]+)\}/g, (_, f: string) => {
    const x = rowValue(ann, set, row, f);
    if (x === undefined || x === null) missing = true;
    return encodeURIComponent(String(x));
  });
  return missing || !url ? null : url;
}

/** Item → rows, built on first use and cached: item k's rows are rows[start[k]..start[k + 1]). */
function itemRows(ann: Annotation): { readonly start: Int32Array; readonly rows: Int32Array } {
  const x = ann as Built;
  if (x.byItem) return x.byItem;
  const start = new Int32Array(ann.itemCount + 1);
  for (let r = 0; r < ann.count; r++) start[ann.item[r] + 1]++;
  for (let k = 0; k < ann.itemCount; k++) start[k + 1] += start[k];
  const rows = new Int32Array(ann.count);
  const next = start.slice(0, ann.itemCount);
  for (let r = 0; r < ann.count; r++) rows[next[ann.item[r]]++] = r;
  x.byItem = { start, rows };
  return x.byItem;
}

/** The item whose `id` is `key`; -1: none. Its Map is built on first use (remapping item loci). */
function itemOf(ann: Annotation, key: NameId): number {
  const x = ann as Built;
  if (!x.itemByKey) {
    x.itemByKey = new Map();
    for (let k = 0; k < ann.itemCount; k++) if (ann.itemKey[k] >= 0) x.itemByKey.set(ann.itemKey[k], k);
  }
  return x.itemByKey.get(key) ?? -1;
}

export const Annotations = {
  build,
  slotOf,
  valueAt,
  itemsOver,
  itemRows,
  itemOf,
  rowValue,
  linkOf,
  beginItems,
  firstVisit,
  /** What set.ts needs to keep, rebuild or orphan an annotation. */
  dependents: {
    build,
    check,
    inputOf: (a: Annotation) => (a as Built).input,
    keysOf: (a: Annotation) => (a as Built).targetKeys,
  },
};

/** Every pair with an endpoint at position `index` of `unit`; side 0: a, 1: b. Allocation-free. */
function touching<C>(
  pairs: Annotation,
  unit: UnitId,
  index: SeqIndex,
  visit: (row: number, side: number, ctx: C) => void,
  ctx: C,
): void {
  const p = (pairs as Built).pairIndex;
  const k = p ? unitAt(p, unit) : -1;
  if (!p || k < 0) return;
  let e = p.unitStart[k];
  const end = p.unitStart[k + 1];
  // Whole-unit endpoints sort first (start -1) and touch every position.
  for (; e < end && p.eStart[e] < 0; e++) visit(p.eRow[e], p.eSide[e], ctx);
  // [lo, hi): from the first entry whose running maximum end passes `index` to the last starting at or before it.
  let lo = e;
  let hi = end;
  while (lo < hi) {
    const m = (lo + hi) >> 1;
    if (p.eStart[m] <= index) lo = m + 1;
    else hi = m;
  }
  hi = lo;
  lo = e;
  for (let b = hi; lo < b; ) {
    const m = (lo + b) >> 1;
    if (p.eRunMax[m] > index) b = m;
    else lo = m + 1;
  }
  for (e = lo; e < hi; e++) if (p.eEnd[e] > index) visit(p.eRow[e], p.eSide[e], ctx);
}

/**
 * Every (unit, position, row) paired with `partner`: the other endpoint of each pair touching it, at every
 * unit that endpoint reaches; index -1 for a whole-unit endpoint. Callers aggregate. Allocation-free.
 */
function partnersOf<C>(
  pairs: Annotation,
  partner: UnitId,
  visit: (unit: UnitId, index: SeqIndex, row: number, ctx: C) => void,
  ctx: C,
): void {
  const p = (pairs as Built).pairIndex;
  const pc = pairs.pairs;
  const k = p ? unitAt(p, partner) : -1;
  if (!p || !pc || k < 0) return;
  for (let e = p.unitStart[k], end = p.unitStart[k + 1]; e < end; e++) {
    const row = p.eRow[e];
    const o = 1 - p.eSide[e];
    const s = o ? pc.bStart[row] : pairs.start[row];
    const t = Math.max(o ? pc.bStop[row] : pairs.stop[row], s + 1);
    for (let j = p.sideStart[row * 2 + o], jn = p.sideStart[row * 2 + o + 1]; j < jn; j++) {
      const u = p.sideUnits[j] as UnitId;
      if (s < 0) visit(u, -1 as SeqIndex, row, ctx);
      else for (let i = s; i < t; i++) visit(u, i as SeqIndex, row, ctx);
    }
  }
}

function unitAt(p: PairIndex, unit: number): number {
  let lo = 0;
  let hi = p.units.length - 1;
  while (lo <= hi) {
    const mid = (lo + hi) >> 1;
    if (p.units[mid] === unit) return mid;
    if (p.units[mid] < unit) lo = mid + 1;
    else hi = mid - 1;
  }
  return -1;
}

/** The units one endpoint stands for: its unit, or every live unit over its molecule. Allocation-free. */
function endUnits<C>(pairs: Annotation, row: number, side: number, visit: (unit: UnitId, ctx: C) => void, ctx: C) {
  const p = (pairs as Built).pairIndex;
  if (!p || !(row >= 0 && row < pairs.count)) return;
  const at = row * 2 + (side ? 1 : 0);
  for (let j = p.sideStart[at], n = p.sideStart[at + 1]; j < n; j++) visit(p.sideUnits[j] as UnitId, ctx);
}

export const Pairs = { touching, partnersOf, endUnits };
