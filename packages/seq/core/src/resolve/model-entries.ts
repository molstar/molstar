/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  type Atom,
  CategoryColumn,
  categoryColumn,
  type Issue,
  issueOf,
  isColumn,
  type Node,
  type Value,
} from '@molstar/core/repr/tree';
import {
  AlignmentRuntime,
  AlignmentSystems,
  type FormatRuntime,
  inlineSeq,
  type RangesField,
  type RangesList,
  SeqRuntime,
  type Shape,
  seqToValue,
  type Tables,
} from '@molstar/seqx-core/data';
import { channelAddress, runtimeFor } from '@molstar/seqx-core/data/channels';
import { type QueryItem } from '@molstar/seqx-core/query';
import { readGrid, readPairs } from './held.js';
import { type Placed } from './model-molecules.js';
import { ResolveIssues } from './issues.js';
import { mapEntry } from './maps-entry.js';
import { type UnitEntry } from './molecules.js';
import type {
  AccessPlan,
  AnnotationEntry,
  BoundSourceEntry,
  MapEntry,
  ResolveConfig,
  SourceBinding,
  SystemEntry,
} from './model.js';

/**
 * An alignment's rows system's address: its node's, then `(rows)` (rows are a system, EC-35: positions are rows in
 * data order); its columns system is the node's own address.
 */
export function rowsSystemOf(alignment: string): string {
  return `${alignment}/(rows)`;
}

const p = (n: Node, key: string): Value | undefined => n.params?.[key];

// the transform's params (the records runtime's: _spec/README.md "Data", Records) ride on the entry
const SourceKeys = [
  'uri',
  'ref',
  'value',
  'index',
  'fields',
  'assembly',
  'records',
  'letters',
  'system',
  'base',
  'ends',
  'types',
] as const;

/** A source's entry from its params, bound to the first runtime that reads its format; none is an issue. */
export function bindSource(
  at: Placed,
  runtimes: readonly FormatRuntime[],
): { binding: SourceBinding; issues: Issue[] } {
  const format = p(at.node, 'format');
  const entry: Record<string, unknown> = { address: at.address, format: typeof format === 'string' ? format : '' };
  for (const k of SourceKeys) if (p(at.node, k) !== undefined) entry[k] = p(at.node, k);
  const bound = entry as unknown as BoundSourceEntry;
  const place = (code: 'seq.resolve.source.format' | 'seq.resolve.source.runtime') => ({
    address: at.address,
    hint: ResolveIssues.codes[code].hint,
  });
  if (typeof format !== 'string' || !format)
    return {
      binding: { address: at.address, entry: bound },
      issues: [
        issueOf(
          ResolveIssues,
          'seq.resolve.source.format',
          'The source names no format',
          place('seq.resolve.source.format'),
        ),
      ],
    };
  // one runtime a format (_spec/README.md "Context", Registry entries): several bind none, whatever their order
  const runtime = runtimeFor(runtimes, format);
  if (runtime && runtime !== 'conflict') return { binding: { address: at.address, entry: bound, runtime }, issues: [] };
  return {
    binding: { address: at.address, entry: bound },
    issues: [
      issueOf(
        ResolveIssues,
        'seq.resolve.source.runtime',
        runtime ? 'Several runtimes read the format of the source' : 'No runtime reads the format of the source',
        place('seq.resolve.source.runtime'),
      ),
    ],
  };
}

/** Inline tables as a `seq` source of one tile over the whole system, named by the system's address. */
export function inlinePlan(
  address: string,
  system: { readonly address: string; readonly length?: number; readonly cyclic: boolean },
  shape: Shape,
  tables: Tables,
): AccessPlan {
  const name = system.address;
  const value = seqToValue(inlineSeq({ name, length: system.length, cyclic: system.cyclic }, tables));
  return { kind: 'inline', source: { address, format: 'seq', value }, runtime: SeqRuntime, system: name, shape };
}

/** A source read by region, the system by the name the source knows it by. */
export function sourcePlan(source: string, system: string, shape: Shape): AccessPlan {
  return { kind: 'source', source, system, shape };
}

export function sequenceEntry(at: Placed, source: string | undefined, known: number | undefined): SystemEntry {
  const letters = p(at.node, 'letters');
  const name = p(at.node, 'name');
  const base = {
    address: at.address,
    kind: 'sequence' as const,
    alphabet: p(at.node, 'alphabet') as SystemEntry['alphabet'],
    cyclic: p(at.node, 'cyclic') === true,
    name: typeof name === 'string' ? name : (at.node.id ?? ''),
    numberings: ['index', 'ordinal', ...declaredNumberings(at.node)],
  };
  if (typeof letters === 'string') {
    const tables: Tables = { letters: { form: 'letters', start: 0, letters } };
    const length = letters.length;
    return { ...base, length, letters, access: inlinePlan(at.address, { ...base, length }, 'letters', tables) };
  }
  return {
    ...base,
    // the length its source's coverage gave (`ResolveOptions.lengths`), else unknown
    ...(source && known !== undefined ? { length: known } : {}),
    ...(source ? { source } : {}),
    access: sourcePlan(source ?? '', base.name, 'letters'),
  };
}

/** The ids of a sequence's `numbering` children, in order (a built-in's id is no declared numbering). */
export function declaredNumberings(n: Node): string[] {
  const ids: string[] = [];
  for (const c of n.children ?? [])
    if (c.kind === 'numbering' && c.id !== undefined && c.id !== 'index' && c.id !== 'ordinal') ids.push(c.id);
  return ids;
}

/** The ids of an exact address (`/data/h` → `data`, `h`). */
export function addressIds(address: string): string[] {
  return address.split('/').filter(Boolean);
}

/** An alignment's two systems and its cells, inline (read through `alignment`) or from a source by region. */
export function alignmentEntry(
  at: Placed,
  source: string | undefined,
  knownColumns?: number,
  knownRows?: number,
): { columns: SystemEntry; rows: SystemEntry; annotation: AnnotationEntry } {
  const rows = p(at.node, 'rows') as { label?: Value; letters?: Value } | undefined;
  const alphabet = p(at.node, 'alphabet') as SystemEntry['alphabet'];
  const reference = p(at.node, 'reference');
  const letters = Array.isArray(rows?.letters) ? (rows.letters as string[]) : undefined;
  const label = Array.isArray(rows?.label) ? (rows.label as string[]) : undefined;
  // the table's other columns are the rows' own fields (_spec/README.md "Render", The rows layer, Row columns)
  const fields = rowsFields(rows as Record<string, Value> | undefined);
  const access: AccessPlan =
    letters && label
      ? {
          kind: 'inline',
          source: {
            address: at.address,
            format: 'alignment',
            value: {
              label,
              letters,
              ...(fields ? { fields } : {}),
              ...(firstOf(rows as Record<string, Value>, letters.length) ?? {}),
            } as unknown as Value,
          },
          runtime: AlignmentRuntime,
          system: AlignmentSystems.rows,
          shape: 'cells',
        }
      : sourcePlan(source ?? '', AlignmentSystems.rows, 'cells');
  const system = (address: string, name: string, length: number | undefined, aliases?: string[]): SystemEntry => ({
    address,
    kind: 'sequence',
    alphabet,
    ...(length !== undefined ? { length } : {}),
    cyclic: false,
    name,
    ...(aliases ? { aliases } : {}),
    numberings: ['index', 'ordinal'],
    ...(source ? { source } : {}),
    access,
  });
  const rowsAddress = rowsSystemOf(at.address);
  return {
    // named by its alignment, `<id>/columns`: a source's records shared by several alignments (one Stockholm file of
    // many families) say whose columns they are; its own source (one alignment) says `columns`
    columns: system(
      at.address,
      `${at.node.id ?? ''}/${AlignmentSystems.columns}`,
      letters ? (letters[0]?.length ?? 0) : source ? knownColumns : undefined,
      [AlignmentSystems.columns],
    ),
    rows: system(rowsAddress, AlignmentSystems.rows, letters ? letters.length : source ? knownRows : undefined),
    annotation: {
      address: at.address,
      shape: 'cells',
      target: at.address,
      fields: ['label', 'first', ...Object.keys(fields ?? {})],
      access,
      rows: rowsAddress,
      ...(typeof reference === 'string' ? { reference } : {}),
    },
  };
}

/** An inline rows table's `first` column, as the rows' first residue numbers, when it is one number a row. */
export function firstOf(rows: Record<string, Value> | undefined, n: number): { first: number[] } | undefined {
  const v = rows?.first;
  const list =
    isColumn(v) && !(v instanceof CategoryColumn)
      ? Array.from(v as ArrayLike<number>)
      : Array.isArray(v)
        ? v
        : undefined;
  return list && list.length === n && list.every((x) => typeof x === 'number')
    ? { first: list as number[] }
    : undefined;
}

/** An inline rows table's columns past `label`, `letters` and `first`, a plain list each; undefined when none. */
export function rowsFields(rows: Record<string, Value> | undefined): Record<string, unknown[]> | undefined {
  let out: Record<string, unknown[]> | undefined;
  for (const key of Object.keys(rows ?? {})) {
    // `first` is the rows' first residue numbers (EC-39), never a field
    if (key === 'label' || key === 'letters' || key === 'first') continue;
    const v = rows![key];
    const list =
      v instanceof CategoryColumn
        ? Array.from(v.codes, (k) => (k < 0 ? null : v.dict[k]))
        : isColumn(v)
          ? Array.from(v as ArrayLike<number>)
          : Array.isArray(v)
            ? v
            : undefined;
    if (list) (out ??= Object.create(null) as Record<string, unknown[]>)[key] = list;
  }
  return out;
}

/** A `sequence` node as a system (the overlay's dynamic sequences are resolved so). */
export function resolveSequenceNode(address: string, node: Node, source: string | undefined): SystemEntry {
  return sequenceEntry({ node, address }, source, undefined);
}

/** A `map` node resolved (the overlay's item maps are built as one): its entry and its issues. */
export function resolveMapNode(
  address: string,
  node: Node,
  from: SystemEntry,
  to: SystemEntry,
  config: ResolveConfig,
  unit?: UnitEntry,
): { map: MapEntry; issues: Issue[] } {
  return mapEntry({ node, address }, from, to, config, unit);
}

/** A `ranges` or `values` node resolved on its target (the overlay's dynamic annotations are resolved so). */
export function resolveAnnotationNode(
  address: string,
  node: Node,
  target: SystemEntry,
  source: string | undefined,
): { annotation: AnnotationEntry; issues: Issue[] } {
  return annotationEntry({ node, address }, target, source);
}

/** A `pairs` or `grid` node on its two systems, held whole. */
export function heldEntry(
  at: Placed,
  x: SystemEntry,
  y: SystemEntry,
): { annotation: AnnotationEntry; issues: Issue[] } {
  const base = { address: at.address, target: x.address, ySystem: y.address };
  if (at.node.kind === 'pairs') {
    const pairs = readPairs(at.node);
    const fields = ['type', ...(pairs.labels ? ['label'] : []), ...Object.keys(pairs.extra ?? {})];
    return { annotation: { ...base, shape: 'pairs', fields, pairs }, issues: [] };
  }
  const { grid, issues } = readGrid(at.address, at.node, x.length, y.length);
  return { annotation: { ...base, shape: 'grid', fields: ['value'], ...(grid ? { grid } : {}) }, issues };
}

const Strands: Readonly<Record<string, number>> = { '+': 1, '-': -1 };

export function annotationEntry(
  at: Placed,
  target: SystemEntry,
  source: string | undefined,
): { annotation: AnnotationEntry; issues: Issue[] } {
  const shape = at.node.kind as 'ranges' | 'values';
  const base = { address: at.address, shape, target: target.address };
  // a values channel (`field`) reads its own address of the source (_spec/README.md "Data", Store, Channels)
  const field = p(at.node, 'field');
  const read = source && shape === 'values' && typeof field === 'string' ? channelAddress(source, field) : source;
  const bySource = (s: Shape) => ({
    annotation: { ...base, fields: [], access: sourcePlan(read ?? '', target.name, s) },
    issues: [],
  });
  if (shape === 'values') {
    const values = p(at.node, 'values');
    if (!isColumn(values)) return bySource('values');
    const tables: Tables = { values: { form: 'positions', start: 0, value: values as Float32Array } };
    return {
      annotation: { ...base, fields: [], access: inlinePlan(at.address, target, 'values', tables) },
      issues: [],
    };
  }
  const rows = p(at.node, 'rows') as Record<string, Value> | undefined;
  if (!rows) return bySource('ranges');
  const { items, table, fields, wrappedUnknown } = inlineRanges(rows, target);
  const issues = wrappedUnknown
    ? [
        issueOf(ResolveIssues, 'seq.resolve.length.unknown', 'A row wraps the origin of a sequence of unknown length', {
          address: at.address,
          hint: ResolveIssues.codes['seq.resolve.length.unknown'].hint,
        }),
      ]
    : [];
  return {
    annotation: { ...base, fields, items, access: inlinePlan(at.address, target, 'ranges', { ranges: table }) },
    issues,
  };
}

/**
 * Inline rows as items and a ranges table: a row that wraps a cyclic origin is two rows of one item (EC-9, EC-10); on
 * a sequence of unknown length it is left out, and said.
 */
/** The columns of inline ranges the format reads itself; the others are a table's custom columns. */
const KnownColumns = new Set(['start', 'end', 'id', 'label', 'category', 'strand']);

/**
 * A column of strings as the normalized spec holds it: a category column, or an array where it has nulls (a null reads
 * as none: no strand, no category); absent, or null whole, none.
 */
export function stringsOf(c: Value | undefined): { get(i: number): string | null } | undefined {
  if (c instanceof CategoryColumn) return c;
  if (!Array.isArray(c)) return undefined;
  return { get: (i: number) => (typeof c[i] === 'string' ? (c[i] as string) : null) };
}

export function inlineRanges(rows: Record<string, Value>, target: SystemEntry) {
  const start = rows.start as Int32Array;
  const end = rows.end as Int32Array;
  const id = stringsOf(rows.id);
  const strand = stringsOf(rows.strand);
  const label = Array.isArray(rows.label) ? (rows.label as (string | null)[]) : undefined;
  const category = stringsOf(rows.category);
  const fieldNames = [label && 'label', category && 'category', strand && 'strand'].filter((f): f is string => !!f);
  const items: QueryItem[] = [];
  const s: number[] = [];
  const e: number[] = [];
  const rowOf: number[] = [];
  const length = target.length;
  let wrappedUnknown = false;
  for (let i = 0; i < start.length; i++) {
    const fields: Record<string, Atom> = {};
    if (label) fields.label = label[i] ?? null;
    if (category) fields.category = category.get(i);
    if (strand) fields.strand = strand.get(i);
    const key = id?.get(i) ?? i;
    let pieces: [number, number][] = [[start[i], end[i]]];
    if (start[i] > end[i]) {
      if (length === undefined) {
        wrappedUnknown = true;
        continue;
      }
      pieces = [
        [start[i], length],
        [0, end[i]],
      ];
    }
    for (const [a, b] of pieces) {
      items.push({ key, start: a, end: b, ...(fieldNames.length ? { fields } : {}) });
      s.push(a);
      e.push(b);
      rowOf.push(i);
    }
  }
  const pick = <T>(f: (i: number) => T) => rowOf.map(f);
  const columns: RangesField[] = [];
  if (label) columns.push({ name: 'label', column: categoryColumn(pick((i) => label[i] ?? null)) });
  if (category) columns.push({ name: 'category', column: categoryColumn(pick((i) => category.get(i))) });
  // a table's custom columns ride along (a theme reads them by name: a score, a consequence), in their types
  for (const [name, c] of Object.entries(rows)) {
    if (KnownColumns.has(name)) continue;
    if (c instanceof CategoryColumn) columns.push({ name, column: categoryColumn(pick((i) => c.get(i))) });
    else if (Array.isArray(c))
      columns.push({
        name,
        column: categoryColumn(pick((i) => (c[i] === null || c[i] === undefined ? null : String(c[i])))),
      });
    else if (ArrayBuffer.isView(c)) {
      const typed = c as unknown as Float64Array;
      columns.push({ name, column: new (typed.constructor as Float64ArrayConstructor)(pick((i) => typed[i])) });
    }
  }
  const table: RangesList = {
    form: 'list',
    start: Int32Array.from(s),
    end: Int32Array.from(e),
    ...(strand ? { strand: Int8Array.from(pick((i) => Strands[strand.get(i) ?? ''] ?? 0)) } : {}),
    ...(id ? { id: categoryColumn(pick((i) => id.get(i))) } : { row: Int32Array.from(rowOf) }),
    ...(columns.length ? { fields: columns } : {}),
  };
  return { items, table, fields: fieldNames, wrappedUnknown };
}
