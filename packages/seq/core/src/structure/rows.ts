/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An annotation's rows as targets on units (_spec/README.md "Structures", Rows): each item (a range's rows, a pair's
// two ends) carried onto the bound units with its fields, for a host to write as its own representation. Inline rows
// at once; a source's from the tiles the store holds at its finest level, the rest pending (never read as empty,
// EC-169). Ungated: a representation is written once, not a frame.

import { CategoryColumn, type Atom, type Issue, issueOf } from '@molstar/core/repr/tree';
import type { DataEntry, DataStore, Want } from '@molstar/seqx-core/data';
import {
  EmptyLocation,
  type Location,
  location,
  type PartInput,
  type QueryModel,
  rangesLocation,
} from '@molstar/seqx-core/query';
import { type AnnotationEntry, overlayQueryModel } from '@molstar/seqx-core/resolve';
import { StructureSyncIssues } from './issues.js';
import { type StructureBinding, type StructureModel, type StructureReads, structureReads } from './model.js';
import { type StructureTarget, toUnits } from './units.js';

/** What the translation reads of the store: its coverage and the tiles it holds. */
export type StoreRead = Pick<DataStore, 'coverage' | 'keysOf' | 'get'>;

export interface HeldOptions {
  /** Asks for what is not held (whoever calls asks: the library holds no consumer). */
  readonly want?: (want: Want) => void;
  /** A molecule unit whole or on its atoms (as `toUnits`). */
  readonly ligand?: 'whole' | 'atoms';
}

export interface UnitRow {
  /** The item's key (its `id`, else its row). */
  readonly key: string | number;
  readonly target: StructureTarget;
  readonly fields: Readonly<Record<string, Atom>>;
}

export interface UnitRowsResult {
  /** Nothing pending. */
  readonly complete: boolean;
  /** What the store does not hold yet, on the annotation's system. */
  readonly pending: Location;
  /** The items reaching a bound unit, in their order. */
  readonly rows: readonly UnitRow[];
  readonly issues: readonly Issue[];
}

/** @internal An item read: its parts on their systems and its fields. */
export interface ReadItem {
  readonly key: string | number;
  readonly parts: PartInput[];
  readonly fields: Record<string, Atom>;
  readonly start: number;
  readonly end: number;
}

/** An annotation's items carried onto the units of the bound complexes (every complex when absent). */
export function unitRows(
  model: StructureModel,
  store: StoreRead,
  annotation: string,
  complexes?: readonly StructureBinding[],
  options: HeldOptions = {},
): UnitRowsResult {
  const r = structureReads(model);
  const read = annotationItems(r, model, store, annotation, options.want);
  if (!read) return noRows(annotation);
  const rows: UnitRow[] = [];
  const issues: Issue[] = [];
  for (const it of read.items) {
    const t = toUnits(model, location({ parts: it.parts }), complexes, {
      limits: null,
      ...(options.ligand ? { ligand: options.ligand } : {}),
    });
    for (const i of t.issues) if (!issues.some((j) => j.code === i.code)) issues.push(i);
    if (t.target.parts.length) rows.push({ key: it.key, target: t.target, fields: it.fields });
  }
  return { complete: !read.pending.parts.length, pending: read.pending, rows, issues };
}

function noRows(annotation: string): UnitRowsResult {
  return {
    complete: true,
    pending: EmptyLocation,
    rows: [],
    issues: [
      issueOf(StructureSyncIssues, 'seq.structure.rows', `${annotation} is no ranges or pairs annotation`, {
        address: annotation,
      }),
    ],
  };
}

/** @internal A ranges or pairs annotation's items, inline or from the store's tiles; undefined when it is neither. */
export function annotationItems(
  r: StructureReads,
  model: StructureModel,
  store: StoreRead,
  address: string,
  want?: (w: Want) => void,
): { items: ReadItem[]; pending: Location; numeric: Set<string> } | undefined {
  const entry: AnnotationEntry | undefined =
    r.model.annotations.get(address) ?? ('overlay' in model ? model.overlay.annotation(address) : undefined);
  if (!entry || (entry.shape !== 'ranges' && entry.shape !== 'pairs')) return undefined;
  const query: QueryModel = 'overlay' in model ? overlayQueryModel(model.model, model.overlay) : r.model;
  const answer = query.annotation(address)?.items();
  const byKey = new Map<string | number, ReadItem>();
  const numeric = new Map<string, boolean>();
  const add = (key: string | number, system: string, start: number, end: number, fields: Record<string, Atom>) => {
    let it = byKey.get(key);
    if (!it) byKey.set(key, (it = { key, parts: [], fields, start, end }));
    it.parts.push({ system, ranges: [start, end] });
    for (const [k, v] of Object.entries(fields))
      if (v !== null) numeric.set(k, (numeric.get(k) ?? true) && typeof v === 'number');
  };
  for (const it of answer?.items ?? []) add(it.key, it.system ?? entry.target, it.start, it.end, { ...it.fields });
  let pending = EmptyLocation;
  const access = entry.access;
  if (answer?.pending && access?.kind === 'source' && access.shape === 'ranges') {
    const length = r.system(entry.target)?.length;
    const held = readHeld(store, access.source, access.system, length, want);
    const seen = new Set<string>();
    for (const e of held.entries) {
      const t = e.tables.ranges;
      if (t?.form !== 'list') continue;
      for (let i = 0; i < t.start.length; i++) {
        const key = t.id?.get(i) ?? t.row?.[i] ?? `${t.start[i]}:${t.end[i]}`;
        const k = `${key}\u0000${t.start[i]}`;
        if (seen.has(k)) continue;
        seen.add(k);
        const fields: Record<string, Atom> = {};
        for (const f of t.fields ?? []) {
          const c = f.column;
          fields[f.name] = c instanceof CategoryColumn ? c.get(i) : ((c as ArrayLike<number>)[i] ?? null);
          numeric.set(f.name, !(c instanceof CategoryColumn));
        }
        add(key, entry.target, t.start[i], t.end[i], fields);
      }
    }
    pending = held.pending.length ? rangesLocation(entry.target, held.pending) : EmptyLocation;
  } else if (answer?.pending) pending = answer.pending;
  const numericFields = new Set([...numeric].filter(([, n]) => n).map(([k]) => k));
  return { items: [...byKey.values()], pending, numeric: numericFields };
}

/**
 * @internal The tiles a source's system holds at its finest level over `[0, length)` (the whole system), and the
 * positions not held: pending, failed or not asked (EC-169). An unknown source or length is pending whole.
 */
export function readHeld(
  store: StoreRead,
  source: string,
  system: string,
  length: number | undefined,
  want?: (w: Want) => void,
): { entries: Extract<DataEntry, { status: 'loaded' }>[]; pending: number[] } {
  const coverage = store.coverage(source);
  if (!coverage || length === undefined) return { entries: [], pending: length === undefined ? [] : [0, length] };
  const w: Want = { source, system, start: 0, end: length, level: coverage.levels[0] ?? 1 };
  want?.(w);
  const entries: Extract<DataEntry, { status: 'loaded' }>[] = [];
  const pending: number[] = [];
  for (const key of store.keysOf(w)) {
    const e = store.get(key);
    if (e?.status === 'loaded') entries.push(e);
    else if (e?.status === 'empty') continue;
    else if (e) pending.push(Math.max(0, e.region.start), Math.min(length, e.region.end));
    else pending.push(0, length);
  }
  return { entries, pending };
}
