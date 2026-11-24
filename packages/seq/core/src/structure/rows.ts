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
import type { Want } from '@molstar/seqx-core/data';
import {
  type LayerSource,
  readTiles,
  type StoreReader,
  type TileContext,
  type TileRead,
} from '@molstar/seqx-core/render';
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

/** What the translation reads of the store: the layers' reader (its coverage, tiles, bound). */
export type StoreRead = StoreReader;

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
  /** Whether the store refused the whole system's tiles as past its bound (its positions in `pending`). */
  readonly tooMany: boolean;
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
  return {
    complete: !read.pending.parts.length && !read.unknown,
    pending: read.pending,
    tooMany: read.tooMany,
    rows,
    issues,
  };
}

function noRows(annotation: string): UnitRowsResult {
  return {
    complete: true,
    pending: EmptyLocation,
    tooMany: false,
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
): { items: ReadItem[]; pending: Location; numeric: Set<string>; tooMany: boolean; unknown: boolean } | undefined {
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
  let tooMany = false;
  let unknown = false;
  const access = entry.access;
  if (answer?.pending && access?.kind === 'source' && access.shape === 'ranges') {
    const length = r.system(entry.target)?.length;
    const held = readHeld(store, { source: access.source, system: access.system, length }, want);
    tooMany ||= held.tooMany;
    unknown ||= held.unknown;
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
  return { items: [...byKey.values()], pending, numeric: numericFields, tooMany, unknown };
}

/** What a held read gives: the loaded tiles, the positions not held, and whether the store's bound refused a want. */
export interface HeldRead {
  readonly entries: TileRead['entries'];
  /** `[start, end, …]` on the system: pending, failed, past the bound or not asked (EC-169). */
  readonly pending: number[];
  /** Past the store's bound on a want's tiles (`seq.store.too_many_tiles`): pending, never read as empty. */
  readonly tooMany: boolean;
  /** The source or the system's length unknown yet: nothing can be held. */
  readonly unknown: boolean;
}

/**
 * @internal The tiles a source's system holds at its finest level over the whole system, through the layers' own
 * reading (`readTiles`, _spec/README.md "Render", Data): a tile not held, failed or past the store's bound is pending
 * (EC-169), never read as empty. `want` asks for the region (no margin: the whole system is asked).
 */
export function readHeld(
  store: StoreReader,
  src: Pick<LayerSource, 'source' | 'system' | 'length'> & Partial<LayerSource>,
  want?: (w: Want) => void,
): HeldRead {
  const length = src.length;
  if (length === undefined) return { entries: [], pending: [], tooMany: false, unknown: true };
  const level = store.coverage(src.source)?.levels[0] ?? 1;
  const ctx: TileContext = { want: (w) => want?.(w), config: { margin: 0 } };
  const full: LayerSource = { dataSystem: src.system, shape: 'ranges', ratio: 1, ...src };
  const read = readTiles(ctx, store, full, [0, length], level);
  return {
    entries: read.entries,
    pending: [...read.pending, ...read.failed],
    tooMany: !!read.tooMany,
    unknown: store.coverage(src.source) === undefined,
  };
}
