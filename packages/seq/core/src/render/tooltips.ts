/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, isColumn, type Value } from '@molstar/core/repr/tree';
import { type MarkdownInline, markdownInlineText, parseMarkdownInline } from '@molstar/core/util/markdown';
import type { RangesList } from '@molstar/seqx-core/data';
import { type IntervalIndex, intervalIndex, meetingInto } from '@molstar/seqx-core/data/intervals';
import { DefaultSpecBounds, type TooltipSchema } from '@molstar/seqx-core/lang';
import { type Location, rangesLocation } from '@molstar/seqx-core/query';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import { MvsSchemaColumns, tableLength } from '@molstar/seqx-core/resolve/tooltip-tables';
import type { ResolvedLayer, ResolvedTooltip, ResolvedView } from '@molstar/seqx-core/view';
import { layerSource } from './layer-data.js';
import type { StoreReader } from './layer-type.js';

// The spec's tooltips at a place (_spec/README.md "Render", Tooltips: `tooltipsAt`): the rows the tooltips of a view
// give a hovered mark, as data (a label, a plain value, the Markdown tree). By data: the hit's item's own row (its data
// node's or its layer's, looked up by key in the store). By range: every row of a table covering the place (a range by
// overlap of one position or more, a cell or a pair by either place), in position order; a `ranges` table through the
// store's entry and its interval index (built once an entry, kept with it), an `address` or MVS table read whole and
// resolved by the caller's `locate`. `views`' tooltips before the view's, each in spec order. A tooltip's one row on
// the mark is plain; several matches are item rows sharing the card's cap, the rest counted (`more`), a row met by
// several parts (a cell's two positions, a pair's two ends) counted once; a tooltip's rows in position order. A warm
// lookup allocates nothing past its matches: a row's text and its row are made once (cached by the table's row and by
// text), an inline table's addresses and locations once a model (read to `tooltipRows`). Values are escaped
// before the Markdown parse (a value adds no mark or link); the parse is cached by text and cleaned after it (controls
// and bidi overrides removed); only Seq*'s own links (`!seq-…`) stay links. A table pending under the mark gives a
// `loading` row (EC-169); a failed one no rows, its tooltip listed in `failed` (EC-174).

/** What a lookup reads: the model, the view resolved (its tooltips), the store, and how an address is located. */
export interface TooltipContext {
  readonly model: StaticModel;
  readonly view: ResolvedView;
  readonly data: StoreReader;
  /** An address as positions over what is held (a viewer's address model, an export's query model); absent: none. */
  locate?(address: string): Location | undefined;
}

/** A stretch of positions the hovered mark covers on a system. */
export interface TooltipPart {
  readonly system: string;
  readonly start: number;
  readonly end: number;
}

/**
 * Where a hover is, as tooltips read it: the hit's layer and item (by data), and the positions under the mark (by range:
 * the pointer's position on it, a cell's two, a pair's two ends, a bin's stretch).
 */
export interface TooltipPlace {
  readonly layer?: ResolvedLayer;
  /** The hit's item: its key and its first position on its data's system (its row is looked up by them). */
  readonly item?: { readonly key: string | number; readonly start: number };
  readonly parts: readonly TooltipPart[];
}

/** A row a tooltip gives: its label (the text itself when the tooltip has none), its plain value, its Markdown tree. */
export interface TooltipRow {
  readonly label: string;
  readonly value: string;
  readonly rich?: readonly MarkdownInline[];
  /** One of several matches (a bin's, several rows over a mark): what the card's cap counts. */
  readonly item?: boolean;
  /** The table under the mark has not arrived. */
  readonly loading?: boolean;
}

/** What a lookup adds to: the rows, the matches past the cap not made into rows, the tooltips whose table failed. */
export interface TooltipOut {
  readonly rows: TooltipRow[];
  more: number;
  readonly failed: string[];
}

/** A formatted text's bound, characters (`SpecBounds.tooltipText`). */
const MaxText = 1000;

/** The rows a place's tooltips give, into `out` (at most `cap` match rows a tooltip; the rest counted in `more`). */
export function tooltipsAt(ctx: TooltipContext, place: TooltipPlace, out: TooltipOut, cap = 8): void {
  const tips = ctx.view.tooltips;
  if (!tips) return;
  for (const t of tips) {
    if (!t.template || t.category) continue;
    if (t.data !== undefined || t.layer !== undefined) byData(ctx, t, place, out);
    else if (t.addresses) byAddresses(ctx, t, place, out, cap);
    else if (t.schema === 'ranges') byRanges(ctx, t, place, out, cap);
    else byTable(ctx, t, place, out, cap);
  }
}

// ---- by data ----

function byData(ctx: TooltipContext, t: ResolvedTooltip, place: TooltipPlace, out: TooltipOut): void {
  const l = place.layer;
  const item = place.item;
  if (!l || !item || !l.data) return;
  if (
    t.layer !== undefined
      ? l.address !== t.layer
      : l.data.address !== t.data && !l.data.address.startsWith(`${t.data}/`)
  )
    return;
  const src = layerSource(ctx.model, l);
  if (!src) return;
  // the item's row in the loaded list tile holding its first position (a layer reads one of its source's levels)
  for (const level of ctx.data.coverage(src.source)?.levels ?? [1])
    for (const k of ctx.data.keysOf({
      source: src.source,
      system: src.system,
      start: item.start,
      end: item.start + 1,
      level,
    })) {
      const e = ctx.data.get(k);
      const table = e?.status === 'loaded' ? e.tables.ranges : undefined;
      if (!table || table.form !== 'list') continue;
      for (let i = 0; i < table.start.length; i++)
        if (table.start[i] === item.start && rowKey(table, i) === item.key) {
          const text = tableText(t, table, i);
          if (text !== undefined) out.rows.push(rowOf(t, text, false));
          return;
        }
    }
}

// ---- by range: a `ranges` table through the store ----

/** Each loaded list table's interval index, built at its first lookup and kept as long as the store keeps the entry. */
const Indices = new WeakMap<RangesList, IntervalIndex>();
/** The matches of one lookup (reused: a warm lookup allocates no buffer). */
let Matches = new Int32Array(16);
/** The stretches a lookup has read (a source's, on its system): a row met again by a later part is counted once. */
const Seen = { source: [] as string[], start: [] as number[], end: [] as number[] };
/** The project's carry of a part, reused when the part is on the table's system. */
const Same = [0, 0];

function indexOf(table: RangesList): IntervalIndex {
  let ix = Indices.get(table);
  if (!ix) Indices.set(table, (ix = intervalIndex(table.start, table.end)));
  return ix;
}

/** Whether a row was met by a stretch read before (the same source): counted there. */
function seenBefore(source: string, s: number, e: number, upTo: number): boolean {
  for (let k = 0; k < upTo; k++) if (Seen.source[k] === source && meets(s, e, Seen.start[k], Seen.end[k])) return true;
  return false;
}

function byRanges(ctx: TooltipContext, t: ResolvedTooltip, place: TooltipPlace, out: TooltipOut, cap: number): void {
  const a = ctx.model.annotations.get(t.address);
  if (!a) return;
  const made = Made.begin();
  let count = 0;
  let pending = false;
  let seen = 0;
  for (const part of place.parts) {
    // the table's plan on the part's system (an assembly's by its sequence), the part carried onto the table's system
    const plan = a.bySystem ? a.bySystem.get(part.system) : a.access;
    const dataSystem = a.bySystem ? part.system : a.target;
    if (!plan) continue;
    const source = plan.kind === 'inline' ? plan.source.address : plan.source;
    let ranges: ArrayLike<number> = Same;
    if (dataSystem === part.system) {
      Same[0] = part.start;
      Same[1] = part.end;
    } else
      ranges =
        ctx.model
          .project(rangesLocation(part.system, [part.start, part.end]), dataSystem)
          .location.parts.find((p) => p.system === dataSystem)?.ranges ?? [];
    const levels = ctx.data.coverage(source)?.levels;
    if (!levels) {
      pending = true;
      continue;
    }
    const before = seen;
    for (let r = 0; r + 1 < ranges.length; r += 2) {
      const s = ranges[r];
      const e = Math.max(ranges[r + 1], s + 1);
      Seen.source[seen] = source;
      Seen.start[seen] = s;
      Seen.end[seen++] = e;
      // the level read: the finest whose tiles here hold rows (a list) and have come; a coarser level's may be
      // summaries, and an inline list is served whole at a coarse level too (what a layer at the fit wanted)
      const level = listLevel(ctx.data, source, plan.system, s, e, levels);
      if (level === undefined) {
        // nothing here has come: pending, or failed where a tile failed
        for (const l of levels)
          for (const k of ctx.data.keysOf({ source, system: plan.system, start: s, end: e, level: l })) {
            const st = ctx.data.get(k)?.status;
            if (st === 'failed') {
              if (!out.failed.includes(t.address)) out.failed.push(t.address);
            } else if (st !== 'empty' && st !== 'loaded') pending = true;
          }
        continue;
      }
      for (const k of ctx.data.keysOf({ source, system: plan.system, start: s, end: e, level })) {
        const entry = ctx.data.get(k);
        if (!entry || entry.status !== 'loaded') {
          if (entry?.status === 'failed' && !out.failed.includes(t.address)) out.failed.push(t.address);
          else if (entry?.status !== 'empty') pending = true;
          continue;
        }
        const table = entry.tables.ranges;
        if (!table || table.form !== 'list') continue;
        if (Matches.length < cap + 1) Matches = new Int32Array(cap + 1);
        const n = meetingInto(indexOf(table), table.start, table.end, s, e, Matches);
        for (let j = 0; j < n; j++) {
          // a row over a tile's edge is in each tile it reaches: counted in the one holding where it meets the place
          const i = j < Matches.length ? Matches[j] : -1;
          if (i < 0) {
            count++;
            continue;
          }
          const at = Math.max(table.start[i], s);
          if (at < entry.region.start || at >= entry.region.end) continue;
          // a row an earlier part met (a cell's two positions, a pair's two ends on one note): counted there
          if (seenBefore(source, table.start[i], table.end[i], before)) continue;
          count++;
          if (made.n >= cap) continue;
          const text = tableText(t, table, i);
          if (text !== undefined) made.add(text, table.start[i]);
        }
      }
    }
  }
  if (pending && !made.n) {
    out.rows.push(loadingRow(t));
    return;
  }
  emit(t, made, count, out);
}

/** The finest of `levels` with a loaded list tile over `[s, e)`; undefined when none has one. */
function listLevel(
  data: StoreReader,
  source: string,
  system: string,
  s: number,
  e: number,
  levels: readonly number[],
): number | undefined {
  for (const level of levels)
    for (const k of data.keysOf({ source, system, start: s, end: e, level })) {
      const entry = data.get(k);
      if (entry?.status === 'loaded' && entry.tables.ranges?.form === 'list') return level;
    }
  return undefined;
}

// ---- by range: addresses, an `address` table, MVS's schemas (read whole) ----

/** A tooltip's rows located, per model: each row's address, its Location (undefined: not resolved yet), its start. */
interface Located {
  readonly addresses: readonly string[];
  readonly locations: (Location | undefined)[];
  /** Each row's first position as located (Infinity until located): the rows' order. */
  readonly starts: number[];
}
const LocatedRows = new WeakMap<StaticModel, Map<string, Located>>();

/** A tooltip's rows located (their addresses made once a model, `addresses` asked then); unresolved ones asked again. */
function located(ctx: TooltipContext, key: string, addresses: () => readonly string[]): Located {
  let byModel = LocatedRows.get(ctx.model);
  if (!byModel) LocatedRows.set(ctx.model, (byModel = new Map()));
  let l = byModel.get(key);
  if (!l) {
    const list = addresses();
    const locations = list.map((a) => ctx.locate?.(a));
    byModel.set(key, (l = { addresses: list, locations, starts: locations.map(firstOf) }));
  } else
    for (let i = 0; i < l.addresses.length; i++)
      if (!l.locations[i]) {
        const loc = (l.locations[i] = ctx.locate?.(l.addresses[i]));
        l.starts[i] = firstOf(loc);
      }
  return l;
}

const firstOf = (loc: Location | undefined) => loc?.parts[0]?.ranges[0] ?? Number.POSITIVE_INFINITY;

function byAddresses(ctx: TooltipContext, t: ResolvedTooltip, place: TooltipPlace, out: TooltipOut, cap: number) {
  const l = located(ctx, t.address, () => t.addresses!);
  const made = Made.begin();
  let count = 0;
  for (let i = 0; i < l.locations.length; i++) {
    const loc = l.locations[i];
    if (!loc || !covers(ctx.model, loc, place)) continue;
    count++;
    const text = made.n < cap ? constantText(t) : undefined;
    if (text !== undefined) made.add(text, l.starts[i]);
  }
  emit(t, made, count, out);
}

function byTable(ctx: TooltipContext, t: ResolvedTooltip, place: TooltipPlace, out: TooltipOut, cap: number) {
  const rows = t.node.params?.rows as Record<string, Value> | undefined;
  // a table from a source is not read whole yet (`seq.resolve.tooltip.source`)
  if (!rows || t.target === undefined || t.source !== undefined) return;
  // read to the bound (`tooltipRows`, `seq.resolve.tooltip.rows` past it)
  const n = Math.min(tableLength(rows), DefaultSpecBounds.tooltipRows);
  const l = located(ctx, `${t.address}#${t.schema}`, () => {
    const addresses: string[] = [];
    for (let i = 0; i < n; i++) addresses.push(rowAddress(t.target!, t.schema!, (c) => cellOf(rows[c], i)) ?? '');
    return addresses;
  });
  const made = Made.begin();
  let count = 0;
  for (let i = 0; i < n; i++) {
    const loc = l.locations[i];
    if (!loc || !covers(ctx.model, loc, place)) continue;
    count++;
    if (made.n >= cap) continue;
    const text = inlineText(t, rows, i);
    if (text !== undefined) made.add(text, l.starts[i]);
  }
  emit(t, made, count, out);
}

/**
 * A row's address over the table's target: an `address` column relative to it, or an MVS schema's columns as an
 * `(elements …)` form (each schema's columns alone, as MVS reads them; a column the form does not read left out).
 */
export function rowAddress(target: string, schema: TooltipSchema, cell: (column: string) => string | undefined) {
  if (schema === 'ranges') return undefined;
  if (schema === 'address') {
    const a = cell('address');
    return a ? `${target}/${a}` : undefined;
  }
  const cols = MvsSchemaColumns[schema].filter((c) => !Unread.has(c));
  let form = '(elements';
  for (const c of cols) {
    const v = cell(c);
    if (v !== undefined && v !== '') form += ` :${c} ${formValue(v)}`;
  }
  return `${target}/${form})`;
}

/** MVS's columns the `elements` form does not read (no property of that name). */
const Unread = new Set(['atom_id', 'auth_atom_id']);

/** A value in a form: a number or a bare word as it is, else quoted. */
const formValue = (v: string) => (/^[A-Za-z0-9_.+-]+$/.test(v) ? v : JSON.stringify(v));

/** Whether a location covers one of the place's parts (on its system, or carried onto it). */
function covers(model: StaticModel, loc: Location, place: TooltipPlace): boolean {
  for (const part of place.parts)
    for (const lp of loc.parts) {
      const ranges =
        lp.system === part.system
          ? [part.start, part.end]
          : (model
              .project(rangesLocation(part.system, [part.start, part.end]), lp.system)
              .location.parts.find((p) => p.system === lp.system)?.ranges ?? []);
      for (let r = 0; r + 1 < ranges.length; r += 2) {
        const s = ranges[r];
        const e = Math.max(ranges[r + 1], s + 1);
        for (let k = 0; k + 1 < lp.ranges.length; k += 2)
          if (lp.ranges[k] < e && Math.max(lp.ranges[k + 1], lp.ranges[k] + 1) > s) return true;
      }
    }
  return false;
}

// ---- rows and texts ----

/** A lookup's distinct texts and their starts (reused: one tooltip's at a time), kept in position order. */
const Made = {
  n: 0,
  texts: [] as string[],
  starts: [] as number[],
  begin() {
    this.n = 0;
    return this;
  },
  /** Adds a text at its start unless made already, by start (a stable insertion: few, `cap` at most). */
  add(text: string, start: number) {
    for (let k = 0; k < this.n; k++) if (this.texts[k] === text) return;
    let k = this.n++;
    for (; k > 0 && this.starts[k - 1] > start; k--) {
      this.texts[k] = this.texts[k - 1];
      this.starts[k] = this.starts[k - 1];
    }
    this.texts[k] = text;
    this.starts[k] = start;
  },
};

/** A tooltip's matches as rows: one a plain row, several item rows (the rest counted). */
function emit(t: ResolvedTooltip, made: typeof Made, count: number, out: TooltipOut) {
  if (!made.n) return;
  const several = count > 1;
  for (let k = 0; k < made.n; k++) out.rows.push(rowOf(t, made.texts[k], several));
  out.more += Math.max(0, count - made.n);
}

/** Rows kept by tooltip and text (a plain one and an item one): a warm card makes none. Bounded a tooltip. */
const Rows = new WeakMap<ResolvedTooltip, Map<string, TooltipRow>>();
const RowsPerTooltip = 256;

function rowOf(t: ResolvedTooltip, text: string, item: boolean): TooltipRow {
  let byText = Rows.get(t);
  if (!byText) Rows.set(t, (byText = new Map()));
  const key = item ? `\u0001${text}` : text;
  let row = byText.get(key);
  if (row) return row;
  const rich = t.format === 'markdown' ? richOf(text) : undefined;
  const value = rich ? markdownInlineText(rich) : clean(text);
  const base = t.label !== undefined ? { label: t.label, value } : { label: value, value: '' };
  row = { ...base, ...(rich && hasMarks(rich) ? { rich } : {}), ...(item ? { item: true } : {}) };
  if (byText.size >= RowsPerTooltip) byText.delete(byText.keys().next().value!);
  byText.set(key, row);
  return row;
}

function loadingRow(t: ResolvedTooltip): TooltipRow {
  let byText = Rows.get(t);
  if (!byText) Rows.set(t, (byText = new Map()));
  let row = byText.get('\u0002');
  if (!row) byText.set('\u0002', (row = { ...rowOf(t, 'loading', false), loading: true }));
  return row;
}

/** Whether a tree has more than plain text (a mark or a link): its `rich` is worth carrying. */
const hasMarks = (nodes: readonly MarkdownInline[]) => nodes.some((n) => n.kind !== 'text');

/** Texts kept by a table's row, a tooltip's (undefined: none; a few thousand rows a table at most, the oldest let go). */
const TableTexts = new WeakMap<object, Map<ResolvedTooltip, Map<number, string | null>>>();
const TextsPerTable = 4096;

function cachedText(key: object, t: ResolvedTooltip, i: number, make: () => string | undefined): string | undefined {
  let byTip = TableTexts.get(key);
  if (!byTip) TableTexts.set(key, (byTip = new Map()));
  let texts = byTip.get(t);
  if (!texts) byTip.set(t, (texts = new Map()));
  const c = texts.get(i);
  if (c !== undefined) return c ?? undefined;
  const text = make();
  if (texts.size >= TextsPerTable) texts.delete(texts.keys().next().value!);
  texts.set(i, text ?? null);
  return text;
}

/** A list table's row's text (cached by the row). */
const tableText = (t: ResolvedTooltip, table: RangesList, i: number) =>
  cachedText(table, t, i, () => formatRow(t, (name) => listField(table, i, name)));

/** An inline table's row's text (cached by the row). */
const inlineText = (t: ResolvedTooltip, rows: Record<string, Value>, i: number) =>
  cachedText(rows, t, i, () => formatRow(t, (name) => cellOf(rows[name], i)));

/** An address's text (no fields: the template alone, cached). */
const constantText = (t: ResolvedTooltip) => cachedText(t, t, -1, () => formatRow(t, () => undefined));

/** A row's text: the template over its fields (`{}` the tooltip's field), values escaped; undefined when one is missing. */
function formatRow(t: ResolvedTooltip, field: (name: string) => string | undefined): string | undefined {
  const md = t.format === 'markdown';
  const text = t.template!.format((name) => {
    const v = field(name || t.field);
    if (v === undefined || v === '') return undefined;
    // a number stays as it is (its format spec reads it; a number has no mark in it), any other value is escaped
    return md && !Numeric.test(v) ? escapeMarkdown(v) : v;
  });
  if (text === undefined || text === '') return undefined;
  return text.length > MaxText ? `${cutAt(text, MaxText - 1)}…` : text;
}

/** A text's first `n` UTF-16 units, one less when the last would split a surrogate pair (a code point kept whole). */
function cutAt(text: string, n: number): string {
  const c = text.charCodeAt(n - 1);
  return text.slice(0, c >= 0xd800 && c <= 0xdbff ? n - 1 : n);
}

/** A value that is a number as text: formatted by the template's spec, never escaped. */
const Numeric = /^[+-]?(\d+\.?\d*|\.\d+)(e[+-]?\d+)?$/i;

/**
 * Markdown's escapable characters, each backslash-escaped, and `&` as its entity (the parse decodes entities, so a
 * value's `&amp;` stays as written): a value is text (`A*02:01` stays as written).
 */
export function escapeMarkdown(v: string): string {
  return v.replace(/[\\`*_{}[\]()#+\-.!<>|~]/g, '\\$&').replace(/&/g, '&amp;');
}

/** C0 and C1 controls (but tab and newline) and the bidi overrides and isolates: removed from a text. */
const Controls = /[\u0000-\u0008\u000b-\u001f\u007f-\u009f‪-‮⁦-⁩]/g;
const clean = (s: string) => s.replace(Controls, '');

/** Parses, keyed by text (64): a warm card parses nothing. */
const Parsed = new Map<string, readonly MarkdownInline[]>();

/** A text's Markdown tree, cleaned after the parse (entities decoded there) and with Seq*'s links alone kept. */
export function richOf(text: string): readonly MarkdownInline[] {
  let r = Parsed.get(text);
  if (r) return r;
  r = cleanTree(parseMarkdownInline(text));
  if (Parsed.size >= 64) Parsed.delete(Parsed.keys().next().value!);
  Parsed.set(text, r);
  return r;
}

/** A link Seq* acts on: its own link grammar (`!seq-focus=…`), through commands; anything else is its text. */
export const isSeqLink = (href: string) => /^!seq-[a-z][a-z-]*=/.test(href);

function cleanTree(nodes: readonly MarkdownInline[]): MarkdownInline[] {
  const out: MarkdownInline[] = [];
  for (const n of nodes) {
    if (n.kind === 'text' || n.kind === 'code') {
      const text = clean(n.text);
      if (text) out.push({ kind: n.kind, text });
    } else if (n.kind === 'break') out.push(n);
    else if (n.kind === 'link') {
      const children = cleanTree(n.children);
      if (isSeqLink(n.href)) out.push({ kind: 'link', href: n.href, children });
      else out.push(...children);
    } else out.push({ kind: n.kind, children: cleanTree(n.children) });
  }
  return out;
}

// ---- reading rows ----

/** A list row's identity: its id, else its row number. */
function rowKey(t: RangesList, i: number): string | number | undefined {
  return t.id ? (t.id.get(i) ?? undefined) : t.row?.[i];
}

/** A field of a list row as text: a named column, else `start`, `end`, `id`; undefined when it has none. */
function listField(t: RangesList, i: number, name: string): string | undefined {
  const f = t.fields?.find((x) => x.name === name);
  if (f) return cellOf(f.column as unknown as Value, i);
  if (name === 'start') return String(t.start[i]);
  if (name === 'end') return String(t.end[i]);
  if (name === 'id') return t.id?.get(i) ?? undefined;
  return undefined;
}

/** A column's cell as text (a category's, a number's, a list's), undefined when absent. */
function cellOf(c: Value | undefined, i: number): string | undefined {
  if (c instanceof CategoryColumn) return c.get(i) ?? undefined;
  if (Array.isArray(c)) {
    const v = c[i];
    return v === null || v === undefined ? undefined : String(v);
  }
  if (isColumn(c) || ArrayBuffer.isView(c)) {
    const v = (c as unknown as ArrayLike<number>)[i];
    return v === undefined || Number.isNaN(v) ? undefined : String(v);
  }
  return undefined;
}

const meets = (s: number, e: number, start: number, end: number) => s < end && Math.max(e, s + 1) > start;
