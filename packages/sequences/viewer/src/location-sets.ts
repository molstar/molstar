/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Location sets (seqspec:format#sections; seqspec:interaction#channels, seqspec:interaction#results;
// seqspec:viewer#actions, LS1, LS2): a section's `highlights` kept as highlight channels from the shown document;
// bookmarks, an inline `ranges` annotation the bookmark actions edit, stepped through `ctx.results`; and
// `results.save`, a results set kept as a location set. Edge rate throughout.

import { type DocNode, type LineStyle, LineStyles, type Op, sectionHighlights } from '@molstar/seq-core/format';
import { type AnnotationInput, KeyKind, type KeyTarget, type Loci, type SeqSet } from '@molstar/seq-core/model';
import { SeqSelection } from '@molstar/seq-core/model/query';
import { type ActionEntry, nodePath } from './actions.js';
import type { SeqViewerContext } from './context.js';
import type { Assembled } from './graph.js';
import type { HighlightManager } from './highlight.js';
import { isCoordinates } from './projection.js';
import type { Viewports } from './viewport.js';

/**
 * The channels of the shown view's sections' `highlights`: `set:<section>/<n>`, one per entry, lit in its section only;
 * each set again only when its entry or the set changes, cleared when its entry goes, never another owner's channel.
 */
export class SetChannels {
  private readonly highlight: HighlightManager;
  private keys = new Map<string, string>();
  private set: SeqSet | null = null;

  constructor(highlight: HighlightManager) {
    this.highlight = highlight;
  }

  update(a: Assembled | null): void {
    const v = a?.view ?? null;
    const fresh = v?.set !== this.set;
    this.set = v?.set ?? null;
    const next = new Map<string, string>();
    for (let s = 0; v && s < v.resolved.sections.count; s++) {
      const section = v.set.strings.get(v.resolved.sections.ref[s]);
      sectionHighlights(v.resolved, v.set, s).forEach((h, n) => {
        const channel = `set:${section}/${n}`;
        // Every item as written: `null` (all) is not `[]` (none), and an id may hold any character.
        const key = JSON.stringify([h.data, h.items, h.style]);
        next.set(channel, key);
        if (fresh || this.keys.get(channel) !== key) this.highlight.set(channel, h.loci, { style: h.style, section });
      });
    }
    for (const c of this.keys.keys()) if (!next.has(c)) this.highlight.clear(c);
    this.keys = next;
  }
}

// ---- Bookmarks ------------------------------------------------------------------------------------------------------

/**
 * `bookmark.add`'s target: what to keep (loci of the shown set), its label, and the set (`<section>.bookmarks`). A
 * `color` is written as the row's `color` column, which a layer shows through `field_color`, and a `line_style` (one
 * of the format's `LineStyles`) as its `line_style` column, which a layer shows through `field_line_style`; left out,
 * no column is written (an old set's rows stay as they are).
 */
export interface BookmarkAdd {
  readonly section: string;
  readonly loci: Loci;
  readonly label?: string;
  readonly color?: string;
  readonly line_style?: LineStyle;
  readonly data?: string;
}
export interface BookmarkItem {
  readonly section: string;
  readonly item: string;
  readonly data?: string;
}
export interface BookmarkStep {
  readonly section: string;
  readonly dir: 1 | -1;
  readonly data?: string;
}
type Row = {
  target: string;
  start: number;
  stop: number;
  id: string;
  label?: string;
  color?: string;
  line_style?: LineStyle;
};
type Run = Parameters<ActionEntry['run']>[2];

const t: KeyTarget = { kind: -1, id: -1, row: -1 };

const setOf = (ctx: SeqViewerContext) => ctx.graph.assembled.value?.set ?? null;
const nodeOf = (ctx: SeqViewerContext, ref: string) => nodePath(ctx.document.shown.value, ref)?.at(-1) ?? null;
const isSection = (ctx: SeqViewerContext, ref: unknown) =>
  typeof ref === 'string' && nodeOf(ctx, ref)?.kind === 'section';
/** The bookmark set's ref: `data`, else `<section>.bookmarks`. */
export const bookmarksOf = (section: string, data?: string) => data ?? `${section}.bookmarks`;
/** An inline rows table of a `ranges` annotation (what the actions write), or null: not one the actions can edit. */
function rowsOf(n: DocNode | null): Row[] | null {
  const p = n?.params;
  return n?.kind === 'annotation' && p?.geometry === 'ranges' && Array.isArray(p.value) ? (p.value as Row[]) : null;
}
/** `toRangesInput`'s columns as the rows form a document keeps (no `score`; a null label left out). */
function rowsFrom(a: AnnotationInput): Row[] {
  const c = a.columns as Record<string, readonly unknown[]>;
  return c.target.map((target, k) => ({
    target: target as string,
    start: c.start[k] as number,
    stop: c.stop[k] as number,
    id: c.id[k] as string,
    ...(c.label[k] != null ? { label: c.label[k] as string } : {}),
  }));
}
/** The item ids of the set's annotation `ref`, in item order; empty when it is not in the set. */
function itemIds(set: SeqSet | null, ref: string): string[] {
  if (!set?.byName(ref, t) || t.kind !== KeyKind.annotation) return [];
  const ann = set.annotations[t.id];
  const out: string[] = [];
  for (let k = 0; ann && k < ann.itemCount; k++) if (ann.itemKey[k] >= 0) out.push(set.strings.get(ann.itemKey[k]));
  return out;
}
/** A patch inserting a ranges annotation under `data` and a `highlights` entry for it on the section. */
function insertSet(ctx: SeqViewerContext, section: string, ref: string, label: string, value: Row[]): Op[] {
  const old = (nodeOf(ctx, section)?.params?.highlights as unknown[] | undefined) ?? [];
  return [
    { op: 'insert', parent: 'data', node: { kind: 'annotation', ref, params: { geometry: 'ranges', label, value } } },
    { op: 'set', ref: section, params: { highlights: [...old, { data: ref }] } },
  ];
}

/** The results set a section's bookmarks step through (`ctx.results`; lit by their band, not again). */
const stepName = (section: string) => `bookmarks:${section}`;
/** The section's bookmarks as that results set, set again when its items changed; their ids. */
function bookmarkResults(ctx: SeqViewerContext, section: string, data?: string): string[] {
  const ref = bookmarksOf(section, data);
  const ids = itemIds(setOf(ctx), ref).filter((id) => ctx.refs.toLoci(ref, { items: [id] }));
  const st = ctx.results.sets.value.get(stepName(section));
  if (!st || st.labels.join('\n') !== ids.join('\n')) {
    const loci = ids.map((id) => ctx.refs.toLoci(ref, { items: [id] }) as Loci);
    ctx.results.set(stepName(section), loci, { labels: ids, section, light: false });
  }
  return ids;
}

/** The bookmark of `section` stepped or gone to last; null: none yet. */
export function bookmarkAt(ctx: SeqViewerContext, section: string): string | null {
  const st = ctx.results.sets.value.get(stepName(section));
  return st && st.index >= 0 ? (st.labels[st.index] ?? null) : null;
}

function add(ctx: SeqViewerContext, a: BookmarkAdd, origin: Run) {
  const set = setOf(ctx) as SeqSet;
  const ref = bookmarksOf(a.section, a.data);
  const node = nodeOf(ctx, ref);
  const old = rowsOf(node) ?? [];
  let n = 1;
  for (const r of old) n = Math.max(n, Number(/^b(\d+)$/.exec(String(r.id))?.[1] ?? 0) + 1);
  const sel = SeqSelection.union(SeqSelection.fromLoci(set, [a.loci]));
  const kept = rowsFrom(SeqSelection.toRangesInput(sel, ref, { ids: [`b${n}`], labels: [a.label ?? null] }));
  const rows = kept.map((r) => ({
    ...r,
    ...(a.color ? { color: a.color } : {}),
    ...(a.line_style ? { line_style: a.line_style } : {}),
  }));
  const ops: Op[] = node
    ? [{ op: 'set', ref, params: { value: [...old, ...rows] } }]
    : insertSet(ctx, a.section, ref, 'Bookmarks', rows);
  void ctx.document.propose({ format: 'seq_view_spec_patch', ops }, 'Add bookmark', origin);
}

function step(ctx: SeqViewerContext, s: BookmarkStep, origin: Run) {
  const at = bookmarkAt(ctx, s.section);
  const ids = bookmarkResults(ctx, s.section, s.data);
  const i = ids.indexOf(at ?? '');
  const next = i < 0 ? (s.dir > 0 ? 0 : ids.length - 1) : (i + s.dir + ids.length) % ids.length;
  ctx.results.go(stepName(s.section), next, { origin });
}
/** One bookmark focused in its section and accented there; the current one for `bookmark.step`. */
function go(ctx: SeqViewerContext, b: BookmarkItem, origin: Run) {
  const i = bookmarkResults(ctx, b.section, b.data).indexOf(b.item);
  if (i >= 0) ctx.results.go(stepName(b.section), i, { origin });
}

/** `results.save`'s target: the results set `name` saved into a section as a location set. */
export interface ResultsSave {
  readonly section: string;
  readonly name: string;
  /**
   * The annotation's ref. Given, it names the kept set: saving again under it replaces that set's rows and label (a
   * query kept twice is one set). Default `<section>.<name>`, numbered on when taken.
   */
  readonly ref?: string;
  readonly label?: string;
}
/** One item per result (`r<k>`, labeled as the results are), by value, and a `highlights` entry. */
function save(ctx: SeqViewerContext, a: ResultsSave, origin: Run) {
  const set = setOf(ctx) as SeqSet;
  const base = a.ref ?? `${a.section}.${a.name.replace(/[^A-Za-z0-9_.-]/g, '-')}`;
  // A ref the caller names is replaced when it is an inline set of ours; otherwise (and by default) one of its own.
  const replaced = a.ref !== undefined && rowsOf(nodeOf(ctx, a.ref)) !== null;
  let ref = base;
  if (!replaced) for (let k = 2; nodeOf(ctx, ref) || set.byName(ref, t); k++) ref = `${base}-${k}`;
  const projector = (ctx.viewports as Viewports).projector;
  const labels = ctx.results.sets.value.get(a.name)?.labels ?? [];
  const rows = ctx.results.results(a.name).flatMap((r, k) => {
    const loci = isCoordinates(r) ? projector.lociOf(r) : [r];
    const sel = SeqSelection.union(SeqSelection.fromLoci(set, loci));
    const opts = { ids: [`r${k + 1}`], labels: [labels[k] ?? null] };
    return SeqSelection.isEmpty(sel) ? [] : rowsFrom(SeqSelection.toRangesInput(sel, ref, opts));
  });
  const label = a.label ?? a.name;
  const ops: Op[] = replaced
    ? [{ op: 'set', ref, params: { label, value: rows } }]
    : insertSet(ctx, a.section, ref, label, rows);
  void ctx.document.propose(
    { format: 'seq_view_spec_patch', ops },
    replaced ? 'Keep results again' : 'Keep results',
    origin,
  );
}

/** `highlight.remove`'s target: entries of a section's `highlights`, by their data refs (several: one patch). */
export interface HighlightRemove {
  readonly section: string;
  readonly data: string | readonly string[];
}
const datasOf = (a: HighlightRemove): readonly string[] => (typeof a.data === 'string' ? [a.data] : a.data);
const highlightsOf = (ctx: SeqViewerContext, section: string) =>
  (nodeOf(ctx, section)?.params?.highlights as { data: string }[] | undefined) ?? [];
/**
 * The entry out of the section's `highlights`; a set the section holds inline (`<section>.…`, a kept result set or its
 * bookmarks) goes with it, one patch, so Undo brings both back. Other data stays: the entry only shows it.
 */
function removeHighlight(ctx: SeqViewerContext, a: HighlightRemove, origin: Run) {
  const datas = datasOf(a);
  const rest = highlightsOf(ctx, a.section).filter((h) => !datas.includes(h.data));
  const ops: Op[] = [{ op: 'set', ref: a.section, params: { highlights: rest } }];
  for (const data of datas)
    if (data.startsWith(`${a.section}.`) && rowsOf(nodeOf(ctx, data)) !== null && !rest.some((h) => h.data === data))
      ops.push({ op: 'remove', ref: data });
  void ctx.document.propose({ format: 'seq_view_spec_patch', ops }, 'Remove highlight', origin);
}

const entry = <T>(name: string, title: string, when: ActionEntry<T>['when'], run: ActionEntry<T>['run']) =>
  ({ name, title, when, run }) as ActionEntry;

/** `bookmark.add`, `bookmark.remove`, `bookmark.step`, `bookmark.go`, `highlight.remove` and `results.save`. */
export const BookmarkActions: readonly ActionEntry[] = [
  entry<BookmarkAdd>(
    'bookmark.add',
    'Bookmark',
    (ctx, a) => {
      const set = setOf(ctx);
      if (!a || !isSection(ctx, a.section) || !set || a.loci?.set !== set) return false;
      if (a.color !== undefined && typeof a.color !== 'string') return false;
      if (a.line_style !== undefined && !LineStyles.includes(a.line_style)) return false;
      const node = nodeOf(ctx, bookmarksOf(a.section, a.data));
      return (node === null || rowsOf(node) !== null) && !SeqSelection.isEmpty(SeqSelection.fromLoci(set, [a.loci]));
    },
    add,
  ),
  entry<BookmarkItem>(
    'bookmark.remove',
    'Remove bookmark',
    (ctx, a) =>
      !!a &&
      isSection(ctx, a.section) &&
      !!rowsOf(nodeOf(ctx, bookmarksOf(a.section, a.data)))?.some((r) => r.id === a.item),
    (ctx, a, origin) => {
      const ref = bookmarksOf(a.section, a.data);
      const value = (rowsOf(nodeOf(ctx, ref)) ?? []).filter((r) => r.id !== a.item);
      void ctx.document.propose(
        { format: 'seq_view_spec_patch', ops: [{ op: 'set', ref, params: { value } }] },
        'Remove bookmark',
        origin,
      );
    },
  ),
  entry<BookmarkStep>(
    'bookmark.step',
    'Next bookmark',
    (ctx, s) =>
      !!s &&
      (s.dir === 1 || s.dir === -1) &&
      isSection(ctx, s.section) &&
      itemIds(setOf(ctx), bookmarksOf(s.section, s.data)).length > 0,
    step,
  ),
  entry<BookmarkItem>(
    'bookmark.go',
    'Go to bookmark',
    (ctx, b) =>
      !!b && isSection(ctx, b.section) && itemIds(setOf(ctx), bookmarksOf(b.section, b.data)).includes(b.item),
    go,
  ),
  entry<HighlightRemove>(
    'highlight.remove',
    'Remove highlight',
    (ctx, a) =>
      !!a && isSection(ctx, a.section) && highlightsOf(ctx, a.section).some((h) => datasOf(a).includes(h.data)),
    removeHighlight,
  ),
  entry<ResultsSave>(
    'results.save',
    'Keep results',
    (ctx, a) => !!a && isSection(ctx, a.section) && !!setOf(ctx) && ctx.results.results(a.name).length > 0,
    save,
  ),
];

/** A section's bookmarks' item ids, in item order (`data`: another set than `<section>.bookmarks`). */
export const bookmarkIds = (ctx: SeqViewerContext, section: string, data?: string): string[] =>
  itemIds(setOf(ctx), bookmarksOf(section, data));
