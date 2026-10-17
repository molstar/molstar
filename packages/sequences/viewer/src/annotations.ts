/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the `Annotations` part shows and does (seqspec:interaction#the-default-ui; DOM: dom/annotations.ts, React:
// react/annotations.tsx): a section's bookmarks set (`bookmarksOf`), or several (the `places` a person chooses between,
// each a set, e.g. a track above the letters and one below), as a list of colored, labeled ranges. A thin client
// of the patch API: Add, Delete and the list are the bookmark actions (`bookmark.add`, `bookmark.remove`), which
// propose a `seq_view_spec_patch` (one undo step each); a row reveals its range (`ctx.reveal`). Edge-rate.

import { type LineStyle, LineStyles } from '@molstar/seq-core/format';
import type { Loci } from '@molstar/seq-core/model';
import { nodePath } from './actions.js';
import type { SeqViewerContext } from './context.js';
import { bookmarksOf } from './location-sets.js';

/** The color a new annotation starts from (the part's swatch). */
export const ANNOTATION_COLOR = '#e76f51';

/** The line style a new annotation starts from: the plain bar or underline. */
export const ANNOTATION_STYLE: LineStyle = 'solid';
/** The styles the picker offers (the format's `LineStyles`). */
export const ANNOTATION_STYLES: readonly LineStyle[] = LineStyles;

/** The glyphs the picker offers a place of markers (`kind` `points`), and the one a new marker starts from. */
export const ANNOTATION_GLYPHS: readonly string[] = ['circle', 'triangle', 'star', 'diamond', 'square'];
export const ANNOTATION_GLYPH = 'circle';

/**
 * A set the panel can add to and lists: its ref (`bookmark.add`'s `data`), the name a person sees, and what its track
 * draws (`kind`): `ranges` (the default) and `frame` take a line style, `points` a glyph (seqspec:interaction#the-default-ui).
 */
export interface AnnotationPlace {
  readonly data: string;
  readonly label: string;
  readonly kind?: 'ranges' | 'points' | 'frame';
}
/** Whether a place's rows are markers (a glyph each, no line style). */
export const placeGlyphs = (place: AnnotationPlace | undefined) => place?.kind === 'points';

/** The sets a panel edits: its `places` when given, else the one set `data` (default `<section>.bookmarks`). */
export const annotationPlaces = (
  section: string,
  data?: string,
  places?: readonly AnnotationPlace[],
): readonly AnnotationPlace[] =>
  places && places.length > 0 ? places : [{ data: bookmarksOf(section, data), label: '' }];

/** An annotation as listed: one item of a set (a row, or several runs of one id). */
export interface AnnotationEntry {
  /** Unique in the list: the id, and with several sets, the set's ref before it. */
  readonly key: string;
  /** The set it is in (`bookmark.remove`'s `data`), and its name: empty when the panel has one set. */
  readonly data: string;
  readonly place: string;
  readonly id: string;
  /** The row's label, else its id. */
  readonly label: string;
  /** The row's `color` column, else the part's default. */
  readonly color: string;
  /** The row's `line_style` column, else `solid`. */
  readonly lineStyle: LineStyle;
  /** The row's `glyph` column (a marker's), else "". */
  readonly glyph: string;
  /** Its positions as people count them (1-based, inclusive), and `+n` for further runs. */
  readonly span: string;
}

/** The states a list of entries and the add control are derived from. */
export const annotationStates = (ctx: SeqViewerContext) => [ctx.document.shown, ctx.selection.current] as const;

type Row = {
  start?: unknown;
  stop?: unknown;
  id?: unknown;
  label?: unknown;
  color?: unknown;
  line_style?: unknown;
  glyph?: unknown;
};

/** One set's items in item order, from the shown document's inline rows; empty when it is not a set the actions edit. */
function entriesOf(ctx: SeqViewerContext, place: AnnotationPlace, keyed: boolean): AnnotationEntry[] {
  const node = nodePath(ctx.document.shown.value, place.data)?.at(-1);
  const rows = node?.kind === 'annotation' && node.params?.geometry === 'ranges' ? node.params.value : null;
  if (!Array.isArray(rows)) return [];
  const by = new Map<string, Row[]>();
  for (const r of rows as Row[]) {
    const id = String(r.id);
    by.set(id, [...(by.get(id) ?? []), r]);
  }
  return [...by].map(([id, runs]) => {
    const first = runs[0];
    const [a, b] = [Number(first.start) + 1, Number(first.stop)];
    return {
      key: keyed ? `${place.data}:${id}` : id,
      data: place.data,
      place: place.label,
      id,
      label: typeof first.label === 'string' && first.label ? first.label : id,
      color: typeof first.color === 'string' && first.color ? first.color : ANNOTATION_COLOR,
      lineStyle:
        typeof first.line_style === 'string' && LineStyles.includes(first.line_style as LineStyle)
          ? (first.line_style as LineStyle)
          : ANNOTATION_STYLE,
      glyph: typeof first.glyph === 'string' ? first.glyph : '',
      span: `${a === b ? a : `${a}–${b}`}${runs.length > 1 ? ` +${runs.length - 1}` : ''}`,
    };
  });
}

/** The sets' items, a set after another in `places` order. */
export const annotationEntries = (ctx: SeqViewerContext, places: readonly AnnotationPlace[]): AnnotationEntry[] =>
  places.flatMap((p) => entriesOf(ctx, p, places.length > 1));

/** A row's accessible names: its place is part of them when the panel has several. */
export const showName = (a: AnnotationEntry) => `Show ${a.label} (${a.span}${a.place ? `, ${a.place}` : ''})`;
export const deleteName = (a: AnnotationEntry) => `Delete ${a.label}${a.place ? ` (${a.place})` : ''}`;

/** The selection's loci an annotation can keep (the sequence's first), or null. */
export function selectedLoci(ctx: SeqViewerContext): Loci | null {
  const sel = ctx.selection.current.value;
  return (sel.seq ?? sel.items ?? sel.rows ?? sel.columns ?? null) as Loci | null;
}

/**
 * What `bookmark.add` takes for the selection and the person's label, color and line style (`solid` writes no column:
 * the layer draws it plain), or a marker's glyph (written as its `glyph` column); null: nothing to add now.
 */
export function annotationTarget(
  ctx: SeqViewerContext,
  section: string,
  label: string,
  color: string,
  data?: string,
  lineStyle: LineStyle = ANNOTATION_STYLE,
  glyph?: string,
) {
  const loci = selectedLoci(ctx);
  const line_style = lineStyle === ANNOTATION_STYLE || glyph ? undefined : lineStyle;
  const target = { section, loci: loci as Loci, label: label.trim() || undefined, color, data, line_style, glyph };
  return loci && ctx.actions.can('bookmark.add', target) ? target : null;
}

/** Adds the selection; false: nothing to add (no selection, or the set cannot take it). */
export function addAnnotation(
  ctx: SeqViewerContext,
  section: string,
  label: string,
  color: string,
  data?: string,
  lineStyle?: LineStyle,
  glyph?: string,
) {
  const target = annotationTarget(ctx, section, label, color, data, lineStyle, glyph);
  return target !== null && ctx.actions.run('bookmark.add', target);
}
/** Deletes one (a patch; Undo brings it back). */
export const removeAnnotation = (ctx: SeqViewerContext, section: string, id: string, data?: string) =>
  ctx.actions.run('bookmark.remove', { section, item: id, data });
/** Brings one into view in every section that shows it, changing no selection, focus or zoom. */
export function revealAnnotation(ctx: SeqViewerContext, section: string, id: string, data?: string) {
  const loci = ctx.refs.toLoci(bookmarksOf(section, data), { items: [id] });
  return loci !== null && ctx.reveal(loci);
}
