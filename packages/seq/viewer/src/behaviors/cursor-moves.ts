/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `cursor` behavior's moves (_spec/README.md "Behaviors", Cursor): internal, not in the behaviors' entry; the
// provider gives them the last frame's layout (`cursorBehavior` in `cursor.ts` is the public face, without rows).
//
// The keyboard cursor: the step keys (`cursor.keys`) move it a position along its section's axis in display order, or
// a line up or down in a layout of lines (a text: the section's layout's step), else to the section above or below
// (its position carried along the view's link); Shift with a step extends the range the anchor is in (from the cursor's
// cell when there is no anchor), up and down by a line; Cmd/Ctrl+Shift with a step starts a new range at the cursor,
// the others kept; Space toggles the cursor's cell; Enter focuses the selection, else the cursor's cell. A cursor moved
// out of view brings its section to it by the layout's least move (`reveal.cursor`). Every move is a command.
//
// On a section of two axes the cursor has a row too: on an alignment's rows ↑ ↓ step display rows (the pinned band to
// the scrolled rows shown, a summary row one stop), Page Up and Page Down a page of rows with the band, Cmd/Ctrl+↑ ↓ the
// first and the last row, ← from the first column onto the row's header and → back; Shift with an arrow extends the
// block, Cmd/Ctrl+Shift starts one (the others kept); Space toggles the cell, on the header the row's clade; Enter
// opens a summary row, on the header focuses the row. Moves are arithmetic on display rows (a row never loaded is
// stepped over by its index); the band scrolls to show the row by the least move. On a grid ↑ ↓ step its y and
// Cmd/Ctrl+↑ ↓ go to its ends, its y shown by the least move.

import {
  EmptyKeys,
  EmptyLocation,
  type Location,
  location,
  rangesLocation,
  unionLocations,
} from '@molstar/seqx-core/query';
import { collapsedSection, type GridLay, type RowsLay, type ViewLayout } from '@molstar/seqx-core/render';
import { type Axis, projectOntoView } from '@molstar/seqx-core/view';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import type { ViewerConfig } from '@molstar/seqx-viewer/config';
import {
  type AddressModel,
  axisBreaks,
  axisTotal,
  type Cursor,
  cellOffset,
  cellsRect,
  elementAt,
  HiddenMessage,
  hiddenLength,
  type SelectionModel,
  type SpecModel,
  type ViewModel,
} from '@molstar/seqx-viewer/models';
import type { CursorBehavior, CursorContext } from './cursor.js';
import { bandFolds, pagedTop, pageRows, revealedTop, setBandFolds } from './rows-band.js';
import { type RowsOf, rowsOf as rowsIn, twoAxesLayer } from './two-axes.js';

/** The cursor less its row: a cell on the section's axis. */
const plain = (section: string, system: string, position: number): Cursor => ({ section, system, position });

/** The cursor on its row's cell (off the header). */
const onCell = (c: Cursor): Cursor => {
  const { part: _, ...rest } = c;
  return rest;
};

/** A gap byte (`-`, `.`, none): the rows layer's (`rows-model.ts` `isGap`, which the render entry does not export). */
const gap = (k: number) => k === 45 || k === 46 || k === 0;

export function createCursor(
  addresses: AddressModel,
  selection: SelectionModel,
  view: ViewModel,
  spec: SpecModel,
  config: ViewerConfig,
  layout: () => ViewLayout | undefined,
): CursorBehavior & { dispose(): void } {
  const keys = config.cursor.keys;
  const sectionOf = (address: string | undefined) => view.resolved.value?.sections.find((s) => s.address === address);
  const rowsOf = (address: string): RowsOf | undefined => rowsIn(layout(), address);
  /** A grid section's y axis (its own system, in display order). */
  const gridY = (address: string): Axis | undefined => sectionOf(address)?.y?.axis;
  /** The band's top now: the model's (the last frame may not have drawn a move just made), else the layout's. */
  const topOf = (r: RowsOf) => view.band(r.track)?.bandTop ?? r.lay.top;

  /**
   * The display row each cursor stands on, by the cursor (a pinned row is in the band and in its place: the copy a
   * step or a click left it on); a cursor with none stands on its row's own (`RowsLay.cursorRow`).
   */
  const standsOn = new WeakMap<Cursor, number>();
  /** Whether display row `i` shows data row `y` (its own row, or a summary row holding it). */
  const shows = (r: RowsOf, i: number, y: number) => {
    const d = r.lay.display;
    if (i < 0 || i >= d.count) return false;
    const e = d.rowAt(i);
    return e === y || (e < 0 && !!d.members.get(i)?.includes(y));
  };
  // a click sets the cursor: it stands on the copy clicked (the rows' last hit)
  const stopClicks = selection.cursor.listen((c) => {
    if (!c || c.yPosition === undefined || standsOn.has(c)) return;
    const r = rowsOf(c.section);
    if (r && c.ySystem === r.lay.rowsSystem && shows(r, r.lay.hitRow, c.yPosition)) standsOn.set(c, r.lay.hitRow);
  });
  /** The display row the cursor stands on: where a step or a click left it, else its row's own (or the nearest). */
  const rowIndex = (r: RowsOf, c: Cursor): number => {
    const y = c.yPosition;
    if (y === undefined || c.ySystem !== r.lay.rowsSystem || !r.lay.display.count) return -1;
    const i = standsOn.get(c);
    return i !== undefined && shows(r, i, y) ? i : r.lay.cursorRow(y);
  };
  /** The cursor on display row `i` (a summary row by its first shown row), its part kept. */
  const atRow = (c: Cursor, r: RowsOf, i: number): Cursor | undefined => {
    const d = r.lay.display;
    if (i < 0 || i >= d.count) return undefined;
    let e = d.rowAt(i);
    if (e < 0) e = d.members.get(i)?.[0] ?? -1;
    if (e < 0) return undefined;
    const next = { ...c, ySystem: r.lay.rowsSystem, yPosition: e };
    standsOn.set(next, i);
    return next;
  };
  /** The first scrolled row shown whole (the band's top), or the last shown whole. */
  const shownRow = (r: RowsOf, which: 'first' | 'last'): number => {
    const l = r.lay;
    const d = l.display;
    const top = topOf(r);
    if (which === 'first') return Math.min(d.count - 1, d.pinned + Math.ceil(top));
    const area = Math.max(0, l.box.height - d.pinned * l.rowScale) / l.rowScale;
    return Math.max(Math.min(d.pinned, d.count - 1), Math.min(d.count - 1, d.pinned + Math.floor(top + area) - 1));
  };
  /** The cursor placed on a section's rows or grid y (entering it): its first row shown, or its last. */
  const enter = (c: Cursor, which: 'first' | 'last'): Cursor => {
    const r = rowsOf(c.section);
    if (r && r.lay.display.count) {
      // the first row shown: the pinned band's first when there is one
      const i = which === 'first' && r.lay.display.pinned ? 0 : shownRow(r, which);
      return atRow(c, r, i) ?? c;
    }
    const y = gridY(c.section);
    const total = y && axisTotal(y);
    if (!y || !total) return c;
    const v = view.get(c.section)?.y;
    const o = which === 'first' ? Math.ceil(v?.start ?? 0) : Math.floor((v?.end ?? total) - 1);
    const p = elementAt(y, Math.max(0, Math.min(total - 1, o)));
    return p === undefined ? c : { ...c, ySystem: y.system, yPosition: p };
  };
  /** Whether a section has a second axis the cursor can stand on (rows, a grid's y). */
  const hasRows = (address: string) => !!rowsOf(address) || !!gridY(address);

  /** The cursor now, or one placed at the first position the section shows (its first row shown too). */
  const current = (ctx: CursorContext): Cursor | undefined => {
    const c = selection.cursor.value;
    if (c && sectionOf(c.section)) return c;
    const section = sectionOf(ctx.section) ?? view.resolved.value?.sections.find((s) => !collapsedSection(s));
    const v = section && view.get(section.address);
    if (!section || !v || v.pending) return undefined;
    const p = elementAt(section.axis, Math.ceil(v.start));
    return p === undefined
      ? undefined
      : enter({ section: section.address, system: section.axis.system, position: p }, 'first');
  };

  /**
   * A section's viewport pans to show its offset `o` when it is out of view. A follower's does not (it would move its
   * leader): the cursor on a follower keeps the rule it had before followers led, revealing only once detached.
   */
  const reveal = (section: string, o: number) => {
    if (!config.reveal.cursor || !view.canMove(section)) return;
    if (view.isFollowing(section) && config.follow.followerMoves !== 'detach') return;
    const v = view.get(section);
    if (!v || v.pending) return;
    view.reveal(section, o, o + 1);
  };
  /** A rows band scrolls to show the cursor's row by the least move (from `top`, the band's now when not given). */
  const revealRow = (c: Cursor, top?: number) => {
    const r = rowsOf(c.section);
    if (!r || !config.reveal.cursor) return;
    const i = rowIndex(r, c);
    if (i < 0) return;
    const next = revealedTop(r.lay, top ?? topOf(r), i);
    if (next !== topOf(r)) view.setBand(r.track, { bandTop: next });
  };
  /**
   * A grid's y shows the cursor's row by the least move: a y of its own (or one following a section) pans by the px the
   * cell is past the plot's edge, a y that shows x's window reveals the offset as x does.
   */
  const revealY = (c: Cursor) => {
    const y = gridY(c.section);
    if (!y || !config.reveal.cursor || c.yPosition === undefined || c.ySystem !== y.system) return;
    if (!view.movesY(c.section)) {
      const o = cellOffset(y, c.yPosition);
      if (o !== undefined) reveal(c.section, o);
      return;
    }
    const s = layout()?.sections.find((x) => x.section.address === c.section);
    const g = s && (twoAxesLayer(s, c.ySystem)?.lay as GridLay | undefined);
    const ay = g?.axisY;
    if (!ay) return;
    const a = ay.x(c.yPosition);
    const b = ay.x(c.yPosition + 1);
    if (!Number.isFinite(a) || !Number.isFinite(b)) return;
    const lo = Math.min(a, b);
    const hi = Math.max(a, b);
    const end = ay.x0 + ay.width;
    const dy = lo < ay.x0 ? lo - ay.x0 : hi > end ? hi - end : 0;
    if (dy) view.panPx(c.section, 0, dy);
  };

  /** The cursor a line up or down within its section (a layout of lines), else undefined. */
  const line = (c: Cursor, dir: 'up' | 'down'): Cursor | undefined => {
    const axis = sectionOf(c.section)?.axis;
    const o = axis && cellOffset(axis, c.position);
    if (!axis || o === undefined) return undefined;
    const next = view.step(c.section, o, dir);
    const p = next < 0 ? undefined : elementAt(axis, next);
    return p === undefined ? undefined : { ...c, position: p };
  };

  /**
   * The cursor `delta` rows on its section's rows (display rows) or its grid's y (display order), clamped to their ends
   * (`clamp`), else undefined past them; undefined on a section of one axis or a cursor with no row there. A step down
   * from the pinned band's last row lands on the first scrolled row shown (the band keeps its scroll).
   */
  const rowStep = (c: Cursor, delta: number, clamp: boolean): Cursor | undefined => {
    const r = rowsOf(c.section);
    if (r && c.ySystem === r.lay.rowsSystem) {
      const d = r.lay.display;
      const i = rowIndex(r, c);
      if (i < 0) return undefined;
      let j = i + delta;
      if (delta === 1 && i === d.pinned - 1 && d.pinned < d.count) j = shownRow(r, 'first');
      return atRow(c, r, clamp ? Math.max(0, Math.min(d.count - 1, j)) : j);
    }
    const y = gridY(c.section);
    const total = y && axisTotal(y);
    const o = y && c.yPosition !== undefined && c.ySystem === y.system ? cellOffset(y, c.yPosition) : undefined;
    if (!y || !total || o === undefined) return undefined;
    const next = clamp ? Math.max(0, Math.min(total - 1, o + delta)) : o + delta;
    const p = next < 0 || next >= total ? undefined : elementAt(y, next);
    return p === undefined ? undefined : { ...c, yPosition: p };
  };
  /** The cursor on the first (`-1`) or the last (`1`) row of its rows or its grid's y. */
  const rowEnd = (c: Cursor, end: -1 | 1): Cursor | undefined => {
    const r = rowsOf(c.section);
    if (r && c.ySystem === r.lay.rowsSystem) return atRow(c, r, end < 0 ? 0 : r.lay.display.count - 1);
    const y = gridY(c.section);
    const total = y && axisTotal(y);
    if (!y || !total || c.ySystem !== y.system) return undefined;
    const p = elementAt(y, end < 0 ? 0 : total - 1);
    return p === undefined ? undefined : { ...c, yPosition: p };
  };
  /** Whether the cursor stands on a row of its section (rows or a grid's y). */
  const onRows = (c: Cursor | undefined): boolean => {
    if (!c || c.ySystem === undefined || c.yPosition === undefined) return false;
    const r = rowsOf(c.section);
    return r ? c.ySystem === r.lay.rowsSystem : c.ySystem === gridY(c.section)?.system;
  };
  /**
   * ↑ ↓: a row on rows (past the ends: the section above or below); from a section's other tracks (no row yet) ↓ onto
   * its rows; else a line, else the section above or below.
   */
  const vertical = (c: Cursor, delta: -1 | 1): Cursor | undefined => {
    if (onRows(c)) return rowStep(c, delta, false) ?? across(c, delta);
    if (delta > 0 && hasRows(c.section)) return enter(c, 'first');
    return line(c, delta < 0 ? 'up' : 'down') ?? across(c, delta);
  };

  /**
   * The cursor a step along its axis, in display order (whatever the direction), clamped to the axis; a step from one
   * piece to the next goes over the break between them (_spec/README.md "Behaviors", Cursor).
   */
  const stepped = (c: Cursor, delta: number): Cursor | undefined => {
    const axis = sectionOf(c.section)?.axis;
    // a section of blocks on several systems: its units are no system's positions (stepping across blocks is open)
    if (axis?.pieces === 'blocks' && axis.system !== c.system) return undefined;
    const total = axis && axisTotal(axis);
    const o = axis && cellOffset(axis, c.position);
    if (!axis || total === undefined || o === undefined) return undefined;
    const next = Math.min(Math.max(o + delta, 0), total - 1);
    const p = elementAt(axis, next);
    return p === undefined ? undefined : { ...c, position: p };
  };
  /** ← on rows: from the first column onto the row's header (a step back elsewhere); on the header it stays. */
  const back = (c: Cursor): Cursor | undefined => {
    if (c.part === 'header') return c;
    const s = stepped(c, -1);
    if (s && s.position === c.position && rowsOf(c.section) && onRows(c)) return { ...c, part: 'header' };
    return s;
  };
  /** → from the header back onto the row's cell (the column it left); a step on elsewhere. */
  const on = (c: Cursor): Cursor | undefined => (c.part === 'header' ? onCell(c) : stepped(c, 1));

  /**
   * What a step from `a` to `b` crossed: the hidden stretch of the break between them, named for the announcement
   * ("crossed a break, 1,240 positions hidden"); empty when it crossed none.
   */
  const crossed = (a: Cursor, b: Cursor): string => {
    const axis = sectionOf(b.section)?.axis;
    if (!axis || a.section !== b.section) return '';
    const oa = cellOffset(axis, a.position);
    const ob = cellOffset(axis, b.position);
    if (oa === undefined || ob === undefined) return '';
    const breaks = axisBreaks(axis);
    let hidden = 0;
    let n = 0;
    for (let i = 0; i + 2 < breaks.length; i += 3) {
      const at = breaks[i];
      if (at > Math.min(oa, ob) && at <= Math.max(oa, ob)) {
        hidden += hiddenLength(breaks[i + 1], breaks[i + 2], axis.length);
        n++;
      }
    }
    if (!n) return '';
    const what = n === 1 ? 'a break' : `${n} breaks`;
    return `, crossed ${what}, ${hidden.toLocaleString('en-US')} positions hidden`;
  };

  /** The cursor shown there, its cell brought into view (its column unless on the header, its row or its y). */
  const moveTo = (c: Cursor, top?: number) => {
    selection.setCursor({ ...c, keyboard: true });
    const now = selection.cursor.value;
    const i = standsOn.get(c);
    if (now && i !== undefined) standsOn.set(now, i);
    if (c.part !== 'header') {
      const axis = sectionOf(c.section)?.axis;
      const o = axis && cellOffset(axis, c.position);
      if (o !== undefined) reveal(c.section, o);
    }
    revealRow(c, top);
    revealY(c);
  };

  /** The cursor on the section `delta` above or below, its position carried along the view's link. */
  const across = (c: Cursor, delta: number): Cursor | undefined => {
    const resolved = view.resolved.value;
    if (!resolved) return undefined;
    const i = resolved.sections.findIndex((s) => s.address === c.section);
    // past a collapsed section (every track hidden): it has no stop for the cursor
    let k = i + delta;
    while (i >= 0 && resolved.sections[k] && collapsedSection(resolved.sections[k])) k += delta;
    const to = resolved.sections[k];
    if (i < 0 || !to) return undefined;
    // onto its rows (or a grid's y) from the side it came in by: the first row shown going down, the last going up
    const which = delta > 0 ? 'first' : 'last';
    if (to.axis.system === c.system) return enter(plain(to.address, c.system, c.position), which);
    const part = projectOntoView(
      spec.model.value,
      resolved,
      rangesLocation(c.system, [c.position, c.position + 1]),
    ).parts.find((p) => p.system === to.axis.system);
    if (part?.ranges.length) return enter(plain(to.address, to.axis.system, part.ranges[0]), which);
    // nothing carries it: the first position the section shows
    const v = view.get(to.address);
    const p = v && !v.pending ? elementAt(to.axis, Math.ceil(v.start)) : undefined;
    return p === undefined ? undefined : enter(plain(to.address, to.axis.system, p), which);
  };

  /**
   * The rows shown from data row `a` to data row `b` in display order (a summary row's all its rows), as data rows:
   * sliced from the rows' axis in display order (`RowsLay.rowsAxis`, made once a display), so a million rows cost what
   * ten do.
   */
  const rowsBetween = (r: RowsOf, a: number, b: number): number[] => {
    const axis = r.lay.rowsAxis();
    const d = r.lay.display;
    const span = (y: number): [number, number] | undefined => {
      const o = cellOffset(axis, y);
      if (o === undefined) return undefined;
      // a summary row's rows are one run in the order, its first shown row first
      const m = d.members.get(r.lay.cursorRow(y));
      return [o, o + (m ? m.length : 1)];
    };
    const sa = span(a);
    const sb = span(b);
    if (!sa || !sb) return [];
    const lo = Math.min(sa[0], sb[0]);
    const hi = Math.max(sa[1], sb[1]);
    const out: number[] = [];
    const ranges = axis.ranges ?? [];
    let acc = 0;
    for (let k = 0; k + 1 < ranges.length && acc < hi; k += 2) {
      const w = ranges[k + 1] - ranges[k];
      const s = Math.max(lo, acc);
      const e = Math.min(hi, acc + w);
      if (e > s) out.push(ranges[k] + s - acc, ranges[k] + e - acc);
      acc += w;
    }
    return out;
  };
  /** The rect of cells from the cursor `from` to `to`: on rows the rows shown between them, on a grid `cellsRect`. */
  const blockRect = (from: Cursor, to: Cursor): Location | undefined => {
    if (from.ySystem === undefined || from.yPosition === undefined || to.yPosition === undefined) return undefined;
    const r = rowsOf(to.section);
    if (!r)
      return cellsRect(
        { ...from, ySystem: from.ySystem, yPosition: from.yPosition },
        { ...to, ySystem: from.ySystem, yPosition: to.yPosition },
      );
    const rows = rowsBetween(r, from.yPosition, to.yPosition);
    if (!rows.length) return undefined;
    const x0 = Math.min(from.position, to.position);
    const x1 = Math.max(from.position, to.position) + 1;
    return location({
      rects: [{ a: { system: to.system, ranges: [x0, x1] }, b: { system: from.ySystem, ranges: rows } }],
    });
  };

  /**
   * The block a Shift or Cmd/Ctrl+Shift arrow grows on a section of two axes: the cell it began on, what it keeps (the
   * others for Cmd/Ctrl+Shift, nothing for Shift) and the selection it made last. Its own, not the model's cell anchor:
   * on rows its rows are in display order, and Cmd/Ctrl+Shift keeps the others while Shift then grows it.
   */
  let block: { readonly from: Cursor; readonly base: Location; made: unknown } | undefined;
  const extendBlock = (c: Cursor, next: Cursor, add: boolean): boolean => {
    if (
      add ||
      !block ||
      block.made !== selection.value.value ||
      block.from.section !== c.section ||
      block.from.ySystem !== c.ySystem
    )
      block = { from: onCell(c), base: add ? selection.value.value.location : EmptyLocation, made: undefined };
    const rect = blockRect(block.from, next);
    if (!rect) return false;
    const loc = unionLocations(block.base, rect);
    selection.set({ elements: EmptyKeys, location: loc, shown: loc });
    block.made = selection.value.value;
    return true;
  };

  const placed = (ctx: CursorContext): true | string => (current(ctx) ? true : 'No section to move the cursor in');

  const move = (
    id: string,
    title: string,
    key: string,
    to: (c: Cursor) => Cursor | undefined,
    select: 'move' | 'extend' | 'add' = 'move',
    horizontal = false,
  ): CommandDefinition<CursorContext> => ({
    id,
    title,
    category: 'Cursor',
    shortcut: {
      key,
      ...(select !== 'move' ? { shift: true } : {}),
      ...(select === 'add' ? { mod: true } : {}),
    },
    ...(select === 'move' ? { scope: 'canvas' } : {}),
    when: placed,
    run: (_, ctx): CommandResult => {
      const c = current(ctx)!;
      // a cursor not yet shown is placed first: its first key shows it where it is
      const next = selection.cursor.value ? to(c) : c;
      if (!next) return { message: '' };
      const header = c.part === 'header' && c.ySystem !== undefined && c.yPosition !== undefined;
      // on the header the column it left is hidden: Shift with ← → moves nothing
      if (select !== 'move' && header && horizontal) return { message: '' };
      if (select !== 'move' && onRows(c) && !header) {
        // a cell of two axes: the block from where it began to the cursor's cell (Cmd/Ctrl+Shift: a new one)
        extendBlock(c, onCell(next), select === 'add');
        moveTo(next);
        return { message: select === 'add' ? 'Added' : 'Extended' };
      }
      if (select !== 'move') {
        // on the header a run of rows (in display order), else positions on the cursor's axis
        const system = header ? c.ySystem! : c.system;
        const p = header ? c.yPosition! : c.position;
        const q = header ? next.yPosition! : next.position;
        const axis = header ? view.rowsAxis(system) : sectionOf(next.section)?.axis;
        const cell = rangesLocation(system, [p, p + 1]);
        const anchor = { system, start: p, end: p + 1 };
        const a = selection.anchor;
        // a new range at the cursor, the others kept; else from the anchor (none on this system: the cursor's cell)
        if (select === 'add') selection.add(cell, anchor);
        else if (!a || a.system !== system) selection.set({ elements: EmptyKeys, location: cell, shown: cell }, anchor);
        selection.extend({ system, position: q }, axis);
      }
      moveTo(next);
      const verb = select === 'move' ? 'Moved' : select === 'add' ? 'Added' : 'Extended';
      return { message: `${verb}${crossed(c, next)}` };
    },
  });

  /** The rows a page key acts on with no cursor shown: the section's under the pointer, else the first band. */
  const pageRowsOf = (ctx: CursorContext): RowsOf | undefined => {
    if (ctx.section !== undefined) {
      const r = rowsOf(ctx.section);
      if (r) return r;
    }
    for (const s of layout()?.sections ?? []) {
      const r = rowsOf(s.section.address);
      if (r) return r;
    }
    return undefined;
  };
  /** A shown cursor's rows, else undefined. */
  const shownOnRows = () => {
    const c = selection.cursor.value;
    return c?.keyboard && onRows(c) ? c : undefined;
  };

  /**
   * Page Down and Page Up: a shown cursor on rows moves a page of rows and the band a page with it (`pagedTop`), the
   * row then shown; with none, the band scrolls a page (the rows' `rows.page-*`); at the band's end the key is the
   * page's.
   */
  const page = (id: string, title: string, dir: 1 | -1): CommandDefinition<CursorContext> => ({
    id,
    title,
    category: 'Cursor',
    scope: 'canvas',
    when: (ctx) => {
      const c = shownOnRows();
      if (c) return rowsOf(c.section) ? true : 'No rows to page';
      const r = pageRowsOf(ctx);
      if (!r) return 'No rows to page';
      return pagedTop(r.lay, topOf(r), dir) !== topOf(r) ? true : 'The rows are at their end';
    },
    run: (_, ctx): CommandResult => {
      const c = shownOnRows();
      const r = c ? rowsOf(c.section) : pageRowsOf(ctx);
      if (!r) return { message: '' };
      const top = pagedTop(r.lay, topOf(r), dir);
      view.setBand(r.track, { bandTop: top });
      if (!c) return { message: dir > 0 ? 'Rows a page down' : 'Rows a page up' };
      const next = rowStep(c, dir * pageRows(r.lay), true);
      if (next) moveTo(next, top);
      return { message: 'Moved' };
    },
  });

  /** Cmd/Ctrl+↑ ↓: a shown cursor on rows or a grid to its first or last row (with none the keys are the page's). */
  const rowsEnd = (id: string, title: string, end: -1 | 1): CommandDefinition<CursorContext> => ({
    id,
    title,
    category: 'Cursor',
    scope: 'canvas',
    when: () => (shownOnRows() ? true : 'No cursor on rows'),
    run: (): CommandResult => {
      const next = rowEnd(selection.cursor.value!, end);
      if (next) moveTo(next);
      return { message: next ? 'Moved' : '' };
    },
  });

  /** The clade Space on a row's header folds or unfolds: a summary row's, else the leaf's parent's (its innermost). */
  const cladeOf = (r: RowsOf, c: Cursor): string | undefined => {
    const l = r.lay;
    const cl = l.clades;
    if (!cl || !l.foldsApply || c.yPosition === undefined) return undefined;
    const i = rowIndex(r, c);
    const e = i >= 0 ? l.display.rowAt(i) : -1;
    if (e < 0 && i >= 0) return cl.name[-1 - e];
    const leaf = cl.rowLeaf[c.yPosition];
    if (leaf < 0) return undefined;
    let best = -1;
    for (let n = 0; n < cl.nodes; n++)
      if (cl.children[n].length && cl.lo[n] <= leaf && leaf < cl.hi[n] && (best < 0 || cl.depth[n] > cl.depth[best]))
        best = n;
    return best < 0 ? undefined : cl.name[best];
  };
  /** Folds or unfolds a clade of the cursor's rows, the row at the band's top kept there. */
  const toggleClade = (r: RowsOf, name: string): string => {
    const now = bandFolds(view, r.track, r.lay);
    const next = now.includes(name) ? now.filter((n) => n !== name) : [...now, name];
    setBandFolds(view, r.lay, r.track, next);
    return next.length > now.length ? 'Folded' : 'Unfolded';
  };
  /** Whether the cursor stands on a summary row (a folded clade's), and the row's layer. */
  const onSummary = (c: Cursor): RowsOf | undefined => {
    const r = rowsOf(c.section);
    if (!r || !onRows(c)) return undefined;
    const i = rowIndex(r, c);
    return i >= 0 && r.lay.display.rowAt(i) < 0 ? r : undefined;
  };

  /**
   * A row's extent on the columns: its first to its last letter where its cells are held or read (a lazy row's tile),
   * else every column (a lazy row not read yet).
   */
  const rowExtent = (l: RowsLay, row: number): Location | undefined => {
    const h = l.held;
    const codes = h ? h.codes.subarray(row * h.columns, (row + 1) * h.columns) : l.lazyCodes.get(row);
    if (!codes) return l.columns > 0 ? rangesLocation(l.system, [0, l.columns]) : undefined;
    let a = -1;
    let b = -1;
    for (let j = 0; j < codes.length; j++) {
      if (gap(codes[j])) continue;
      if (a < 0) a = j;
      b = j;
    }
    return a < 0 ? undefined : rangesLocation(l.system, [a, b + 1]);
  };

  const commands: CommandDefinition<CursorContext>[] = [
    move('cursor.prev', 'Cursor back, or onto the row’s header', keys.prev, back),
    move('cursor.next', 'Cursor on', keys.next, on),
    move('cursor.up', 'Cursor a row or a line up, or to the section above', keys.up, (c) => vertical(c, -1)),
    move('cursor.down', 'Cursor a row or a line down, or to the section below', keys.down, (c) => vertical(c, 1)),
    move(
      'cursor.extend-up',
      'Extend the selection a row or a line up',
      keys.up,
      (c) => (onRows(c) ? rowStep(c, -1, true) : line(c, 'up')),
      'extend',
    ),
    move(
      'cursor.extend-down',
      'Extend the selection a row or a line down',
      keys.down,
      (c) => (onRows(c) ? rowStep(c, 1, true) : line(c, 'down')),
      'extend',
    ),
    move('cursor.extend-prev', 'Extend the selection back', keys.prev, (c) => stepped(c, -1), 'extend', true),
    move('cursor.extend-next', 'Extend the selection on', keys.next, (c) => stepped(c, 1), 'extend', true),
    move('cursor.add-prev', 'Start a new range back from the cursor', keys.prev, (c) => stepped(c, -1), 'add', true),
    move('cursor.add-next', 'Start a new range on from the cursor', keys.next, (c) => stepped(c, 1), 'add', true),
    move(
      'cursor.add-up',
      'Start a new block or range a row or a line up from the cursor',
      keys.up,
      (c) => (onRows(c) ? rowStep(c, -1, true) : line(c, 'up')),
      'add',
    ),
    move(
      'cursor.add-down',
      'Start a new block or range a row or a line down from the cursor',
      keys.down,
      (c) => (onRows(c) ? rowStep(c, 1, true) : line(c, 'down')),
      'add',
    ),
    page('cursor.page-down', 'Cursor and rows a page down', 1),
    page('cursor.page-up', 'Cursor and rows a page up', -1),
    rowsEnd('cursor.first-row', 'Cursor to the first row', -1),
    rowsEnd('cursor.last-row', 'Cursor to the last row', 1),
    {
      id: 'cursor.toggle',
      title: 'Toggle the cursor’s cell in the selection, or its row’s clade',
      category: 'Cursor',
      shortcut: { key: ' ' },
      scope: 'canvas',
      when: (ctx) => (selection.cursor.value && current(ctx) ? true : 'No cursor'),
      run: (_, ctx): CommandResult => {
        const c = current(ctx)!;
        const r = rowsOf(c.section);
        if (c.part === 'header' && r && onRows(c)) {
          // the row's clade folds or unfolds (its summary row's, else its parent's); no tree: the row in the selection
          const name = cladeOf(r, c);
          if (name !== undefined) return { message: toggleClade(r, name) };
          selection.toggle({ system: c.ySystem!, position: c.yPosition! }, view.rowsAxis(c.ySystem!));
          return { message: 'Toggled' };
        }
        if (onRows(c) && c.ySystem !== undefined && c.yPosition !== undefined) {
          // a summary row's cell: its rows, added, else (all in already) taken out
          if (onSummary(c)) {
            const rect = blockRect(c, c);
            if (rect && !selection.add(rect)) selection.subtract(rect);
          } else selection.selectCell({ ...c, ySystem: c.ySystem, yPosition: c.yPosition }, 'toggle');
          return { message: 'Toggled' };
        }
        selection.toggle({ system: c.system, position: c.position }, sectionOf(c.section)?.axis);
        return { message: 'Toggled' };
      },
    },
    {
      id: 'cursor.focus',
      title: 'Focus the selection, else the cursor; open a folded clade',
      category: 'Cursor',
      shortcut: { key: 'Enter' },
      scope: 'canvas',
      when: () => (selection.cursor.value || selection.value.value.shown.parts.length ? true : 'Nothing to focus'),
      run: (): CommandResult => {
        const c = selection.cursor.value;
        // a summary row opens (its clade unfolds), the cursor kept on its first row
        const s = c && onSummary(c);
        if (c && s?.lay.clades) {
          const e = s.lay.display.rowAt(rowIndex(s, c));
          return { message: toggleClade(s, s.lay.clades.name[-1 - e]) === 'Unfolded' ? 'Opened' : '' };
        }
        // the selection, else the cursor's cell (on the header its row's extent), focused as an address is: the
        // leading sections move (a cursor on a follower focuses its leader), each to the hull of what it shows of it
        const shown = selection.value.value.shown;
        const r = c && rowsOf(c.section);
        const header = c?.part === 'header' && r && c.yPosition !== undefined;
        const loc = shown.parts.length
          ? shown
          : header
            ? rowExtent(r.lay, c.yPosition!)
            : c
              ? rangesLocation(c.system, [c.position, c.position + 1])
              : undefined;
        if (!loc) return { message: '' };
        const res = addresses.focusLocation(loc);
        const focused = addresses.movedToBreak ? HiddenMessage : 'Focused';
        return { message: res === 'moved' ? focused : res === 'pending' ? 'Kept until the view is placed' : '' };
      },
    },
  ];

  return { commands, dispose: stopClicks };
}
