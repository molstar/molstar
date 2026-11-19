/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `cursor` behavior (_spec/README.md "Behaviors"): the keyboard cursor (the interaction grammar): the step keys
// (`cursor.keys`) move it a position along its section's axis in display order, or a line up or down in a layout of
// lines (a text: the section's layout's step), else to the section above or below (its position carried along the
// view's link); Shift with a step extends the range the anchor is in (from the cursor's cell when there is no anchor),
// up and down by a line; Cmd/Ctrl+Shift with a step starts a new range at the cursor, the others kept; Space toggles
// the cursor's cell; Enter focuses the selection, else the cursor's cell. A cursor moved out of view brings its section
// to it by the layout's least move (`reveal.cursor`). Every move is a command.
//
// On a section of two axes the cursor has a row too (grids-and-rows, The keyboard cursor): on an alignment's rows ↑ ↓
// step display rows (the pinned band to the scrolled rows, a summary row one stop), Page Up and Page Down a page of
// rows, Cmd/Ctrl+↑ ↓ the first and the last row, ← from the first column onto the row's header and → back; Shift with
// an arrow extends the block (the rect from where it began); Space toggles the cell, on the header the row's clade;
// Enter opens a summary row. Moves are arithmetic on display rows (a row never loaded is stepped over by its index);
// the band scrolls to show the row by the least move. On a grid ↑ ↓ step its y and Cmd/Ctrl+↑ ↓ go to its ends.

import { EmptyElements, type Location, location, rangesLocation, subtractLocations } from '@molstar/seqx-core/query';
import type { RowsLay, ViewLayout } from '@molstar/seqx-core/render';
import { type Axis, projectOntoView } from '@molstar/seqx-core/view';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
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
import { bandFolds, pagedTop, pageRows, revealedTop, setBandFolds } from './rows-band.js';

/** The section a key acts on when there is no cursor yet (the one under the pointer, else the first). */
export interface CursorContext {
  readonly section?: string;
}

export interface CursorBehavior {
  readonly commands: readonly CommandDefinition<CursorContext>[];
}

/** A section's rows layer as the last frame laid it out, and its track's address (where its band's state is kept). */
interface RowsOf {
  readonly lay: RowsLay;
  readonly track: string;
}

/** The cursor less its row: a cell on the section's axis. */
const plain = (section: string, system: string, position: number): Cursor => ({ section, system, position });

/** The cursor on its row's cell (off the header). */
const onCell = (c: Cursor): Cursor => {
  const { part: _, ...rest } = c;
  return rest;
};

export function cursorBehavior(
  addresses: AddressModel,
  selection: SelectionModel,
  view: ViewModel,
  spec: SpecModel,
  config: ViewerConfig = DefaultViewerConfig,
  layout: () => ViewLayout | undefined = () => undefined,
): CursorBehavior {
  const keys = config.cursor.keys;
  const sectionOf = (address: string | undefined) => view.resolved.value?.sections.find((s) => s.address === address);

  /** A section's rows layer (one a section: the cursor has one row order), from the last frame. */
  const rowsOf = (address: string): RowsOf | undefined => {
    for (const s of layout()?.sections ?? []) {
      if (s.section.address !== address) continue;
      for (const band of s.lays)
        for (const l of band ?? [])
          if (l.type.name === 'rows') return { lay: l.lay as RowsLay, track: l.key.slice(0, l.key.lastIndexOf('/')) };
    }
    return undefined;
  };
  /** A grid section's y axis (its own system, in display order). */
  const gridY = (address: string): Axis | undefined => sectionOf(address)?.y?.axis;

  /** The band's top now: the model's (the last frame may not have drawn a move just made), else the layout's. */
  const topOf = (r: RowsOf) => view.band(r.track)?.bandTop ?? r.lay.top;

  /**
   * The display row a step left the cursor on (a pinned row is in the band and in its place: the one it stands on),
   * kept while it still shows the cursor's data row.
   */
  let lastRow = -1;
  /** The display row the cursor's data row is on: where a step left it, else its own (its summary row's), else the nearest. */
  const rowIndex = (r: RowsOf, c: Cursor): number => {
    const y = c.yPosition;
    const d = r.lay.display;
    if (y === undefined || c.ySystem !== r.lay.rowsSystem || !d.count) return -1;
    if (lastRow >= 0 && lastRow < d.count) {
      const e = d.rowAt(lastRow);
      if (e === y || (e < 0 && d.members.get(lastRow)?.includes(y))) return lastRow;
    }
    const i = d.displayRow(y);
    return i >= 0 ? i : d.nearestRow(y);
  };
  /** The cursor on display row `i` (a summary row by its first shown row), its part kept. */
  const atRow = (c: Cursor, r: RowsOf, i: number): Cursor | undefined => {
    const d = r.lay.display;
    if (i < 0 || i >= d.count) return undefined;
    let e = d.rowAt(i);
    if (e < 0) e = d.members.get(i)?.[0] ?? -1;
    if (e < 0) return undefined;
    lastRow = i;
    return { ...c, ySystem: r.lay.rowsSystem, yPosition: e };
  };
  /** The first row shown (the pinned band's first, else the band's top), or the last shown whole. */
  const shownRow = (r: RowsOf, which: 'first' | 'last'): number => {
    const l = r.lay;
    const d = l.display;
    const top = topOf(r);
    if (which === 'first') return d.pinned > 0 ? 0 : Math.min(d.count - 1, Math.ceil(top));
    const area = Math.max(0, l.box.height - d.pinned * l.rowScale) / l.rowScale;
    return Math.max(Math.min(d.pinned, d.count - 1), Math.min(d.count - 1, d.pinned + Math.floor(top + area) - 1));
  };
  /** The cursor placed on a section's rows or grid y (entering it), `which` row; unchanged on a section of one axis. */
  const enter = (c: Cursor, which: 'first' | 'last'): Cursor => {
    const r = rowsOf(c.section);
    if (r && r.lay.display.count) return atRow(c, r, shownRow(r, which)) ?? c;
    const y = gridY(c.section);
    const total = y && axisTotal(y);
    if (!y || !total) return c;
    const v = view.get(c.section)?.y;
    const o = which === 'first' ? Math.ceil(v?.start ?? 0) : Math.floor((v?.end ?? total) - 1);
    const p = elementAt(y, Math.max(0, Math.min(total - 1, o)));
    return p === undefined ? c : { ...c, ySystem: y.system, yPosition: p };
  };

  /** The cursor now, or one placed at the first position the section shows (its first row shown too). */
  const current = (ctx: CursorContext): Cursor | undefined => {
    const c = selection.cursor.value;
    if (c && sectionOf(c.section)) return c;
    const section = sectionOf(ctx.section) ?? view.resolved.value?.sections[0];
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
    const from = top ?? topOf(r);
    const next = revealedTop(r.lay, from, i);
    if (next !== topOf(r)) view.setBand(r.track, { bandTop: next });
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
   * (`clamp`), else undefined past them; undefined on a section of one axis or a cursor with no row there.
   */
  const rowStep = (c: Cursor, delta: number, clamp: boolean): Cursor | undefined => {
    const r = rowsOf(c.section);
    if (r && c.ySystem === r.lay.rowsSystem) {
      const n = r.lay.display.count;
      const i = rowIndex(r, c);
      if (i < 0) return undefined;
      const j = clamp ? Math.max(0, Math.min(n - 1, i + delta)) : i + delta;
      return atRow(c, r, j);
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

  /** The cursor shown there, its cell brought into view (its column unless on the header, its row on rows). */
  const moveTo = (c: Cursor, top?: number) => {
    selection.setCursor({ ...c, keyboard: true });
    if (c.part !== 'header') {
      const axis = sectionOf(c.section)?.axis;
      const o = axis && cellOffset(axis, c.position);
      if (o !== undefined) reveal(c.section, o);
    }
    revealRow(c, top);
  };

  /** The cursor on the section `delta` above or below, its position carried along the view's link. */
  const across = (c: Cursor, delta: number): Cursor | undefined => {
    const resolved = view.resolved.value;
    if (!resolved) return undefined;
    const i = resolved.sections.findIndex((s) => s.address === c.section);
    const to = resolved.sections[i + delta];
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
   * The block an extend on a section of two axes grows (Shift with an arrow): its first cell (its display row on
   * rows), and the selection it made last; a selection changed since starts a new block at the cursor.
   */
  let block: { readonly from: Cursor; readonly row: number; made: unknown } | undefined;
  /** The rect from the block's first cell to `c`: on rows the rows shown between them in display order, as data rows. */
  const blockRect = (from: Cursor, fromRow: number, c: Cursor): Location | undefined => {
    if (from.ySystem === undefined || from.yPosition === undefined || c.yPosition === undefined) return undefined;
    const r = rowsOf(c.section);
    if (!r)
      return cellsRect(
        { ...from, ySystem: from.ySystem, yPosition: from.yPosition },
        { ...c, ySystem: from.ySystem, yPosition: c.yPosition },
      );
    const d = r.lay.display;
    const i = rowIndex(r, c);
    if (i < 0 || fromRow < 0) return undefined;
    const lo = Math.min(i, fromRow);
    const hi = Math.max(i, fromRow);
    const rows: number[] = [];
    // data order shows display row i as data row i: one range, however many rows (arithmetic)
    if (d.identity) rows.push(lo, hi + 1);
    else
      for (let k = lo; k <= hi; k++) {
        const e = d.rowAt(k);
        if (e >= 0) rows.push(e, e + 1);
        else for (const m of d.members.get(k) ?? []) rows.push(m, m + 1);
      }
    const x0 = Math.min(from.position, c.position);
    const x1 = Math.max(from.position, c.position) + 1;
    return location({
      rects: [{ a: { system: c.system, ranges: [x0, x1] }, b: { system: from.ySystem, ranges: rows } }],
    });
  };
  /** Extends the block to `next` (from the cursor `c`'s cell when no block runs). */
  const extendBlock = (c: Cursor, next: Cursor): boolean => {
    const r = rowsOf(c.section);
    if (
      !block ||
      block.made !== selection.value.value ||
      block.from.section !== c.section ||
      block.from.ySystem !== c.ySystem
    )
      block = { from: onCell(c), row: r ? rowIndex(r, c) : -1, made: undefined };
    const rect = blockRect(block.from, block.row, next);
    if (!rect) return false;
    selection.set({ elements: EmptyElements, location: rect, shown: rect });
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
      if (select === 'extend' && onRows(c) && c.part !== 'header') {
        // a cell of two axes: the block from where the extend began to the cursor's cell
        extendBlock(c, onCell(next));
        moveTo(next);
        return { message: 'Extended' };
      }
      if (select !== 'move') {
        // on the header a run of rows (in display order), else positions on the cursor's axis
        const header = c.part === 'header' && c.ySystem !== undefined && c.yPosition !== undefined;
        const system = header ? c.ySystem! : c.system;
        const p = header ? c.yPosition! : c.position;
        const q = header ? next.yPosition! : next.position;
        const axis = header ? view.rowsAxis(system) : sectionOf(next.section)?.axis;
        const cell = rangesLocation(system, [p, p + 1]);
        const anchor = { system, start: p, end: p + 1 };
        const a = selection.anchor;
        // a new range at the cursor, the others kept; else from the anchor (none on this system: the cursor's cell)
        if (select === 'add') selection.add(cell, anchor);
        else if (!a || a.system !== system)
          selection.set({ elements: EmptyElements, location: cell, shown: cell }, anchor);
        selection.extend({ system, position: q }, axis);
      }
      moveTo(next);
      const verb = select === 'move' ? 'Moved' : select === 'add' ? 'Added' : 'Extended';
      return { message: `${verb}${crossed(c, next)}` };
    },
  });

  /** A rows command on a shown cursor's row (Page Up and Down, the first and the last row). */
  const rowsMove = (
    id: string,
    title: string,
    shortcut: { key: string; mod?: boolean },
    to: (c: Cursor) => Cursor | undefined,
    top?: (c: Cursor) => number | undefined,
    rowsOnly = false,
  ): CommandDefinition<CursorContext> => ({
    id,
    title,
    category: 'Cursor',
    shortcut,
    scope: 'canvas',
    // with no cursor shown the keys are the rows' own (Page Up and Down scroll the band)
    when: () => {
      const c = selection.cursor.value;
      if (!c?.keyboard || !onRows(c)) return 'No cursor on rows';
      return rowsOnly && !rowsOf(c.section) ? 'No cursor on rows' : true;
    },
    run: (): CommandResult => {
      const c = selection.cursor.value!;
      const next = to(c);
      if (!next) return { message: '' };
      moveTo(next, top?.(c));
      return { message: 'Moved' };
    },
  });

  /** Page Down and Page Up: the cursor a page of rows on, the band a page with it (`pagedTop`), the row shown. */
  const page = (dir: 1 | -1) => ({
    to: (c: Cursor) => {
      const r = rowsOf(c.section);
      return r ? rowStep(c, dir * pageRows(r.lay), true) : undefined;
    },
    top: (c: Cursor) => {
      const r = rowsOf(c.section);
      if (!r) return undefined;
      const top = pagedTop(r.lay, topOf(r), dir);
      view.setBand(r.track, { bandTop: top });
      return top;
    },
  });
  const pageDown = page(1);
  const pageUp = page(-1);

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
  /** The cursor's cell as a rect on a section of two axes (a summary row's: its rows), else undefined. */
  const cellRect = (c: Cursor): Location | undefined => {
    if (!onRows(c) || c.ySystem === undefined || c.yPosition === undefined) return undefined;
    const r = rowsOf(c.section);
    return r
      ? blockRect(c, rowIndex(r, c), c)
      : cellsRect(
          { ...c, ySystem: c.ySystem, yPosition: c.yPosition },
          { ...c, ySystem: c.ySystem, yPosition: c.yPosition },
        );
  };

  /** A held row's letters on the columns (its first to its last letter), undefined when not held or all gaps. */
  const rowLetters = (l: RowsLay, row: number): Location | undefined => {
    const h = l.held;
    if (!h || row < 0 || row >= h.rows) return undefined;
    const base = row * h.columns;
    let a = -1;
    let b = -1;
    for (let j = 0; j < h.columns; j++) {
      const k = h.codes[base + j];
      // a gap byte (`-`, `.`, none: the rows layer's)
      if (k === 45 || k === 46 || k === 0) continue;
      if (a < 0) a = j;
      b = j;
    }
    return a < 0 ? undefined : rangesLocation(l.system, [a, b + 1]);
  };

  const commands: CommandDefinition<CursorContext>[] = [
    move('cursor.prev', 'Cursor back, or onto the row’s header', keys.prev, back),
    move('cursor.next', 'Cursor on', keys.next, on),
    move('cursor.up', 'Cursor a row or a line up, or to the section above', keys.up, (c) =>
      onRows(c) ? (rowStep(c, -1, false) ?? across(c, -1)) : (line(c, 'up') ?? across(c, -1)),
    ),
    move('cursor.down', 'Cursor a row or a line down, or to the section below', keys.down, (c) =>
      onRows(c) ? (rowStep(c, 1, false) ?? across(c, 1)) : (line(c, 'down') ?? across(c, 1)),
    ),
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
    move('cursor.extend-prev', 'Extend the selection back', keys.prev, (c) => stepped(c, -1), 'extend'),
    move('cursor.extend-next', 'Extend the selection on', keys.next, (c) => stepped(c, 1), 'extend'),
    move('cursor.add-prev', 'Start a new range back from the cursor', keys.prev, (c) => stepped(c, -1), 'add'),
    move('cursor.add-next', 'Start a new range on from the cursor', keys.next, (c) => stepped(c, 1), 'add'),
    rowsMove('cursor.page-down', 'Cursor a page of rows down', { key: 'PageDown' }, pageDown.to, pageDown.top, true),
    rowsMove('cursor.page-up', 'Cursor a page of rows up', { key: 'PageUp' }, pageUp.to, pageUp.top, true),
    rowsMove('cursor.first-row', 'Cursor to the first row', { key: keys.up, mod: true }, (c) => rowEnd(c, -1)),
    rowsMove('cursor.last-row', 'Cursor to the last row', { key: keys.down, mod: true }, (c) => rowEnd(c, 1)),
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
        const rect = cellRect(c);
        if (rect) {
          // the cell's rect: out of the selection when all of it is in, else added
          const sel = selection.value.value.location;
          const inside =
            subtractLocations(location({ rects: rect.rects }), location({ rects: sel.rects })).rects.length === 0;
          if (inside) selection.subtract(rect);
          else selection.add(rect);
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
        const r = c && rowsOf(c.section);
        if (c && r && onRows(c)) {
          const i = rowIndex(r, c);
          const e = i >= 0 ? r.lay.display.rowAt(i) : 0;
          if (e < 0 && r.lay.clades)
            return { message: toggleClade(r, r.lay.clades.name[-1 - e]) === 'Unfolded' ? 'Opened' : '' };
        }
        // the selection, else the cursor's cell (on the header its row's letters), focused as an address is: the
        // leading sections move (a cursor on a follower focuses its leader), each to the hull of what it shows of it
        const shown = selection.value.value.shown;
        const header = c?.part === 'header' && r && c.yPosition !== undefined;
        const loc = shown.parts.length
          ? shown
          : header
            ? rowLetters(r.lay, c.yPosition!)
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

  return { commands };
}
