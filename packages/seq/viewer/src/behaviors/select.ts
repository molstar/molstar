/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `select` behavior (_spec/README.md "Behaviors"): clicks and drags to the selection model through the bindings
// (the interaction grammar: click sets, Shift+click extends the anchor's range, Cmd/Ctrl+click toggles, a click on empty
// space or Esc clears; a drag on data selects a range, Cmd/Ctrl+drag adds it, Alt/Option+drag removes it, Shift+drag
// extends the anchor's range). The granularity is the hit's: an item is an element, a cell a position. Every gesture is
// a command too (`selection.add` and `selection.subtract` over a range or an address), so a key, a menu, a host and a
// test select alike. A click with a position moves the keyboard cursor there. On a grid's cell (two axes, not an
// alignment's rows) the same gestures select rects: a click a one-cell rect, Shift+click the rect from the anchor
// cell, Cmd/Ctrl+click toggles the cell, a drag a rect, Cmd/Ctrl+drag adds one, Alt+drag subtracts one.

import { EmptyKeys, type Location, location, rangesLocation } from '@molstar/seqx-core/query';
import { type Axis, blockPositionAt, blockSpan, blockUnitOf, type ResolvedSection } from '@molstar/seqx-core/view';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import type { BindingTable, ModifierState } from '@molstar/seqx-viewer/bindings';
import {
  type AddressModel,
  axisTotal,
  type DragMode,
  elementAt,
  type Emphasized,
  type GridCell,
  type Hit,
  type HitTarget,
  hitRange,
  hitTarget,
  positionAt,
  type SelectionModel,
  type ViewModel,
} from '@molstar/seqx-viewer/models';
import { type BlocksCell, blocksRect, type RectBlock } from './blocks-rect.js';

/** Positions `[start, end)` on a system: what `selection.add` and `selection.subtract` take besides an address. */
export interface SelectRange {
  readonly system: string;
  readonly start: number;
  readonly end: number;
}

/**
 * The section a gesture is on and what it hit (none: empty space); `selection.add` and `selection.subtract` take a
 * `range` or an `address` instead of a hit.
 */
export interface SelectContext {
  readonly section?: string;
  readonly hit?: Hit;
  readonly range?: SelectRange;
  readonly address?: string;
  /** `selection.copyLetters`: the stored letters of the selection's span, not what is shown. */
  readonly stored?: boolean;
}

/** The binding actions that drag-select, and what each does to the selection (the model's `DragMode`). */
export const DragSelectModes: Readonly<Record<string, DragMode>> = {
  'drag-select': 'replace',
  'drag-add': 'add',
  'drag-subtract': 'subtract',
  'drag-extend': 'extend',
};

/** Whether a binding's action is a drag-select of some mode. */
export function isDragSelect(action: string | undefined): boolean {
  return action !== undefined && Object.hasOwn(DragSelectModes, action);
}

export interface SelectBehavior {
  /** A click on `hit` (undefined: empty space) through the bindings; true when a binding took it. */
  click(section: string | undefined, hit: Hit | undefined, modifiers: ModifierState): boolean;
  /**
   * A drag begins: true when the bindings make it a drag-select here (of any mode: replace, add, subtract, extend)
   * and something selectable was hit.
   */
  beginDrag(section: string, hit: Hit | undefined, modifiers: ModifierState, target: HitTarget): boolean;
  /**
   * The drag moved to `x` (px from the section's data box's left edge; past it, clamped to the axis), or to the display
   * unit `offset` the placement has under the pointer (a text's in reading order, the nearest past a line's end): the
   * range from where it began, snapped to whole items when it began on one and `hit` is an item, applied by the drag's
   * mode.
   */
  dragTo(section: string, x: number, hit: Hit | undefined, offset?: number): boolean;
  endDrag(): void;
  /** Cancels a drag-select: the selection before it comes back. */
  cancelDrag(): void;
  readonly commands: readonly CommandDefinition<SelectContext>[];
}

/**
 * `addresses`: what `selection.add` and `selection.subtract` resolve an `address` with (positions, elements by their
 * extents); without it they take ranges and hits only.
 */
export function selectBehavior(
  selection: SelectionModel,
  view: ViewModel,
  bindings: BindingTable,
  addresses?: AddressModel,
): SelectBehavior {
  /**
   * The axis a hit's positions are on and its first offset in the section: the section's (spans are built in its
   * display order: a range over a cyclic origin is two), or on a section of blocks on several systems (its axis its own
   * units, on no system) the hit block's own.
   */
  const sectionOf = (section: string | undefined) =>
    section === undefined ? undefined : view.resolved.value?.sections.find((x) => x.address === section);
  const placeOf = (
    section: string | undefined,
    hit?: Hit,
  ): { axis: Axis; offset: number; block: boolean } | undefined => {
    const s = section === undefined ? undefined : view.resolved.value?.sections.find((x) => x.address === section);
    if (!s) return undefined;
    const b = hit?.block !== undefined && s.axis.system !== hit.system ? s.blocks?.[hit.block] : undefined;
    return b ? { axis: b.axis, offset: b.offset ?? 0, block: true } : { axis: s.axis, offset: 0, block: false };
  };
  // a rows hit's positions are rows, along the rows layer's display order (Shift+click extends in it)
  const axisOf = (section: string | undefined, hit?: Hit): Axis | undefined =>
    (hit ? view.rowsAxis(hit.system) : undefined) ?? placeOf(section, hit)?.axis;

  /**
   * A grid's cell: a hit with a position on two axes, its y no alignment's rows, and no extent (a matrix's cell is its
   * item, or bin, selected as a ranges click selects: _spec/README.md "Behaviors").
   */
  const cellOf = (hit: Hit | undefined): GridCell | undefined =>
    hit &&
    !hit.element &&
    !hit.extent &&
    hit.ySystem !== undefined &&
    hit.yPosition !== undefined &&
    hit.position !== undefined &&
    !view.rowsAxis(hit.ySystem)
      ? { system: hit.system, position: hit.position, ySystem: hit.ySystem, yPosition: hit.yPosition }
      : undefined;

  /**
   * A click with a position moves the cursor to it, on the section clicked: a cell of two axes with its row (a rows
   * cell, a grid's; a summary row's cell its first row), a row's label or a summary row's to its header (the column the
   * cursor was on kept, else the first shown).
   */
  const follow = (ctx: SelectContext) => {
    const h = ctx.hit;
    const section = ctx.section;
    if (!h || section === undefined) return;
    // a summary row's hit names its clade: the cursor stands on one of its rows (the row shows them all)
    const summary = h.target === 'summary-row' ? h.extent?.parts[0] : undefined;
    const ySystem = summary ? summary.system : h.ySystem;
    const yPosition = summary ? summary.ranges[0] : h.yPosition;
    const y = ySystem !== undefined && yPosition !== undefined ? Math.floor(yPosition) : undefined;
    const header = h.target === 'row' || (summary !== undefined && h.position === undefined);
    if (header && y !== undefined) {
      const s = view.resolved.value?.sections.find((x) => x.address === section);
      const now = selection.cursor.value;
      const v = view.get(section);
      const first = s && v && !v.pending ? elementAt(s.axis, Math.ceil(v.start)) : undefined;
      const position = now?.section === section && now.system === s?.axis.system ? now.position : first;
      if (!s || position === undefined) return;
      selection.setCursor({ section, system: s.axis.system, position, ySystem, yPosition: y, part: 'header' });
      return;
    }
    if (h.position === undefined || Number.isNaN(h.position)) return;
    selection.setCursor({
      section,
      system: h.system,
      position: Math.floor(h.position),
      ...(y !== undefined ? { ySystem, yPosition: y } : {}),
    });
  };

  const gesture = (
    id: string,
    title: string,
    act: (hit: Hit | undefined, axis: Axis | undefined) => boolean,
    done: string,
  ): CommandDefinition<SelectContext> => ({
    id,
    title,
    category: 'Selection',
    when: (ctx) => (ctx.hit ? true : 'Nothing to select here'),
    run: (_, ctx): CommandResult => {
      const cell = cellOf(ctx.hit);
      const how = id === 'selection.extend' ? 'extend' : id === 'selection.toggle' ? 'toggle' : 'replace';
      const changed = cell ? selection.selectCell(cell, how) : act(ctx.hit, axisOf(ctx.section, ctx.hit));
      follow(ctx);
      return { message: changed ? done : '' };
    },
  });

  /** What `selection.add` or `selection.subtract` acts on: the address's positions, the range, else the hit's. */
  const positionsOf = (ctx: SelectContext): { loc: Location; anchor?: SelectRange } | undefined => {
    if (ctx.address !== undefined) {
      const loc = addresses?.locate(ctx.address).location;
      return loc && (loc.parts.length || loc.rects.length) ? { loc } : undefined;
    }
    const r = ctx.range ?? (ctx.hit ? hitRange(ctx.hit, axisOf(ctx.section, ctx.hit)) : undefined);
    if (!r || !(r.end > r.start)) return undefined;
    return { loc: rangesLocation(r.system, [r.start, r.end]), anchor: r };
  };

  const positions = (
    id: string,
    title: string,
    act: (what: { loc: Location; anchor?: SelectRange }) => boolean,
    done: string,
  ): CommandDefinition<SelectContext> => ({
    id,
    title,
    category: 'Selection',
    when: (ctx) =>
      ctx.address !== undefined && !addresses
        ? 'No addresses to resolve'
        : positionsOf(ctx)
          ? true
          : 'No range, address or position given',
    run: (_, ctx): CommandResult => {
      const what = positionsOf(ctx);
      return { message: what && act(what) ? done : '' };
    },
  });

  const commands: CommandDefinition<SelectContext>[] = [
    gesture('selection.set', 'Select', (h, a) => selection.select(h, a), 'Selected'),
    gesture('selection.extend', 'Extend the selection', (h, a) => selection.extend(h, a), 'Extended'),
    gesture('selection.toggle', 'Toggle in the selection', (h, a) => selection.toggle(h, a), 'Toggled'),
    positions('selection.add', 'Add to the selection', (w) => selection.add(w.loc, w.anchor), 'Added'),
    positions('selection.subtract', 'Remove from the selection', (w) => selection.subtract(w.loc), 'Removed'),
    {
      id: 'selection.clear',
      title: 'Clear the selection',
      category: 'Selection',
      shortcut: { key: 'Escape' },
      scope: 'canvas',
      when: () =>
        selection.value.value.elements.keys.length || selection.value.value.shown.parts.length
          ? true
          : 'Nothing selected',
      run: (): CommandResult => {
        selection.clear();
        return { message: 'Cleared' };
      },
    },
  ];

  /**
   * The drag: its section and system, and the axis its positions are read on (a block's: its offset in the section);
   * on a section of blocks its block, the anchor's unit, the mode and the selection as it began (a drag across blocks).
   */
  let drag:
    | {
        section: string;
        system: string;
        axis: Axis;
        offset: number;
        block: boolean;
        across?: { index: number; unit: number; mode: DragMode; before: Emphasized };
      }
    | undefined;
  /** A rect drag on a grid: its section and its cell's two systems. */
  let rect: { section: string; system: string; ySystem: string } | undefined;
  /**
   * A rect drag on a grid over blocks (grids-and-rows' grid over a complex): the blocks, the first cell, the mode and the
   * selection as it began; each move sets the rect between the two cells, split at the blocks' boundaries.
   */
  let over:
    | {
        section: string;
        blocks: readonly RectBlock[];
        from: BlocksCell;
        mode: string;
        before: Emphasized;
      }
    | undefined;
  const blocksCell = (hit: Hit): BlocksCell | undefined =>
    hit.block !== undefined && hit.yBlock !== undefined && hit.position !== undefined && hit.yPosition !== undefined
      ? { block: hit.block, position: hit.position, yBlock: hit.yBlock, yPosition: hit.yPosition }
      : undefined;

  return {
    commands,
    click(section, hit, modifiers) {
      const action = bindings.resolve('click', modifiers, hitTarget(hit));
      const c = action === undefined ? undefined : commands.find((x) => x.id === action);
      if (!c) return false;
      const ctx: SelectContext = { ...(section !== undefined ? { section } : {}), ...(hit ? { hit } : {}) };
      if (c.when?.(ctx) === true || !c.when) c.run({}, ctx, NeverAborted);
      return true;
    },
    beginDrag(section, hit, modifiers, target) {
      const action = bindings.resolve('drag', modifiers, target);
      if (action === undefined || !isDragSelect(action) || !hit) return false;
      const cell = cellOf(hit);
      if (cell) {
        // a grid over a section of blocks: the rect is split at the chains' boundaries (_spec/README.md "Behaviors",
        // A rect over blocks)
        const s = view.resolved.value?.sections.find((x) => x.address === section);
        const from = s?.layout === 'grid' && s.blocks ? blocksCell(hit) : undefined;
        if (s?.blocks && from) {
          const blocks = s.blocks.map((b) => ({ system: b.axis.system, axis: b.axis, offset: b.offset ?? 0 }));
          over = { section, blocks, from, mode: DragSelectModes[action], before: selection.value.value };
          return true;
        }
        selection.beginRectDrag(cell, DragSelectModes[action]);
        rect = { section, system: cell.system, ySystem: cell.ySystem };
        return true;
      }
      const at = placeOf(section, hit);
      const before = selection.value.value;
      if (!at || !selection.beginDrag(hit, at.axis, DragSelectModes[action])) return false;
      drag = { section, system: hit.system, axis: at.axis, offset: at.offset, block: at.block };
      // on a section of blocks the anchor's place: a drag past its block selects the shown positions between
      const blocks = at.block ? sectionOf(section)?.blocks : undefined;
      const unit =
        blocks && hit.block !== undefined && hit.position !== undefined
          ? blockUnitOf(blocks[hit.block], Math.floor(hit.position))
          : undefined;
      if (blocks && unit !== undefined)
        drag.across = { index: hit.block!, unit, mode: DragSelectModes[action], before };
      return true;
    },
    dragTo(section, x, hit, at) {
      if (over) {
        const to = hit && over.section === section && cellOf(hit) ? blocksCell(hit) : undefined;
        const r = to && blocksRect(over.blocks, over.from, to);
        if (!r) return false;
        const was = selection.value.value;
        // from the selection as the drag began: replaced, added to, or taken out of
        if (over.mode === 'add' || over.mode === 'subtract') {
          selection.set(over.before);
          if (over.mode === 'add') selection.add(r);
          else selection.subtract(r);
        } else selection.set({ elements: EmptyKeys, location: r, shown: r });
        return selection.value.value !== was;
      }
      if (rect) {
        // past the plot the cell stays the last one under the pointer
        const cell = cellOf(hit);
        if (rect.section !== section || !cell || cell.system !== rect.system || cell.ySystem !== rect.ySystem)
          return false;
        return selection.rectDragTo(cell.position, cell.yPosition);
      }
      const d = drag;
      if (!d || d.section !== section) return false;
      const v = view.get(section);
      const axis = d.axis;
      const total = axisTotal(axis);
      if (!v || v.pending || total === undefined || v.width <= 0) return false;
      const unit = at ?? v.start + (Math.min(Math.max(x, 0), v.width - 1e-6) * (v.end - v.start)) / v.width;
      // past the anchor's block: the shown positions of each block between, a part a system (genome coordinates)
      const s = d.across ? sectionOf(section) : undefined;
      if (d.across && s?.blocks) {
        const loc = acrossBlocks(s, d.across.index, d.across.unit, unit);
        if (loc) {
          const was = selection.value.value;
          if (d.across.mode === 'add' || d.across.mode === 'subtract') {
            selection.set(d.across.before);
            if (d.across.mode === 'add') selection.add(loc);
            else selection.subtract(loc);
          } else selection.set({ elements: EmptyKeys, location: loc, shown: loc });
          return selection.value.value !== was;
        }
      }
      // the placement's offset (a ring's by its angle), else the x's; a block's offsets from its first (past its ends:
      // its end cells, never another block's positions)
      const offset = unit - d.offset;
      let cell = d.block ? Math.min(Math.max(Math.floor(offset), 0), total - 1) : Math.floor(offset);
      // a window over a cyclic origin: its units past the total are the next turn's
      if (cell >= total && view.wraps(section)) cell -= total;
      // the middle of the cell at the offset: an offset on a piece's end names no position past it
      const position = positionAt(axis, cell + 0.5);
      if (position === undefined) return false;
      // round a ring, or across a window's origin: the shorter arc between the ends (`view.arcs`)
      return selection.dragTo(position, hit && hit.system === d.system ? hit : undefined, axis, view.arcs(section));
    },
    endDrag() {
      drag = undefined;
      rect = undefined;
      over = undefined;
      selection.endDrag();
      selection.endRectDrag();
    },
    cancelDrag() {
      drag = undefined;
      if (over) selection.set(over.before);
      over = undefined;
      if (rect) selection.endRectDrag(true);
      rect = undefined;
      selection.cancelDrag();
    },
  };
}

const NeverAborted = new AbortController().signal;

/** The block holding a section unit (the last starting at or before it). */
function blockIndexOf(blocks: NonNullable<ResolvedSection['blocks']>, u: number): number {
  let k = 0;
  while (k + 1 < blocks.length && (blocks[k + 1].offset ?? 0) <= u) k++;
  return k;
}

/**
 * A drag across blocks (genome coordinates, "Selection in display order"): from the anchor's cell (block `from`, its
 * unit) to the cell at `unit`, the shown positions of each block between in display order: the anchor's block from its
 * cell to its edge, the blocks between whole, the last from its edge to the pointer's cell, a part a system (ranges
 * ascending on it). Undefined while the pointer is in the anchor's block (the block's own drag).
 */
export function acrossBlocks(
  section: ResolvedSection,
  from: number,
  anchor: number,
  unit: number,
): Location | undefined {
  const blocks = section.blocks;
  if (!blocks?.length) return undefined;
  const to = blockIndexOf(blocks, unit);
  if (to === from) return undefined;
  // the two cells' units: each end's whole cell (a scaled block's cell is its scale wide)
  const cellOf = (k: number, u: number) => {
    const b = blocks[k];
    const p = blockPositionAt(b, Math.min(Math.max(u, b.offset ?? 0), (b.offset ?? 0) + blockSpan(b) - 1e-9), 'start');
    return p === undefined ? undefined : blockUnitOf(b, Math.floor(p));
  };
  const a = cellOf(from, anchor);
  const c = cellOf(to, unit);
  if (a === undefined || c === undefined) return undefined;
  const scaleOf = (k: number) => blocks[k].scale ?? 1;
  const lo = Math.min(a, c);
  const hi = to > from ? c + scaleOf(to) : a + scaleOf(from);
  const byStart = to > from ? a : c;
  const parts = new Map<string, number[]>();
  for (let k = Math.min(from, to); k <= Math.max(from, to); k++) {
    const b = blocks[k];
    const o = b.offset ?? 0;
    const s0 = Math.max(k === Math.min(from, to) ? byStart : o, o, lo);
    const s1 = Math.min(o + blockSpan(b), hi);
    if (!(s1 > s0)) continue;
    const p0 = blockPositionAt(b, s0, 'start');
    const p1 = blockPositionAt(b, s1, 'end');
    if (p0 === undefined || p1 === undefined) continue;
    const r = parts.get(b.axis.system) ?? [];
    r.push(Math.floor(Math.min(p0, p1) + 1e-6), Math.ceil(Math.max(p0, p1) - 1e-6));
    parts.set(b.axis.system, r);
  }
  return location({ parts: [...parts].map(([system, ranges]) => ({ system, ranges })) });
}
