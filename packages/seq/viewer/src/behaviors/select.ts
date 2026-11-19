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
// test select alike. A click with a position moves the keyboard cursor there.

import { type Location, rangesLocation } from '@molstar/seqx-core/query';
import type { Axis } from '@molstar/seqx-core/view';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import type { BindingTable, ModifierState } from '@molstar/seqx-viewer/bindings';
import {
  type AddressModel,
  type DragMode,
  type Hit,
  type HitTarget,
  hitRange,
  hitTarget,
  positionAt,
  type SelectionModel,
  type ViewModel,
} from '@molstar/seqx-viewer/models';

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
  /** A section's axis: spans are built in its display order (a range over a cyclic origin is two). */
  const axisOf = (section: string | undefined): Axis | undefined =>
    section === undefined ? undefined : view.resolved.value?.sections.find((s) => s.address === section)?.axis;

  /** A click with a position moves the cursor to it, on the section clicked. */
  const follow = (ctx: SelectContext) => {
    const h = ctx.hit;
    if (!h || ctx.section === undefined || h.position === undefined || Number.isNaN(h.position)) return;
    selection.setCursor({ section: ctx.section, system: h.system, position: Math.floor(h.position) });
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
      const changed = act(ctx.hit, axisOf(ctx.section));
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
    const r = ctx.range ?? (ctx.hit ? hitRange(ctx.hit, axisOf(ctx.section)) : undefined);
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

  let drag: { section: string; system: string } | undefined;

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
      if (!selection.beginDrag(hit, axisOf(section), DragSelectModes[action])) return false;
      drag = { section, system: hit.system };
      return true;
    },
    dragTo(section, x, hit, at) {
      const d = drag;
      if (!d || d.section !== section) return false;
      const v = view.get(section);
      const axis = view.resolved.value?.sections.find((s) => s.address === section)?.axis;
      if (!v || v.pending || !axis || v.width <= 0) return false;
      const offset = at ?? v.start + (Math.min(Math.max(x, 0), v.width - 1e-6) * (v.end - v.start)) / v.width;
      // the middle of the cell at the offset: an offset on a piece's end names no position past it
      const position = positionAt(axis, Math.floor(offset) + 0.5);
      if (position === undefined) return false;
      return selection.dragTo(position, hit && hit.system === d.system ? hit : undefined, axis);
    },
    endDrag() {
      drag = undefined;
      selection.endDrag();
    },
    cancelDrag() {
      drag = undefined;
      selection.cancelDrag();
    },
  };
}

const NeverAborted = new AbortController().signal;
