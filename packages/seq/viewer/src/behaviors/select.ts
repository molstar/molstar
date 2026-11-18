/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `select` behavior (_spec/README.md "Behaviors"): clicks and drags to the selection model through the bindings
// (the interaction grammar: click sets, Shift+click extends from the anchor, Cmd/Ctrl+click toggles, a click on empty
// space or Esc clears, a drag on data selects a range). The granularity is the hit's: an item is an element, a cell a
// position. Every gesture is a command too, so a key, a menu and a test select alike. A click with a position moves the
// keyboard cursor there.

import type { Axis } from '@molstar/seqx-core/view';
import type { CommandDefinition, CommandResult } from '@molstar/ui/commands/types';
import type { BindingTable, ModifierState } from '@molstar/seqx-viewer/bindings';
import {
  type Hit,
  type HitTarget,
  hitTarget,
  positionAt,
  type SelectionModel,
  type ViewModel,
} from '@molstar/seqx-viewer/models';

/** The section a gesture is on, and what it hit (none: empty space). */
export interface SelectContext {
  readonly section?: string;
  readonly hit?: Hit;
}

export interface SelectBehavior {
  /** A click on `hit` (undefined: empty space) through the bindings; true when a binding took it. */
  click(section: string | undefined, hit: Hit | undefined, modifiers: ModifierState): boolean;
  /** A drag begins: true when the bindings make it a drag-select here and something selectable was hit. */
  beginDrag(section: string, hit: Hit | undefined, modifiers: ModifierState, target: HitTarget): boolean;
  /**
   * The drag moved to `x` (px from the section's data box's left edge; past it, clamped to the axis): the selection
   * is the range from where it began, snapped to whole items when it began on one and `hit` is an item.
   */
  dragTo(section: string, x: number, hit: Hit | undefined): boolean;
  endDrag(): void;
  readonly commands: readonly CommandDefinition<SelectContext>[];
}

export function selectBehavior(selection: SelectionModel, view: ViewModel, bindings: BindingTable): SelectBehavior {
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

  const commands: CommandDefinition<SelectContext>[] = [
    gesture('selection.set', 'Select', (h, a) => selection.select(h, a), 'Selected'),
    gesture('selection.extend', 'Extend the selection', (h, a) => selection.extend(h, a), 'Extended'),
    gesture('selection.toggle', 'Toggle in the selection', (h, a) => selection.toggle(h, a), 'Toggled'),
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
      if (bindings.resolve('drag', modifiers, target) !== 'drag-select' || !hit) return false;
      if (!selection.beginDrag(hit, axisOf(section))) return false;
      drag = { section, system: hit.system };
      return true;
    },
    dragTo(section, x, hit) {
      const d = drag;
      if (!d || d.section !== section) return false;
      const v = view.get(section);
      const axis = view.resolved.value?.sections.find((s) => s.address === section)?.axis;
      if (!v || v.pending || !axis || v.width <= 0) return false;
      const offset = v.start + (Math.min(Math.max(x, 0), v.width - 1e-6) * (v.end - v.start)) / v.width;
      // the middle of the cell at the offset: an offset on a piece's end names no position past it
      const position = positionAt(axis, Math.floor(offset) + 0.5);
      if (position === undefined) return false;
      return selection.dragTo(position, hit && hit.system === d.system ? hit : undefined, axis);
    },
    endDrag() {
      drag = undefined;
      selection.endDrag();
    },
  };
}

const NeverAborted = new AbortController().signal;
