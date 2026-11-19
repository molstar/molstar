/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SelectActions } from '@molstar/seqx-viewer/bindings';
import { cellOffset } from '@molstar/seqx-viewer/models';
import type { BehaviorProvider } from '@molstar/seqx-viewer/viewer-spec';
import { type AddressingBehavior, addressingBehavior } from './addressing.js';
import { type CopyWrite, copyLettersCommand } from './copy.js';
import type { CursorBehavior } from './cursor.js';
import { createCursor } from './cursor-moves.js';
import { type EditingBehavior, editingBehavior } from './editing.js';
import { type HoverBehavior, hoverBehavior } from './hover.js';
import { type NavigateBehavior, navigateBehavior } from './navigate.js';
import { type NavigatorBehavior, navigatorBehavior } from './navigator.js';
import { type RowsBehavior, RowsActions, rowsBehavior } from './rows.js';
import { RowsRailGesture } from './rows-rail.js';
import { type SelectBehavior, selectBehavior } from './select.js';
import { type ThemingBehavior, themingBehavior } from './theming.js';
import { type WindowBehavior, windowBehavior } from './window.js';

// The behaviors' providers (_spec/README.md "Composition"): what a viewer spec lists, by name.

export const HoverProvider: BehaviorProvider<undefined, HoverBehavior> = {
  name: 'hover',
  create: ({ hover, bindings }) => hoverBehavior(hover, bindings),
};

export const AddressingProvider: BehaviorProvider<undefined, AddressingBehavior> = {
  name: 'addressing',
  create: ({ addresses }) => addressingBehavior(addresses),
};

export const NavigateProvider: BehaviorProvider<undefined, NavigateBehavior> = {
  name: 'navigate',
  // a step zoom with no pointer is about the keyboard cursor on its section
  create: ({ view, bindings, config, selection }) =>
    navigateBehavior(view, bindings, config, (section) => {
      const c = selection.cursor.value;
      const axis = c && view.resolved.value?.sections.find((s) => s.address === section)?.axis;
      if (!c || !axis || axis.system !== c.system || !view.groupOf(section).includes(c.section)) return undefined;
      const o = cellOffset(axis, c.position);
      return o === undefined ? undefined : o + 0.5;
    }),
};

/** The select behavior's params: where `selection.copyLetters` writes (the platform's clipboard by default). */
export interface SelectParams {
  readonly write?: CopyWrite;
}

export const SelectProvider: BehaviorProvider<SelectParams | undefined, SelectBehavior> = {
  name: 'select',
  // the selection's actions only with the behavior that runs them: else its bindings are reported, never dead
  actions: SelectActions,
  // the selection's gestures, and copying its letters (_spec/README.md "Behaviors", Copy)
  create: (ctx, params) => {
    const b = selectBehavior(ctx.selection, ctx.view, ctx.bindings, ctx.addresses);
    return { ...b, commands: [...b.commands, copyLettersCommand(ctx, params?.write)] };
  },
};

export const CursorProvider: BehaviorProvider<undefined, CursorBehavior> = {
  name: 'cursor',
  create: ({ addresses, selection, view, spec, config, layout }) =>
    createCursor(addresses, selection, view, spec, config, layout),
};

export const EditingProvider: BehaviorProvider<undefined, EditingBehavior> = {
  name: 'editing',
  create: ({ spec, view, selection }) => editingBehavior(spec, { view, selection }),
};

/** A window's drags and its slider (a ring's, a navigator's): a host may leave it out (_spec/README.md "Behaviors"). */
export const WindowProvider: BehaviorProvider<undefined, WindowBehavior> = {
  name: 'window',
  create: ({ view, config }) => windowBehavior(view, config),
};

/**
 * A navigator's gestures on its own window behavior (_spec/README.md "Behaviors", Navigator): a host may leave it out,
 * and its navigators and rails are then drawn with no targets.
 */
export const NavigatorProvider: BehaviorProvider<undefined, NavigatorBehavior> = {
  name: 'navigator',
  create: ({ view, config }) => navigatorBehavior(view, config),
};

export const LegendProvider: BehaviorProvider<undefined, ThemingBehavior> = {
  name: 'legend',
  create: ({ spec, addresses, legends }) => themingBehavior(spec, addresses, legends),
};

/** An alignment's rows bands: folds, pins, hidden rows, order, the rows' zoom and scroll (_spec/README.md "Behaviors"). */
export const RowsProvider: BehaviorProvider<undefined, RowsHandle> = {
  name: 'rows',
  actions: RowsActions,
  create: ({ spec, view, layout, config }) =>
    Object.assign(rowsBehavior(spec, view, layout), {
      railGesture: () => new RowsRailGesture(view, config),
    }),
};

/** The rows behavior's handle: the rows bands' policy, and the rows rail's drags (a viewer with no rows has none). */
export interface RowsHandle extends RowsBehavior {
  /** A new gesture on a rows rail (the input keeps one). */
  railGesture(): RowsRailGesture;
}
