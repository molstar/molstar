/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { cellOffset } from '@molstar/seqx-viewer/models';
import type { BehaviorProvider } from '@molstar/seqx-viewer/viewer-spec';
import { type AddressingBehavior, addressingBehavior } from './addressing.js';
import { type CursorBehavior, cursorBehavior } from './cursor.js';
import { type EditingBehavior, editingBehavior } from './editing.js';
import { type HoverBehavior, hoverBehavior } from './hover.js';
import { type NavigateBehavior, navigateBehavior } from './navigate.js';
import { type SelectBehavior, selectBehavior } from './select.js';

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

export const SelectProvider: BehaviorProvider<undefined, SelectBehavior> = {
  name: 'select',
  create: ({ selection, view, bindings, addresses }) => selectBehavior(selection, view, bindings, addresses),
};

export const CursorProvider: BehaviorProvider<undefined, CursorBehavior> = {
  name: 'cursor',
  create: ({ addresses, selection, view, spec, config }) => cursorBehavior(addresses, selection, view, spec, config),
};

export const EditingProvider: BehaviorProvider<undefined, EditingBehavior> = {
  name: 'editing',
  create: ({ spec, view, selection }) => editingBehavior(spec, { view, selection }),
};
