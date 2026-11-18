/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

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
  create: ({ view, bindings, config }) => navigateBehavior(view, bindings, config),
};

export const SelectProvider: BehaviorProvider<undefined, SelectBehavior> = {
  name: 'select',
  create: ({ selection, view, bindings }) => selectBehavior(selection, view, bindings),
};

export const CursorProvider: BehaviorProvider<undefined, CursorBehavior> = {
  name: 'cursor',
  create: ({ selection, view, spec, config }) => cursorBehavior(selection, view, spec, config),
};

export const EditingProvider: BehaviorProvider<undefined, EditingBehavior> = {
  name: 'editing',
  create: ({ spec }) => editingBehavior(spec),
};
