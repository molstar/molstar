/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { BehaviorProvider } from '@molstar/seqx-viewer/viewer-spec';
import { type HoverBehavior, hoverBehavior } from './hover.js';
import { type NavigateBehavior, navigateBehavior } from './navigate.js';

// The behaviors' providers (_spec/README.md "Composition"): what a viewer spec lists, by name.

export const HoverProvider: BehaviorProvider<undefined, HoverBehavior> = {
  name: 'hover',
  create: ({ hover, bindings }) => hoverBehavior(hover, bindings),
};

export const NavigateProvider: BehaviorProvider<undefined, NavigateBehavior> = {
  name: 'navigate',
  create: ({ view, bindings, config }) => navigateBehavior(view, bindings, config),
};
