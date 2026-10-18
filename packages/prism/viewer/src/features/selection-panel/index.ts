/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ViewerExtension } from '../../spec.js';
import { selectionPanelShell } from './shell.js';

// The `selection-panel` feature (viewer-spec.md §3.1): the Selection tab: what is selected, and its slot's sections
// (style the selection, make a component of it; Measure and Superpose with their features). Selecting, focusing and
// the selection bar are the core's; the tab's actions are the scene's API, the `structure` feature.

export const SelectionPanelFeature: ViewerExtension = {
  id: 'selection-panel',
  title: 'Selection panel',
  builtIn: true,
  requires: ['structure'],
  shell: (ctx) => selectionPanelShell(ctx),
};
