/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSceneSource } from '../../sources/scene.js';
import { createSceneLookModel } from '../../sources/look.js';
import { createPresetsModel, createQuickStylesModel } from '../../sources/structure-tab.js';
import type { ViewerExtension } from '../../spec.js';
import { structurePanelShell } from './shell.js';

// The `structure-panel` feature (viewer-spec.md §3.1): the Structure panel: the scene tree and its inspectors (a
// structure's contents, a component's representations, a representation's quick settings and form, a map's levels and
// transfer function), Preset and Styles, the look's scene defaults, and their commands. Its models (`viewer.views.scene`, …) are made headless
// too: the agent's account reads the scene. On the scene's API, the `structure` feature.

export const StructurePanelFeature: ViewerExtension = {
  id: 'structure-panel',
  title: 'Structure panel',
  builtIn: true,
  requires: ['structure'],
  views: ({ scope, plugin, person, openCommand, sequence }) => ({
    scene: createSceneSource(scope, plugin, person, openCommand, sequence),
    presets: createPresetsModel(scope, plugin, person),
    look: createSceneLookModel(scope, plugin, person),
    quickStyles: createQuickStylesModel(person),
  }),
  shell: (ctx) => structurePanelShell(ctx),
};
