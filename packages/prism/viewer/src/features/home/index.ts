/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createHomeActions } from '../../sources/home.js';
import type { ViewerExtension } from '../../spec.js';
import { homeShell } from './shell.js';

// The `home` feature (viewer-spec.md §3.1): Home, at the head of the rail: the registry's loading actions
// (`viewer.views.homeActions`), the features' import sections (the loaders'), files, the session, the records kept
// earlier, and the offer of a structure's annotations (with the `annotations` feature).

export const HomeFeature: ViewerExtension = {
  id: 'home',
  title: 'Home',
  builtIn: true,
  views: ({ scope, plugin, person }) => ({ homeActions: createHomeActions(scope, plugin, person) }),
  shell: (ctx) => homeShell(ctx),
};
