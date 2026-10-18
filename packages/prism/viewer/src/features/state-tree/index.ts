/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ApiResult } from '../../api/area.js';
import { operationCommand } from '../../commands/operations.js';
import type { ViewerExtension } from '../../spec.js';
import { StateApi } from './operations.js';

declare module '../../api.js' {
  // The feature's area (present while it is in the spec, as by default).
  interface ExtensionApis {
    state: StateApi;
  }
}

// The `state-tree` feature (viewer-spec.md §3.1): the State tree tab, and the `state` area over the plugin's
// data tree: `get`, `remove`, `visibility`, `update`, `action`. Its views are a chunk of their own.

export const StateTreeFeature: ViewerExtension = {
  id: 'state-tree',
  title: 'State tree',
  builtIn: true,
  api: { state: (ctx) => new StateApi(ctx) },
  commands: (api) => {
    const state = api.state;
    if (!state) return [];
    return [
      operationCommand(api, state.visibility, {
        category: 'State',
        icon: 'eye',
        alias: ['show cell', 'hide cell', 'state visibility'],
      }),
      operationCommand(api, state.remove, {
        category: 'State',
        icon: 'trash-2',
        alias: ['delete cell', 'remove cell', 'state remove'],
      }),
      {
        id: 'state-tree.open',
        title: 'Show the state tree',
        category: 'State',
        icon: 'list-tree',
        alias: ['state tree', 'state', 'cells', 'transforms'],
        // Opens the tab: a view, not a change.
        async run() {
          const r = await api.view.open({ tab: 'state-tree.tab.tree' });
          if (ApiResult.isFailure(r)) throw new Error(r.message);
          return { message: 'Opened the State tree tab.', ops: 'view.open({ tab: "state-tree.tab.tree" })' };
        },
      },
    ];
  },
  ui: async (ctx) => (await import('./ui.js')).stateTreeUI(ctx),
};

export type { StateCellDetails, StateCellInfo } from './operations.js';
