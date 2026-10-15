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
import { SnapshotApi } from './operations.js';

declare module '../../api.js' {
  // The `snapshot` area (operations.md §2): the plugin's snapshot manager as operations.
  interface ExtensionApis {
    snapshot: SnapshotApi;
  }
}

// The `snapshots` feature (viewer-spec.md §3.1): the `snapshot` area, the Snapshots tab and the playback card. The views
// are a chunk of their own.

export const SnapshotsFeature: ViewerExtension = {
  id: 'snapshots',
  title: 'Snapshots',
  builtIn: true,
  api: { snapshot: (ctx) => new SnapshotApi(ctx) },
  commands: (api) => {
    const s = api.snapshot;
    if (!s) return [];
    return [
      operationCommand(api, s.add, {
        category: 'Snapshots',
        icon: 'plus',
        alias: ['snapshot', 'save view', 'save state', 'bookmark'],
        action: 'Add',
      }),
      operationCommand(api, s.apply, { category: 'Snapshots', icon: 'images', alias: ['restore', 'go to snapshot'] }),
      operationCommand(api, s.next, {
        category: 'Snapshots',
        icon: 'chevron-right',
        alias: ['next snapshot'],
      }),
      operationCommand(api, s.previous, {
        category: 'Snapshots',
        icon: 'chevron-left',
        alias: ['previous snapshot', 'back'],
      }),
      operationCommand(api, s.play, { category: 'Snapshots', icon: 'play', alias: ['slideshow', 'tour', 'story'] }),
      operationCommand(api, s.stop, { category: 'Snapshots', icon: 'pause', alias: ['pause', 'stop slideshow'] }),
      {
        id: 'snapshots.open',
        title: 'Show the snapshots',
        category: 'Snapshots',
        icon: 'images',
        alias: ['snapshots', 'states'],
        // Opens the tab: a view, not a change.
        async run() {
          const r = await api.view.open({ tab: 'snapshots.tab.snapshots' });
          if (ApiResult.isFailure(r)) throw new Error(r.message);
          return { message: 'Opened the Snapshots tab.', ops: 'view.open({ tab: "snapshots" })' };
        },
      },
    ];
  },
  ui: async (ctx) => (await import('./ui.js')).snapshotsUI(ctx),
};
