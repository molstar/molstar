/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ApiResult, type OperationFn } from '../../api/area.js';
import { operationCommand } from '../../commands/operations.js';
import type { ViewerExtension } from '../../spec.js';
import {
  type ListParams,
  type OpenParams,
  type RemoteStateList,
  type RemoveParams,
  remoteStatesOperations,
  type UploadParams,
} from './operations.js';

declare module '../../api/session.js' {
  // What the feature adds to the core's `session` area (present while it is in the spec).
  interface SessionApi {
    readonly remoteList?: OperationFn<typeof ListParams, RemoteStateList>;
    readonly remoteOpen?: OperationFn<typeof OpenParams>;
    readonly remoteUpload?: OperationFn<typeof UploadParams>;
    readonly remoteRemove?: OperationFn<typeof RemoveParams>;
  }
}

// The `remote-states` feature (viewer-spec.md §3.1, parity.md §1): the state server's list, open, upload and remove as
// `session.remote.*`, and the Remote states card. The spec lists it when `ui.showRemoteStates` is on. Its views are a
// chunk of their own.

export const RemoteStatesFeature: ViewerExtension = {
  id: 'remote-states',
  title: 'Remote states',
  builtIn: true,
  contributes: () => remoteStatesOperations(),
  commands: (api) => {
    const { remoteList, remoteOpen, remoteUpload, remoteRemove } = api.session;
    if (!remoteList || !remoteOpen || !remoteUpload || !remoteRemove) return [];
    return [
      operationCommand(api, remoteList, {
        category: 'Remote states',
        icon: 'cloud',
        alias: ['state server', 'list states'],
      }),
      operationCommand(api, remoteOpen, {
        category: 'Remote states',
        icon: 'cloud',
        alias: ['open remote state', 'fetch state'],
      }),
      operationCommand(api, remoteUpload, {
        category: 'Remote states',
        icon: 'cloud-upload',
        alias: ['save to server', 'share state'],
        action: 'Upload',
      }),
      operationCommand(api, remoteRemove, {
        category: 'Remote states',
        icon: 'trash-2',
        alias: ['delete remote state'],
      }),
      {
        id: 'remote-states.open',
        title: 'Show the remote states',
        category: 'Remote states',
        icon: 'cloud',
        alias: ['remote states', 'state server'],
        // Opens the card: a view, not a change.
        async run() {
          const r = await api.view.open({ card: 'remote-states.overlay.states' });
          if (ApiResult.isFailure(r)) throw new Error(r.message);
          return {
            message: 'Opened the Remote states card.',
            ops: 'view.open({ card: "remote-states.overlay.states" })',
          };
        },
      },
    ];
  },
  ui: async (ctx) => (await import('./ui.js')).remoteStatesUI(ctx),
};

export type { RemoteStateEntry, RemoteStateList } from './operations.js';
