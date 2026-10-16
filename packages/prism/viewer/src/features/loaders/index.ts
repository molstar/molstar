/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { OperationFn } from '../../api/area.js';
import type { Opened, SessionOpenParams } from '../../api/session.js';
import { operationCommand } from '../../commands/operations.js';
import type { ViewerExtension } from '../../spec.js';
import {
  type DataParams,
  type EntryParams,
  type Loaded,
  loadersOperations,
  type MvsParams,
  type SvsParams,
  type TrajectoryParams,
  type UrlParams,
  type VolumeParams,
} from './operations.js';

declare module '../../api/load.js' {
  // What the feature adds to the core's `load` area (present while it is in the spec, as by default).
  interface LoadApi {
    readonly entry: OperationFn<typeof EntryParams, Loaded>;
    readonly url: OperationFn<typeof UrlParams, Loaded>;
    readonly data: OperationFn<typeof DataParams, Loaded>;
    readonly mvs: OperationFn<typeof MvsParams, Loaded>;
    readonly svs: OperationFn<typeof SvsParams, Loaded>;
    /** `session.open` under the loaders' name (the app's `snapshot-url`, `snapshot-id`). */
    readonly snapshot: OperationFn<typeof SessionOpenParams, Opened>;
    readonly trajectory: OperationFn<typeof TrajectoryParams, Loaded>;
    readonly volume: OperationFn<typeof VolumeParams, Loaded>;
  }
}

// The `loaders` feature (viewer-spec.md §3.1, plan 3.2): the database and URL loaders as operations of the core's
// `load` area (`entry`, `url`, `data`, `mvs`, `svs`, `snapshot`, `trajectory`, `volume`), their commands, and the URL's load
// parameters (`runUrlLoads`, which the app calls once the viewer is ready; app.md §4). It has no views: Home's own
// "Download Structure" (every source with its options) is the form for an id, and a second one beside it did the same;
// `load.entry` stays for the command menu, tools and scripts. `load.svs` shows a SeqViewSpec in the viewer's sequence
// (sequence.md §6).

export const LoadersFeature: ViewerExtension = {
  id: 'loaders',
  title: 'Loaders',
  builtIn: true,
  contributes: () => loadersOperations(),
  commands: (api) => {
    const load = api.load;
    if (
      !load.entry ||
      !load.url ||
      !load.data ||
      !load.mvs ||
      !load.svs ||
      !load.snapshot ||
      !load.trajectory ||
      !load.volume
    ) {
      return [];
    }
    return [
      operationCommand(api, load.entry, {
        category: 'Load',
        icon: 'download',
        action: 'Download',
        alias: ['pdb', 'download', 'fetch', 'afdb', 'alphafold', 'emdb', 'model archive', 'by id'],
      }),
      operationCommand(api, load.url, {
        category: 'Load',
        icon: 'link',
        action: 'Load',
        alias: ['open url', 'import url', 'download url'],
      }),
      operationCommand(api, load.data, {
        category: 'Load',
        icon: 'clipboard',
        action: 'Load',
        alias: ['paste', 'raw data', 'text', 'pdb text', 'from clipboard'],
      }),
      operationCommand(api, load.mvs, {
        category: 'Load',
        icon: 'folder-open',
        action: 'Load',
        alias: ['molviewspec', 'story', 'mvsj', 'mvsx'],
      }),
      operationCommand(api, load.svs, {
        category: 'Load',
        icon: 'folder-open',
        action: 'Load',
        alias: ['seqviewspec', 'svs', 'svsj', 'sequence document', 'sequence view'],
      }),
      operationCommand(api, load.snapshot, {
        category: 'Load',
        icon: 'folder-open',
        action: 'Load',
        alias: ['session', 'state', 'molj', 'molx'],
      }),
      operationCommand(api, load.trajectory, {
        category: 'Load',
        icon: 'trajectory',
        action: 'Load',
        alias: ['md', 'molecular dynamics', 'xtc', 'dcd', 'trr', 'coordinates', 'topology'],
      }),
      operationCommand(api, load.volume, {
        category: 'Load',
        icon: 'volume',
        action: 'Load',
        alias: ['map', 'density', 'mrc', 'ccp4', 'cube', 'dsn6', 'em map'],
      }),
    ];
  },
};

export {
  type UrlLoad,
  type UrlLoadResult,
  LoadUrlOptions,
  LoadUrlParameters,
  loadsFromUrl,
  parseSearch,
  runUrlLoads,
} from './url.js';
