/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ViewerExtension } from '../../spec.js';
import { MoveApi } from './operations.js';

declare module '../../api.js' {
  // The `move` area (move.md): handles in 3D to move and turn a model (a structure, a particle list, a shape, a volume).
  // Optional: not in the default spec.
  interface ExtensionApis {
    move: MoveApi;
  }
}

// The `move` feature (viewer-spec.md §3.1, move.md): the Selection tab's Move section, handles in 3D on a model (the
// gizmo library's transform handles), each drag one transform operation of its kind (`structure.transform`,
// `particles.transform`, `shape.transform`, `volume.transform`). Opt-in: an app adds it to its spec (the Prism app
// does); its views are a chunk of their own.

export const MoveSectionId = 'move.section.move';

export const MoveFeature: ViewerExtension = {
  id: 'move',
  title: 'Move models',
  requires: ['structure'],
  api: { move: (ctx) => new MoveApi(ctx) },
  ui: async (ctx) => (await import('./ui.js')).moveUI(ctx),
};
