/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ViewerExtension } from '../../spec.js';
import { MoveApi } from './operations.js';

declare module '../../api.js' {
  // The `move` area (structure-move.md): handles in 3D to move and turn a structure. Optional: not in the default spec.
  interface ExtensionApis {
    move: MoveApi;
  }
}

// The `structure-move` feature (viewer-spec.md §3.1, structure-move.md): the Selection tab's Move section, handles in
// 3D on a structure (the gizmo library's transform handles), each drag one `structure.transform`. Opt-in: an app adds
// it to its spec (the Prism app does); its views are a chunk of their own.

export const MoveSectionId = 'structure-move.section.move';

export const StructureMoveFeature: ViewerExtension = {
  id: 'structure-move',
  title: 'Move structures',
  requires: ['structure'],
  api: { move: (ctx) => new MoveApi(ctx) },
  ui: async (ctx) => (await import('./ui.js')).structureMoveUI(ctx),
};
