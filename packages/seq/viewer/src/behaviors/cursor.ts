/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `cursor` behavior (_spec/README.md "Behaviors", Cursor): the keyboard cursor's commands. Its moves are in
// `cursor-moves.ts` (internal); the provider gives them the last frame's layout, which its rows and a grid's y read.

import type { CommandDefinition } from '@molstar/ui/commands/types';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import type { AddressModel, SelectionModel, SpecModel, ViewModel } from '@molstar/seqx-viewer/models';
import { createCursor } from './cursor-moves.js';

/** The section a key acts on when there is no cursor yet (the one under the pointer, else the first). */
export interface CursorContext {
  readonly section?: string;
}

export interface CursorBehavior {
  readonly commands: readonly CommandDefinition<CursorContext>[];
}

/** The cursor's commands over the models alone (no frame: a section's rows and a grid's y are not reached). */
export function cursorBehavior(
  addresses: AddressModel,
  selection: SelectionModel,
  view: ViewModel,
  spec: SpecModel,
  config: ViewerConfig = DefaultViewerConfig,
): CursorBehavior {
  return createCursor(addresses, selection, view, spec, config, () => undefined);
}
