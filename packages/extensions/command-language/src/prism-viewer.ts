/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ViewerExtension } from '@molstar/prism-viewer/spec';
import { VmdLinesPresetEntry } from './plugin/vmd-preset.js';
import { commandLanguageCommands, commandPrefixes } from './prism-viewer/commands.js';
import { ScriptApi } from './prism-viewer/operations.js';

// PyMOL and VMD commands in the Prism viewer (_spec README): an opt-in extension, as Mol*'s others. It brings the
// `script` area (`script.run`: a submission as one step), the command menu's `pymol `, `vmd ` and `!` prefixes, and
// VMD's default look as a representation preset. Eager: the menu's rows parse each keystroke, which a lazy chunk could
// not do before it loaded (a prefix's rows are made synchronously).

export const CommandLanguageExtension: ViewerExtension = {
  id: 'command-language',
  title: 'PyMOL and VMD commands',
  plugin: { registry: [VmdLinesPresetEntry] },
  api: { script: (ctx) => new ScriptApi(ctx) },
  commands: (api) => commandLanguageCommands(api),
  commandPrefixes,
};
