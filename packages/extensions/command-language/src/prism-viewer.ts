/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ViewerExtension } from '@molstar/prism-viewer/spec';
import { ApiResult } from '@molstar/prism-viewer/api/area';
import { VmdLinesPresetEntry } from './plugin/vmd-preset.js';
import { commandLanguageCommands, commandPrefixes } from './prism-viewer/commands.js';
import { ScriptApi } from './prism-viewer/operations.js';
import { ConsoleTab } from './prism-viewer/parts.js';
import { CommandLanguageSkill } from './prism-viewer/skill.js';

// PyMOL and VMD commands in the Prism viewer (_spec README): an opt-in extension, as Mol*'s others. It brings the
// `script` area (`script.run`: a submission as one step), the command menu's `pymol `, `vmd ` and `!` prefixes, and
// VMD's default look as a representation preset. Eager: the menu's rows parse each keystroke, which a lazy chunk could
// not do before it loaded (a prefix's rows are made synchronously).

export const CommandLanguageExtension: ViewerExtension = {
  id: 'command-language',
  title: 'PyMOL and VMD commands',
  plugin: { registry: [VmdLinesPresetEntry] },
  api: { script: (ctx) => new ScriptApi(ctx) },
  commands: (api) => [
    ...commandLanguageCommands(api),
    {
      id: 'script.console',
      title: 'Open the command console (PyMOL, VMD)',
      category: 'Commands',
      icon: 'terminal',
      alias: ['console', 'pymol console', 'vmd console', 'tk console'],
      // Opens the console's tab: a view, not a change.
      async run() {
        // The console's view is a chunk of its own: right after the viewer starts its tab may not be there yet.
        let r = await api.view.open({ tab: ConsoleTab });
        for (let i = 0; i < 20 && ApiResult.isFailure(r); i++) {
          await new Promise((done) => setTimeout(done, 100));
          r = await api.view.open({ tab: ConsoleTab });
        }
        if (ApiResult.isFailure(r)) throw new Error(r.message);
        return { message: 'Opened the console.', ops: `view.open({ tab: "${ConsoleTab}" })` };
      },
    },
  ],
  commandPrefixes,
  skills: [CommandLanguageSkill],
  // The console's view: a chunk of its own, loaded with the views.
  ui: async (ctx) => (await import('./prism-viewer/console.js')).consoleUI(ctx),
};
