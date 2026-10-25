/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView } from '@molstar/prism-ui';
import type { BuiltInShellContext } from '../../spec.js';
import { StructureCommandGroups, structureCommands } from '../../ui/structure-commands.js';
import { StructureTab } from '../../views/structure/structure-tab.js';

/** The Structure tab (the core's id, `structure`), and Preset's and Styles' commands in their groups. */
export function structurePanelShell({ shell, core, own, openCommandMenu }: BuiltInShellContext): () => void {
  const { views } = core;
  const undo = [
    shell.registerTabs(
      own('tab', [
        {
          id: 'structure',
          compactBar: true,
          title: 'Structure',
          icon: 'structure',
          order: 1,
          render: domView(({ ctx }) => (
            <StructureTab
              shell={ctx.shell}
              scene={views.scene}
              look={views.look}
              // Preset and Styles open the command menu in their groups (`structure-commands.ts`).
              openCommands={(group) => openCommandMenu(group)}
            />
          )),
        },
      ]),
    ),
    core.registerCommands(structureCommands(views.presets, views.quickStyles)),
    core.managers.commands.registerGroups(StructureCommandGroups),
  ];
  return () => {
    for (const u of undo.reverse()) u();
  };
}
