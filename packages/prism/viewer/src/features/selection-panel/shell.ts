/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { BuiltInShellContext } from '../../spec.js';
import { selectionSummaryModel } from '../../ui/interaction.js';
import { selectionSections, selectionTab } from './panels.js';

/** The Selection tab (the core's id, `selection`) and its own sections. */
export function selectionPanelShell({ shell, core, person, own }: BuiltInShellContext): () => void {
  const { sources } = core;
  const undo = [
    shell.registerTabs(own('tab', [selectionTab(selectionSummaryModel(sources.interaction, person))])),
    shell.registerSections(
      own(
        'section',
        selectionSections(person, sources.interaction.selection, sources.layers, (ref) => {
          // A component made of the selection, shown in the Structure panel where there is one.
          if (!('scene' in core.views)) return;
          shell.show('structure');
          // The tree follows the state on its own schedule: tried until the row is there (a second at most).
          let tries = 0;
          const reveal = () => {
            if (core.views.scene.reveal?.(ref) || ++tries > 20) return;
            setTimeout(reveal, 50);
          };
          reveal();
        }),
      ),
    ),
  ];
  return () => {
    for (const u of undo.reverse()) u();
  };
}
