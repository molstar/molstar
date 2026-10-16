/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView, HelpPanel } from '@molstar/prism-ui';
import type { BuiltInShellContext } from '../../spec.js';
import { CanvasBindings } from '../../ui/bindings.js';
import { HelpProse } from '../../views/help/help-prose.js';

/** Help (the core's id, `help`) and its command. */
export function helpShell({ shell, core, own }: BuiltInShellContext): () => void {
  const undo = [
    shell.registerTabs(
      own('tab', [
        {
          id: 'help',
          title: 'Help',
          icon: 'circle-help',
          order: 2.8,
          place: 'foot',
          render: domView(() => (
            <HelpPanel commands={core.managers.commands} bindings={CanvasBindings} intro={() => <HelpProse />} />
          )),
        },
      ]),
    ),
    core.registerCommands([
      {
        id: 'help.show',
        title: 'Keyboard and mouse',
        alias: ['help', 'shortcuts', 'keys', 'bindings', 'controls'],
        category: 'Help',
        icon: 'circle-help',
        shortcut: { key: '?' },
        run: () => {
          shell.show('help');
          return { message: '' };
        },
      },
    ]),
  ];
  return () => {
    for (const u of undo.reverse()) u();
  };
}
