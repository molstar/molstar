/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandDefinition } from '@molstar/prism-ui';

// The ways into the Inspector, without the Inspector: its page id and the commands that open it. The page itself is
// `@molstar/prism-assistant/inspector` (`inspectorPage`), so the Assistant's tab and card do not pull it in.

/** The Inspector's page (`ShellModel.openPage`), opened at a focus: `trace`, `tools`, `limits`, `turn:<id>`, `tool:<name>`. */
export const InspectorPageId = 'inspector';

/** The Inspector's views, as the ways in name them. */
export type InspectorFocus = 'trace' | 'tools' | 'limits';

/**
 * The commands that open the Inspector at a view (Agent trace on ⌥T, Agent tools, Agent limits), for an app's command
 * registry; `open` opens the page (`shell.openPage(InspectorPageId, focus)`).
 */
export function inspectorCommands(open: (focus: InspectorFocus) => void): CommandDefinition<any>[] {
  const command = (
    focus: InspectorFocus,
    title: string,
    alias: readonly string[],
    extra: Partial<CommandDefinition<any>> = {},
  ): CommandDefinition<any> => ({
    id: `assistant.${focus}`,
    title,
    category: 'Assistant',
    alias,
    ...extra,
    run: () => {
      open(focus);
      return { message: '' };
    },
  });
  return [
    command('trace', 'Agent trace', ['trace', 'inspector', 'agent', 'steps', 'what the agent did'], {
      icon: 'list-tree',
      shortcut: { key: 't', alt: true },
    }),
    command('tools', 'Agent tools', ['tools', 'skills', 'inspector'], { icon: 'wrench' }),
    command('limits', 'Agent limits', ['limits', 'budget', 'context window', 'thinking'], { icon: 'gauge' }),
  ];
}
