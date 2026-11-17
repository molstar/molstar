/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandDefinition } from '../commands/types.js';
import { mainNode } from './layout.js';
import type { ShellModel } from './shell-model.js';

/**
 * The shell's own commands (shell §3.9, headless §0.1): what its controls do, so the command menu, a key and an agent
 * do the same. An app registers them with its command manager (`manager.register(shellCommands(shell))` or its own
 * equivalent), in its category.
 */
export function shellCommands<Ctx = unknown>(shell: ShellModel, category = 'View'): CommandDefinition<Ctx>[] {
  return [
    {
      id: 'view.expand',
      title: 'Expand the main area',
      category,
      icon: 'maximize-2',
      alias: ['expand', 'focus mode'],
      run: () => {
        shell.toggleExpanded();
        return { message: shell.expanded.value ? 'Expanded: Escape to come back' : 'Back from expanded' };
      },
    },
    {
      id: 'view.fullscreen',
      title: 'Fullscreen',
      category,
      icon: 'fullscreen',
      alias: ['full screen'],
      when: () => shell.fullscreenAvailable.value,
      run: () => {
        shell.requestFullscreen();
        return {};
      },
    },
    {
      id: 'view.toggle-rail',
      title: 'Show or hide the views',
      category,
      icon: 'chevrons-right',
      run: () => {
        shell.toggleRail();
        return { message: shell.railShown.value ? 'Views shown' : 'Views hidden' };
      },
    },
    {
      id: 'view.maximize-pane',
      title: 'Maximize the pane',
      category,
      icon: 'maximize',
      alias: ['maximize', 'zoom pane'],
      // Only with several panes in the main area (shell §2.1).
      when: () => (mainNode(shell.layout.value).root ? true : 'There is one pane'),
      run: () => {
        const node = shell.paneNode(shell.activePane.value) ?? 'content';
        shell.toggleMaximize(node);
        return { message: shell.maximized.value ? 'Maximized: Escape to come back' : 'All panes' };
      },
    },
    {
      id: 'view.grid-panes',
      title: 'Arrange the panes as a grid',
      category,
      icon: 'layout-grid',
      alias: ['grid', 'tile panes'],
      when: () => (mainNode(shell.layout.value).root ? true : 'There is one pane'),
      run: () => {
        shell.arrangeGrid();
        return { message: 'Panes arranged as a grid' };
      },
    },
    {
      id: 'view.close-pane',
      title: 'Close the pane',
      category,
      icon: 'x',
      // The app's content is no pane to close; a view beside it is.
      when: () =>
        shell.activePane.value !== 'content' && shell.paneNode(shell.activePane.value) ? true : 'Nothing to close',
      run: () => {
        shell.close(shell.activePane.value);
        shell.setActivePane('content');
        return { message: 'Pane closed' };
      },
    },
    {
      id: 'view.reset-layout',
      title: 'Reset the layout',
      category,
      icon: 'list-restart',
      run: () => {
        shell.resetLayout();
        return { message: 'Layout reset' };
      },
    },
  ];
}
