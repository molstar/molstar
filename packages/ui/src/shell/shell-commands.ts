/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandDefinition } from '../commands/types.js';
import { findView, mainNode } from './layout.js';
import type { ShellModel } from './shell-model.js';

/**
 * The shell's own commands (shell §3.9, headless §0.1): what its controls do, so the command menu, a key and an agent
 * do the same. An app registers them with its command manager (`manager.register(shellCommands(shell))` or its own
 * equivalent), in its category.
 */
export function shellCommands<Ctx = unknown>(
  shell: ShellModel,
  category = shell.labels.commandCategory,
): CommandDefinition<Ctx>[] {
  // The shell's words (shell §6): the app's, merged in the model.
  const t = shell.labels;
  // Several panes in the main area, or the active one in a window (a detached stage maximizes over the main area).
  const panes = () =>
    mainNode(shell.layout.value).root || findView(shell.layout.value, shell.activePane.value)?.window
      ? true
      : t.onePane;
  return [
    {
      id: 'view.expand',
      title: t.expandTitle,
      category,
      icon: 'maximize-2',
      alias: ['expand', 'focus mode'],
      run: () => {
        shell.toggleExpanded();
        return { message: shell.expanded.value ? t.expandedMessage : t.notExpandedMessage };
      },
    },
    {
      id: 'view.fullscreen',
      title: t.fullscreenTitle,
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
      title: t.toggleRailTitle,
      category,
      icon: 'chevrons-right',
      run: () => {
        shell.toggleRail();
        return { message: shell.railShown.value ? t.viewsShownMessage : t.viewsHiddenMessage };
      },
    },
    {
      id: 'view.maximize-pane',
      title: t.maximizePaneTitle,
      category,
      icon: 'maximize',
      alias: ['maximize', 'zoom pane'],
      // Only with several panes (shell §2.1), or the active pane in a window.
      when: panes,
      run: () => {
        // The port's rule: the active pane's region node, its ref in a window, else the content.
        shell.toggleMaximize(shell.maximizeTarget(shell.activePane.value));
        return { message: shell.maximized.value ? t.maximizedMessage : t.allPanesMessage };
      },
    },
    {
      id: 'view.grid-panes',
      title: t.gridTitle,
      category,
      icon: 'layout-grid',
      alias: ['grid', 'tile panes'],
      when: () => (mainNode(shell.layout.value).root ? true : t.onePane),
      run: () => {
        shell.arrangeGrid();
        return { message: t.gridMessage };
      },
    },
    {
      id: 'view.close-pane',
      title: t.closePaneTitle,
      category,
      icon: 'x',
      // The app's content is no pane to close; a view beside it is.
      when: () =>
        shell.activePane.value !== 'content' && shell.paneNode(shell.activePane.value) ? true : t.nothingToClose,
      run: () => {
        // The active pane falls back to the content as its ref leaves (shell §2.1).
        shell.close(shell.activePane.value);
        return { message: t.paneClosedMessage };
      },
    },
    {
      id: 'view.reset-layout',
      title: t.resetLayoutTitle,
      category,
      icon: 'list-restart',
      run: () => {
        shell.resetLayout();
        return { message: t.layoutResetMessage };
      },
    },
  ];
}
