/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandDefinition } from '../commands/types.js';
import type { HistoryModel, UndoResult } from './history-model.js';

/**
 * Undo and Redo as commands (history §1, headless §0.1): the menu, ⌘Z and an agent do what the History view's buttons
 * do. An app registers them with its command manager, in its category. A refusal fails the command with its reason.
 */
export function historyCommands<C, Ctx = unknown>(
  history: HistoryModel<C>,
  category = 'Edit',
): CommandDefinition<Ctx>[] {
  const settle = async (p: Promise<UndoResult>, verb: string) => {
    const r = await p;
    if (!r.ok) throw new Error(r.reason);
    return { message: `${verb}: ${r.step.label}` };
  };
  return [
    {
      id: 'edit.undo',
      title: 'Undo',
      category,
      icon: 'undo-2',
      shortcut: { key: 'z', mod: true },
      when: () => history.undoBlocked.value ?? (history.canUndo.value ? true : 'Nothing to undo'),
      run: () => settle(history.undo(), 'Undone'),
    },
    {
      id: 'edit.redo',
      title: 'Redo',
      category,
      icon: 'redo-2',
      shortcut: { key: 'z', mod: true, shift: true },
      when: () => (history.canRedo.value ? true : 'Nothing to redo'),
      run: () => settle(history.redo(), 'Redone'),
    },
  ];
}
