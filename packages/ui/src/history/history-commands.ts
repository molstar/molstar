/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandDefinition, KeyChord } from '../commands/types.js';
import type { HistoryModel, UndoResult } from './history-model.js';

/** The commands' ids: a view reads their keys from the app's command manager by these (history §4). */
export const HistoryCommandIds = { undo: 'edit.undo', redo: 'edit.redo' } as const;

/** Undo's and Redo's own keys (⌘Z and ⇧⌘Z on macOS, Ctrl elsewhere): before the app's or the person's keymap. */
export const HistoryShortcuts: { readonly undo: KeyChord; readonly redo: KeyChord } = {
  undo: { key: 'z', mod: true },
  redo: { key: 'z', mod: true, shift: true },
};

/** The commands' words (history §6); the refusals are the model's (`HistoryModel.labels`). */
export interface HistoryCommandLabels {
  readonly undo: string;
  readonly redo: string;
  /** What a done Undo says ("Undone: Set color"). */
  readonly undone: (step: string) => string;
  readonly redone: (step: string) => string;
}

export const DefaultHistoryCommandLabels: HistoryCommandLabels = {
  undo: 'Undo',
  redo: 'Redo',
  undone: (step) => `Undone: ${step}`,
  redone: (step) => `Redone: ${step}`,
};

export interface HistoryCommandsOptions {
  /** The commands' category in the app's menu ('Edit'). */
  readonly category?: string;
  readonly labels?: Partial<HistoryCommandLabels>;
}

/**
 * Undo and Redo as commands (history §1, headless §0.1): the menu, ⌘Z and an agent do what the History view's buttons
 * do. An app registers them with its command manager, in its category. A refusal fails the command with its reason;
 * while a run is in flight both are greyed with the model's `running` reason.
 */
export function historyCommands<C, Ctx = unknown>(
  history: HistoryModel<C>,
  options: HistoryCommandsOptions = {},
): CommandDefinition<Ctx>[] {
  const category = options.category ?? 'Edit';
  const l = { ...DefaultHistoryCommandLabels, ...options.labels };
  const settle = async (p: Promise<UndoResult>, said: (step: string) => string) => {
    const r = await p;
    if (!r.ok) throw new Error(r.reason);
    return { message: said(r.step.label) };
  };
  return [
    {
      id: HistoryCommandIds.undo,
      title: l.undo,
      category,
      icon: 'undo-2',
      shortcut: HistoryShortcuts.undo,
      when: () => history.undoBlocked.value?.reason ?? true,
      run: () => settle(history.undo(), l.undone),
    },
    {
      id: HistoryCommandIds.redo,
      title: l.redo,
      category,
      icon: 'redo-2',
      shortcut: HistoryShortcuts.redo,
      when: () =>
        history.running.value > 0
          ? history.labels.running
          : history.canRedo.value
            ? true
            : history.labels.nothingToRedo,
      run: () => settle(history.redo(), l.redone),
    },
  ];
}
