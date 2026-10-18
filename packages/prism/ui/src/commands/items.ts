/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Menu entries by command id (commands.md §7): an entry takes its label, icon, key and whether it can run from the
// command, and runs through the manager (its form opens where it has params), so a menu does not repeat an action's
// logic and a removed command takes its entries with it. DOM-free; `CommandButton` is the reactive button.

import type { IconName } from '../components/icon/icon.js';
import { formatChord, isMac } from './keys.js';
import type { CommandManager } from './manager.js';
import type { CommandValues } from './types.js';

/** A menu entry made from a command: fits the context menu's and the canvas toolbar's menu items. */
export interface CommandItem {
  id: string;
  label: string;
  icon?: IconName;
  /** The command's key now (after the keymap), as text: `⌘K`. */
  shortcut?: string;
  /** A short explanation (the caller's). */
  hint?: string;
  /** Why it cannot run now (its `when`). */
  disabled?: string;
  run(): void;
}

export interface CommandItemOptions {
  /** In place of the command's title (a menu's shorter word). */
  label?: string;
  hint?: string;
  /** Values it opens with (the row it is for: `{ component: ref }`); a command with params opens its form at them. */
  values?: CommandValues;
  /** The entry's id; the command's by default. */
  id?: string;
  /** Told when the command refused to open (a run in flight, `when`). */
  onBlocked?(reason: string): void;
  mac?: boolean;
}

/**
 * A menu entry for a command, as it is now (a menu is built when it opens); none when no such command is registered
 * (removed by an override, or its feature is not there), so callers filter: `[...].filter((x) => !!x)`.
 */
export function commandItem(
  manager: CommandManager<any>,
  id: string,
  options: CommandItemOptions = {},
): CommandItem | undefined {
  const command = manager.get(id);
  if (!command) return undefined;
  const chord = manager.shortcut(command);
  const blocked = manager.blocked(command);
  return {
    id: options.id ?? id,
    label: options.label ?? command.title,
    ...(command.icon ? { icon: command.icon } : {}),
    ...(chord ? { shortcut: formatChord(chord, options.mac ?? isMac()) } : {}),
    ...(options.hint ? { hint: options.hint } : {}),
    ...(blocked ? { disabled: blocked } : {}),
    run: () => {
      const r = manager.open(id, undefined, options.values);
      if (r.kind === 'blocked') options.onBlocked?.(r.reason);
    },
  };
}

/** Entries for several commands, those that are registered. */
export function commandItems(
  manager: CommandManager<any>,
  entries: readonly (string | ({ command: string } & CommandItemOptions))[],
): CommandItem[] {
  return entries
    .map((e) => (typeof e === 'string' ? commandItem(manager, e) : commandItem(manager, e.command, e)))
    .filter((x): x is CommandItem => !!x);
}
