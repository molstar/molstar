/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Changing keys (commands §6.6): what Help lists, its conflicts, and a remap recorded, cleared or reset through the
// keymap, with what each did said once. DOM-free: the view (`help-view.tsx`) renders it, a test or an agent drives it.

import { derive, type ReadonlyState, State } from '@molstar/core/reactive';
import { chordId, formatChord } from './keys.js';
import type { KeyConflict } from './keymap.js';
import type { CommandManager } from './manager.js';
import type { CommandDefinition, KeyChord } from './types.js';

/** A mouse or touch binding as Help lists it: what it does, and the input (keys and a gesture, in order). */
export interface InputBinding {
  readonly action: string;
  /** ['Shift', 'Click'], ['Right drag']. */
  readonly input: readonly string[];
  /** Another way to do the same (a trackpad gesture, a touch). */
  readonly or?: readonly string[];
}

/** Bindings that belong together ("Camera", "Selection"): an app (or a mode of it) gives its own. */
export interface InputBindingGroup {
  readonly title: string;
  readonly bindings: readonly InputBinding[];
}

/** What changing keys says (commands §6.6): strings, or functions for what names a command or a key. */
export interface KeysLabels {
  /** A key recorded: the command, its new keys as shown, the commands it now shares them with. */
  readonly recorded: (command: string, keys: string, sharedWith: readonly string[]) => string;
  readonly cleared: (command: string) => string;
  /** Back to its default key (none: undefined). */
  readonly reset: (command: string, keys: string | undefined) => string;
  readonly resetAll: string;
}

export const DefaultKeysLabels: KeysLabels = {
  recorded: (command, keys, shared) =>
    `${command}: ${keys}${shared.length ? `, also used by ${shared.join(' and ')}` : ''}.`,
  cleared: (command) => `${command}: no key.`,
  reset: (command, keys) => `${command}: ${keys ?? 'no key'} again.`,
  resetAll: 'Every key is back to its default.',
};

export interface KeysModelOptions {
  readonly labels?: Partial<KeysLabels>;
}

/** The keys Help lists and changes: the commands' keys after the keymap, conflicts, and the edit in progress. */
export class KeysModel<Ctx> {
  readonly #manager: CommandManager<Ctx>;
  readonly #labels: KeysLabels;
  readonly #editing = new State(false);
  readonly #recording = new State<string | undefined>(undefined);
  readonly #said = new State('');

  /** Commands with a key, those the person remapped or cleared (to reset them), and the one being recorded. */
  readonly listed: ReadonlyState<readonly CommandDefinition<Ctx>[]>;
  /** Commands on one key in one scope (commands §6.4). */
  readonly conflicts: ReadonlyState<readonly KeyConflict<Ctx>[]>;
  /** Commands with no key now, by title: what a key can be given to. */
  readonly keyless: ReadonlyState<readonly CommandDefinition<Ctx>[]>;

  constructor(manager: CommandManager<Ctx>, options: KeysModelOptions = {}) {
    this.#manager = manager;
    this.#labels = { ...DefaultKeysLabels, ...options.labels };
    const keymap = manager.keymap;
    this.listed = derive([manager.commands, keymap.layers, this.#recording], (all, _, rec) =>
      all.filter((c) => manager.shortcut(c) || keymap.isRemapped(c.id) || c.id === rec),
    );
    this.conflicts = derive([manager.commands, keymap.layers], (all) => keymap.conflicts(all), sameList);
    this.keyless = derive([manager.commands, keymap.layers], (all) =>
      all.filter((c) => !manager.shortcut(c)).sort((a, b) => a.title.localeCompare(b.title)),
    );
  }

  get manager(): CommandManager<Ctx> {
    return this.#manager;
  }
  /** Changing keys: each listed command offers Change, Clear and Reset. */
  get editing(): ReadonlyState<boolean> {
    return this.#editing.asReadonly;
  }
  /** The command whose new key is being recorded. */
  get recording(): ReadonlyState<string | undefined> {
    return this.#recording.asReadonly;
  }
  /**
   * What the last change did, to be said once (a polite live region): the same text again (Reset all twice) toggles a
   * trailing no-break space, so the region is told and says it again (as `StatusModel.spoken`).
   */
  get said(): ReadonlyState<string> {
    return this.#said.asReadonly;
  }

  setEditing(on: boolean): void {
    this.#editing.set(on);
    this.#recording.set(undefined);
  }

  /** Waits for the new key of a command (one without a key too). */
  record(id: string): void {
    if (this.#manager.get(id)) this.#recording.set(id);
  }

  /** Stops waiting, changing nothing. */
  cancel(): void {
    this.#recording.set(undefined);
  }

  /**
   * The recorded key for the command being recorded (or `id`): remapped and saved, said with the commands it now
   * shares the key with.
   */
  recorded(chord: KeyChord, id = this.#recording.value): void {
    const c = id ? this.#manager.get(id) : undefined;
    this.#recording.set(undefined);
    if (!c) return;
    const keymap = this.#manager.keymap;
    keymap.remap(c.id, chord);
    const others = keymap.sharing(c.id, chord, this.#manager.commands.value);
    this.#say(
      this.#labels.recorded(
        c.title,
        formatChord(chord),
        others.map((o) => o.title),
      ),
    );
  }

  /** Takes the command's key away (the person's `null`). */
  clear(id: string): void {
    const c = this.#manager.get(id);
    if (!c) return;
    this.#manager.keymap.remap(id, null);
    this.#say(this.#labels.cleared(c.title));
  }

  /** Back to the app's or the command's own key. */
  reset(id: string): void {
    const c = this.#manager.get(id);
    if (!c) return;
    const keymap = this.#manager.keymap;
    keymap.reset(id);
    const d = keymap.defaultShortcut(c);
    this.#say(this.#labels.reset(c.title, d ? formatChord(d) : undefined));
  }

  /** Every key back to its default. */
  resetAll(): void {
    this.#manager.keymap.resetAll();
    this.#say(this.#labels.resetAll);
  }

  #say(text: string) {
    this.#said.set(text === this.#said.value ? `${text}\u00a0` : text);
  }
}

/** Conflicts as a list: the same chords on the same commands are the same list (nothing is told). */
function sameList<Ctx>(a: readonly KeyConflict<Ctx>[], b: readonly KeyConflict<Ctx>[]): boolean {
  return (
    a.length === b.length &&
    a.every(
      (x, i) =>
        x.scope === b[i].scope &&
        chordId(x.chord) === chordId(b[i].chord) &&
        x.commands.length === b[i].commands.length &&
        x.commands.every((c, j) => c === b[i].commands[j]),
    )
  );
}
