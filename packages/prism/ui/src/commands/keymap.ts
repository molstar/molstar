/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The keymap (commands.md §6): which key runs which command, in three layers. The person's remaps over the app's
// overrides over the commands' own `shortcut`s; the person's are saved through a store the host gives (the viewer: its
// preferences). DOM-free.

import { type ObservableState, observableState } from '../reactive/observable-state.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';
import { chordId } from './keys.js';
import type { CommandDefinition, KeyBindings, KeyChord, KeymapStore } from './types.js';

export interface KeymapOptions {
  /** The app's keys (the viewer's `spec.ui.keybindings`): replace the commands' own. */
  app?: KeyBindings;
  /** Where the person's remaps are kept; without one they last as long as the keymap. */
  store?: KeymapStore | undefined;
}

/** Commands on one key in one scope: only the first that can run takes it. */
export interface KeyConflict<Ctx> {
  chord: KeyChord;
  scope: string;
  commands: readonly CommandDefinition<Ctx>[];
}

/** Where a command's key works: `global` unless it names a scope. */
export function scopeOf(command: { scope?: string }): string {
  return command.scope && command.scope !== 'global' ? command.scope : 'global';
}

/** A stored value read back: only chords (an object with a string `key`) and `null`, by string ids. */
function readBindings(value: unknown): KeyBindings {
  if (!value || typeof value !== 'object' || Array.isArray(value)) return {};
  const out: Record<string, KeyChord | null> = {};
  for (const [id, v] of Object.entries(value)) {
    if (v === null) out[id] = null;
    else if (v && typeof v === 'object' && typeof (v as KeyChord).key === 'string' && (v as KeyChord).key) {
      const c = v as KeyChord;
      out[id] = {
        key: c.key,
        ...(c.mod ? { mod: true } : {}),
        ...(c.alt ? { alt: true } : {}),
        ...(c.shift ? { shift: true } : {}),
      };
    }
  }
  return out;
}

export class Keymap {
  readonly #store: KeymapStore | undefined;
  readonly #layers: ObservableState<{ app: KeyBindings; user: KeyBindings }>;

  constructor(options: KeymapOptions = {}) {
    this.#store = options.store;
    let user: KeyBindings = {};
    try {
      user = readBindings(options.store?.load());
    } catch {
      user = {};
    }
    this.#layers = observableState<{ app: KeyBindings; user: KeyBindings }>({
      app: readBindings(options.app ?? {}),
      user,
    });
  }

  /** The app's and the person's layers; emits on every change (Help and the menu redraw the keys). */
  get layers(): ReadonlyState<{ readonly app: KeyBindings; readonly user: KeyBindings }> {
    return this.#layers.asReadonly;
  }

  /** Replaces the app's layer (the viewer's spec). */
  setApp(bindings: KeyBindings): void {
    this.#layers.set({ ...this.#layers.value, app: readBindings(bindings) });
  }

  /** The command's key now: the person's, else the app's, else its own; none when a layer took it away. */
  shortcut(command: Pick<CommandDefinition<unknown>, 'id' | 'shortcut'>): KeyChord | undefined {
    const { app, user } = this.#layers.value;
    if (Object.hasOwn(user, command.id)) return user[command.id] ?? undefined;
    if (Object.hasOwn(app, command.id)) return app[command.id] ?? undefined;
    return command.shortcut;
  }

  /** The key the command has without the person's remap (what Reset goes back to). */
  defaultShortcut(command: Pick<CommandDefinition<unknown>, 'id' | 'shortcut'>): KeyChord | undefined {
    const { app } = this.#layers.value;
    if (Object.hasOwn(app, command.id)) return app[command.id] ?? undefined;
    return command.shortcut;
  }

  /** Did the person remap (or clear) this command's key? */
  isRemapped(id: string): boolean {
    return Object.hasOwn(this.#layers.value.user, id);
  }

  /** The person's key for a command (`null`: none), saved. */
  remap(id: string, chord: KeyChord | null): void {
    this.#setUser({ ...this.#layers.value.user, [id]: chord });
  }

  /** Back to the app's or the command's own key. */
  reset(id: string): void {
    if (!this.isRemapped(id)) return;
    const { [id]: _, ...rest } = this.#layers.value.user;
    this.#setUser(rest);
  }

  resetAll(): void {
    this.#setUser({});
  }

  /**
   * Commands sharing a key in one scope (commands.md §6): a scoped command over a global one on the same key is not a
   * conflict (it wins inside its scope, by design).
   */
  conflicts<Ctx>(commands: readonly CommandDefinition<Ctx>[]): KeyConflict<Ctx>[] {
    const by = new Map<string, { chord: KeyChord; scope: string; commands: CommandDefinition<Ctx>[] }>();
    for (const c of commands) {
      const chord = this.shortcut(c);
      if (!chord) continue;
      const scope = scopeOf(c);
      const k = `${scope} ${chordId(chord)}`;
      const at = by.get(k);
      if (at) at.commands.push(c);
      else by.set(k, { chord, scope, commands: [c] });
    }
    return [...by.values()].filter((x) => x.commands.length > 1);
  }

  /** The other commands that a key would share with this one, in its scope. */
  sharing<Ctx>(id: string, chord: KeyChord, commands: readonly CommandDefinition<Ctx>[]): CommandDefinition<Ctx>[] {
    const self = commands.find((c) => c.id === id);
    const scope = self ? scopeOf(self) : 'global';
    const want = chordId(chord);
    return commands.filter((c) => {
      if (c.id === id || scopeOf(c) !== scope) return false;
      const k = this.shortcut(c);
      return !!k && chordId(k) === want;
    });
  }

  #setUser(user: KeyBindings): void {
    this.#layers.set({ ...this.#layers.value, user });
    try {
      this.#store?.save(user);
    } catch {
      // A store that cannot save keeps nothing; the remap holds for this session.
    }
  }
}
