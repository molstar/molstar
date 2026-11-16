/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Which key runs which command (commands §6): the keymap's three layers, and the per-key index a keydown is looked up
// in. DOM-free: targets are read through `getAttribute`, `tagName`, `parentNode` where they have them.

import { type ReadonlyState, State } from '@molstar/core/reactive';
import { chordId, chordMatches, codeKey, isPlainChord, type KeyEventLike, typesNothing } from './keys.js';
import type { CommandDefinition, KeyBindings, KeyChord, KeymapStore } from './types.js';

/** The attribute that names a key scope on an element around the focus (commands §6.2). */
export const ScopeAttribute = 'data-msui-scope';

export interface KeymapOptions {
  /** The app's keys: replace the commands' own. */
  app?: KeyBindings;
  /** Where the person's remaps are kept; without one they last as long as the keymap. */
  store?: KeymapStore | undefined;
}

export interface KeymapLayers {
  readonly app: KeyBindings;
  readonly user: KeyBindings;
}

/** Commands on one key in one scope: only the first that can run takes it. */
export interface KeyConflict<Ctx> {
  chord: KeyChord;
  scope: string;
  commands: readonly CommandDefinition<Ctx>[];
}

type Keyed = Pick<CommandDefinition<unknown>, 'id' | 'shortcut'>;

/** Where a command's key works: `global` unless it names a scope. */
export function scopeOf(command: { scope?: string }): string {
  return command.scope && command.scope !== 'global' ? command.scope : 'global';
}

/** A stored value read back: only chords (an object with a non-empty string `key`) and `null`, by string ids. */
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

/** The keymap (commands §6.3): the person's remaps over the app's over the commands' own `shortcut`s. */
export class Keymap {
  readonly #store: KeymapStore | undefined;
  readonly #layers: State<KeymapLayers>;

  constructor(options: KeymapOptions = {}) {
    this.#store = options.store;
    let user: KeyBindings = {};
    try {
      user = readBindings(options.store?.load());
    } catch {
      user = {};
    }
    this.#layers = new State<KeymapLayers>({ app: readBindings(options.app ?? {}), user });
  }

  /** The app's and the person's layers; a new value on every change (keys shown anywhere follow it). */
  get layers(): ReadonlyState<KeymapLayers> {
    return this.#layers.asReadonly;
  }

  /** Replaces the app's layer. */
  setApp(bindings: KeyBindings): void {
    this.#layers.set({ ...this.#layers.value, app: readBindings(bindings) });
  }

  /** The command's key now: the person's, else the app's, else its own; none when a layer took it away. */
  shortcut(command: Keyed): KeyChord | undefined {
    const { app, user } = this.#layers.value;
    if (Object.hasOwn(user, command.id)) return user[command.id] ?? undefined;
    if (Object.hasOwn(app, command.id)) return app[command.id] ?? undefined;
    return command.shortcut;
  }

  /** The key the command has without the person's remap (what Reset goes back to). */
  defaultShortcut(command: Keyed): KeyChord | undefined {
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

  /** Commands sharing a key in one scope; a scoped command over a global one is not a conflict (commands §6.4). */
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

// --- The per-key index (commands §6.1) ------------------------------------------------------------

interface KeyEntry<Ctx> {
  command: CommandDefinition<Ctx>;
  chord: KeyChord;
  scope: string;
  /** Its place in the list of commands: the order among hits of one rank. */
  order: number;
}

/** What a keydown found: the command to open, and its chord. */
export interface KeyHit<Ctx> {
  command: CommandDefinition<Ctx>;
  chord: KeyChord;
}

/**
 * The commands by their key (`chord.key` in lower case, what `chordMatches` compares the event's key and code with),
 * and the app's and person's chords the same way; made again only when the list of commands or the layers change.
 * Why: every keydown in the app comes here, and a lookup must check only the commands on its key and allocate nothing.
 */
export class KeyIndex<Ctx> {
  readonly #keymap: Keymap;
  #of: readonly CommandDefinition<Ctx>[] | undefined;
  #layers: KeymapLayers | undefined;
  #byKey = new Map<string, KeyEntry<Ctx>[]>();
  #bound = new Map<string, KeyChord[]>();

  constructor(keymap: Keymap) {
    this.#keymap = keymap;
  }

  #at(commands: readonly CommandDefinition<Ctx>[]): void {
    const layers = this.#keymap.layers.value;
    if (this.#of === commands && this.#layers === layers) return;
    const byKey = new Map<string, KeyEntry<Ctx>[]>();
    for (let order = 0; order < commands.length; order++) {
      const command = commands[order];
      const chord = this.#keymap.shortcut(command);
      if (chord) push(byKey, chord.key.toLowerCase(), { command, chord, scope: scopeOf(command), order });
    }
    const bound = new Map<string, KeyChord[]>();
    for (const layer of [layers.user, layers.app])
      for (const id in layer) {
        const chord = layer[id];
        if (chord) push(bound, chord.key.toLowerCase(), chord);
      }
    this.#of = commands;
    this.#layers = layers;
    this.#byKey = byKey;
    this.#bound = bound;
  }

  /** Does a chord of the app's or the person's layer match the event (its command may be one not made yet)? */
  bound(commands: readonly CommandDefinition<Ctx>[], e: KeyEventLike, mac: boolean): boolean {
    this.#at(commands);
    const key = e.key.toLowerCase();
    const code = codeKey(e, key);
    return (
      matchesAny(this.#bound.get(key), e, mac) || (code !== undefined && matchesAny(this.#bound.get(code), e, mac))
    );
  }

  /**
   * The command the event runs (commands §6.1–6.2): among those on its key that the target allows (a text field takes
   * only ⌘/Ctrl chords without Alt; a control that steps with arrows keeps plain keys) and whose scope is around it,
   * the nearest scope first, global last, then the first registered; the first that can run, else the nearest.
   */
  find(
    commands: readonly CommandDefinition<Ctx>[],
    e: KeyEventLike,
    mac: boolean,
    blocked: (c: CommandDefinition<Ctx>) => string | undefined,
  ): KeyHit<Ctx> | undefined {
    this.#at(commands);
    const key = e.key.toLowerCase();
    const code = codeKey(e, key);
    const first = this.#byKey.get(key);
    const second = code !== undefined && code !== key ? this.#byKey.get(code) : undefined;
    // Most keys bind nothing: the target and its path are read only for a key that does.
    if (!first && !second) return undefined;
    const path = e.composedPath?.();
    const target = path?.length ? path[0] : e.target;
    const typing = typesText(target);
    const owns = !typing && usesKeys(target);
    let hit: KeyEntry<Ctx> | undefined;
    let hitRank = 0;
    let runnable: KeyEntry<Ctx> | undefined;
    let runnableRank = 0;
    for (let pass = 0; pass < 2; pass++) {
      const bucket = pass === 0 ? first : second;
      if (!bucket) continue;
      for (const at of bucket) {
        const { chord } = at;
        if ((typing && !typesNothing(chord)) || (owns && isPlainChord(chord)) || !chordMatches(chord, e, mac)) continue;
        const rank = at.scope === 'global' ? Number.POSITIVE_INFINITY : scopeRank(path, target, at.scope);
        if (rank < 0) continue;
        if (!hit || before(rank, at, hitRank, hit)) {
          hit = at;
          hitRank = rank;
        }
        if ((!runnable || before(rank, at, runnableRank, runnable)) && !blocked(at.command)) {
          runnable = at;
          runnableRank = rank;
        }
      }
    }
    const found = runnable ?? hit;
    return found ? { command: found.command, chord: found.chord } : undefined;
  }
}

function push<T>(map: Map<string, T[]>, key: string, value: T): void {
  const bucket = map.get(key);
  if (bucket) bucket.push(value);
  else map.set(key, [value]);
}

function matchesAny(chords: readonly KeyChord[] | undefined, e: KeyEventLike, mac: boolean): boolean {
  if (chords) for (const chord of chords) if (chordMatches(chord, e, mac)) return true;
  return false;
}

/** Is `a` at `rank` nearer than `b` at `bRank`: a nearer scope, else registered first. */
function before(rank: number, a: KeyEntry<unknown>, bRank: number, b: KeyEntry<unknown>): boolean {
  return rank < bRank || (rank === bRank && a.order < b.order);
}

type Node = {
  tagName?: string;
  type?: string;
  isContentEditable?: boolean;
  getAttribute?(name: string): string | null;
  parentNode?: unknown;
  host?: unknown;
};

function attribute(node: unknown, name: string): string | null {
  const n = node as Node | null | undefined;
  return n && typeof n.getAttribute === 'function' ? n.getAttribute(name) : null;
}

/**
 * How near the scope is to the target (commands §6.2): the index of the first element naming it on the event's composed
 * path (across shadow roots), or, without a path, among the target and its ancestors; -1 when it is not around.
 */
function scopeRank(path: readonly unknown[] | undefined, target: unknown, scope: string): number {
  if (path?.length) {
    for (let i = 0; i < path.length; i++) if (attribute(path[i], ScopeAttribute) === scope) return i;
    return -1;
  }
  let node = target as Node | null | undefined;
  for (let i = 0; node && i < 256; i++) {
    if (attribute(node, ScopeAttribute) === scope) return i;
    node = (node.parentNode ?? node.host) as Node | null | undefined;
  }
  return -1;
}

const NotText = new Set(['button', 'checkbox', 'radio', 'range', 'color', 'file', 'image', 'reset', 'submit']);

/** Does the target take typed text (a text input, a textarea, a select, an editable element)? */
function typesText(target: unknown): boolean {
  if (!target || typeof target !== 'object') return false;
  const t = target as Node;
  if (t.isContentEditable) return true;
  const editable = attribute(t, 'contenteditable');
  if (editable !== null && editable !== 'false') return true;
  const tag = t.tagName?.toUpperCase();
  if (tag === 'TEXTAREA' || tag === 'SELECT') return true;
  return tag === 'INPUT' && !NotText.has((t.type ?? 'text').toLowerCase());
}

const KeyRoles = new Set(['slider', 'spinbutton', 'radio', 'tab', 'menuitem', 'option', 'treeitem', 'listbox']);

/** Does the target move with the arrow keys itself (a range or radio input, a widget with a role that does)? */
function usesKeys(target: unknown): boolean {
  if (!target || typeof target !== 'object') return false;
  const t = target as Node;
  if (t.tagName?.toUpperCase() === 'INPUT') {
    const type = (t.type ?? '').toLowerCase();
    if (type === 'range' || type === 'radio') return true;
  }
  const role = attribute(t, 'role');
  return role !== null && KeyRoles.has(role);
}
