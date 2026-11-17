/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What is registered (commands §4.1–4.2): the commands, the groups they browse in, who registered each, and the
// commands made when wanted. The manager holds one and adds what runs.

import { batch, type ReadonlyState, State } from '@molstar/core/reactive';
import type { CommandDefinition, CommandGroup } from './types.js';

/**
 * Commands made when they are wanted, not at every change of what they are made from (commands §4.2). `make` runs on
 * the first read after `stale()`, or in a microtask while the lists are wanted (`want`). A make that throws is logged,
 * not thrown to the read.
 */
export interface LazyCommands<Ctx> {
  /** The commands and groups now; `undefined` keeps what was made last (nothing they list changed). */
  make(): { commands: readonly CommandDefinition<Ctx>[]; groups: readonly CommandGroup[] } | undefined;
}

export interface LazyRegistration {
  /** What the commands are made from changed: made again when next wanted. */
  stale(): void;
  /** Takes them away. */
  dispose(): void;
}

interface Lazy<Ctx> {
  readonly source: LazyCommands<Ctx>;
  stale: boolean;
  /** The ids it made last: a `get` of one of them makes it again when stale. */
  ids: ReadonlySet<string>;
  drop?: () => void;
}

/** An index over a list, made once per list (a `find` per read was the list's length per call). */
function indexById<T extends { id: string }>(all: readonly T[]): Map<string, T> {
  const byId = new Map<string, T>();
  for (const x of all) if (!byId.has(x.id)) byId.set(x.id, x);
  return byId;
}

export class CommandRegistry<Ctx> {
  readonly #commands = new State<readonly CommandDefinition<Ctx>[]>([]);
  readonly #groups = new State<readonly CommandGroup[]>([]);
  /** Each command's registration (a token per `register` call): a run reports only while it lasts. */
  readonly #owners = new Map<CommandDefinition<Ctx>, object>();
  readonly #lazy = new Set<Lazy<Ctx>>();
  /** Holders of `want`: stale lazy commands are made in a microtask, not by the reads. */
  #wanted = 0;
  /** How many lazy contributions are stale: what a keydown checks, a number rather than a walk of the set. */
  #staleCount = 0;
  #making = false;
  #soon = false;
  #commandIndex: { of: readonly CommandDefinition<Ctx>[]; byId: Map<string, CommandDefinition<Ctx>> } | undefined;
  #groupIndex: { of: readonly CommandGroup[]; byId: Map<string, CommandGroup> } | undefined;

  /** Told before commands go (their run is cancelled, their invocation closed). */
  readonly #onRemove: (gone: ReadonlySet<CommandDefinition<Ctx>>, owner: object) => void;

  constructor(onRemove: (gone: ReadonlySet<CommandDefinition<Ctx>>, owner: object) => void) {
    this.#onRemove = onRemove;
  }

  get commands(): ReadonlyState<readonly CommandDefinition<Ctx>[]> {
    return this.#commands.asReadonly;
  }

  get groups(): ReadonlyState<readonly CommandGroup[]> {
    return this.#groups.asReadonly;
  }

  /** The registration a command came with, while it lasts. */
  ownerOf(command: CommandDefinition<Ctx>): object | undefined {
    return this.#owners.get(command);
  }

  /**
   * Adds commands, as one change of the list; the returned function takes them away again (one change too). An id
   * already registered, or repeated in `defs`, throws, naming it, and nothing is added.
   */
  register(defs: readonly CommandDefinition<Ctx>[]): () => void {
    const owner = {};
    const taken = new Set(this.#commands.value.map((c) => c.id));
    for (const d of defs) {
      if (taken.has(d.id)) throw new Error(`Command '${d.id}' is already registered.`);
      taken.add(d.id);
    }
    for (const d of defs) this.#owners.set(d, owner);
    this.#commands.set([...this.#commands.value, ...defs]);
    let done = false;
    return () => {
      if (done) return;
      done = true;
      // A set, not `defs.includes` per command (taking D away from C commands was C × D).
      const gone = new Set<CommandDefinition<Ctx>>();
      for (const d of defs) if (this.#owners.get(d) === owner) gone.add(d);
      if (!gone.size) return;
      this.#onRemove(gone, owner);
      for (const d of gone) this.#owners.delete(d);
      this.#commands.set(this.#commands.value.filter((c) => !gone.has(c)));
    };
  }

  /** Adds groups, as one change; the returned function takes them away. A taken or repeated id throws, adding nothing. */
  registerGroups(groups: readonly CommandGroup[]): () => void {
    const taken = new Set(this.#groups.value.map((g) => g.id));
    for (const g of groups) {
      if (taken.has(g.id)) throw new Error(`Command group '${g.id}' is already registered.`);
      taken.add(g.id);
    }
    if (groups.length) this.#groups.set([...this.#groups.value, ...groups]);
    let done = false;
    return () => {
      if (done || !groups.length) return;
      done = true;
      const these = new Set(groups);
      this.#groups.set(this.#groups.value.filter((g) => !these.has(g)));
    };
  }

  /** Adds commands made when wanted: stale from the start, made by the first read that needs them (commands §4.2). */
  registerLazy(source: LazyCommands<Ctx>): LazyRegistration {
    const lazy: Lazy<Ctx> = { source, stale: true, ids: new Set() };
    this.#lazy.add(lazy);
    this.#staleCount++;
    if (this.#wanted) this.#makeSoon();
    return {
      stale: () => {
        if (!this.#lazy.has(lazy)) return;
        if (!lazy.stale) this.#staleCount++;
        lazy.stale = true;
        if (this.#wanted) this.#makeSoon();
      },
      dispose: () => {
        if (!this.#lazy.delete(lazy)) return;
        if (lazy.stale) this.#staleCount--;
        lazy.drop?.();
      },
    };
  }

  /** The lists are wanted until the returned release: stale lazy commands are made now, then in a microtask per change. */
  want(): () => void {
    this.#wanted++;
    this.#make();
    let released = false;
    return () => {
      if (released) return;
      released = true;
      this.#wanted--;
    };
  }

  /** Before a read of the lists: while they are wanted the microtask makes them, so a read inside a render does not. */
  fresh(): void {
    if (!this.#wanted) this.#make();
  }

  /** Is a lazy contribution stale (a bound key may name one of its commands)? */
  get anyStale(): boolean {
    return this.#staleCount > 0;
  }

  get(id: string): CommandDefinition<Ctx> | undefined {
    // Made only when this id may be a stale lazy one: a `get` of another command makes nothing.
    if (this.#mayMake(id)) this.#make();
    return this.byId().get(id);
  }

  /** The commands by id, made once per list. */
  byId(): Map<string, CommandDefinition<Ctx>> {
    const all = this.#commands.value;
    if (this.#commandIndex?.of !== all) this.#commandIndex = { of: all, byId: indexById(all) };
    return this.#commandIndex.byId;
  }

  /** The groups by id, made once per list. */
  groupsById(): Map<string, CommandGroup> {
    const all = this.#groups.value;
    if (this.#groupIndex?.of !== all) this.#groupIndex = { of: all, byId: indexById(all) };
    return this.#groupIndex.byId;
  }

  /** A group's path from the top level to it (itself last); empty for none or an unknown id. */
  groupPath(id: string | undefined): CommandGroup[] {
    this.fresh();
    const byId = this.groupsById();
    const path: CommandGroup[] = [];
    let g = id === undefined ? undefined : byId.get(id);
    // A cycle of parents stops at the first group seen again.
    while (g && !path.includes(g)) {
      path.push(g);
      g = g.parent === undefined ? undefined : byId.get(g.parent);
    }
    return path.reverse();
  }

  /** What a group holds (or the top level, for none): its subgroups by `order`, and its commands in order. */
  children(id: string | undefined): { groups: CommandGroup[]; commands: CommandDefinition<Ctx>[] } {
    this.fresh();
    const known = this.groupsById();
    const parentOf = (g: string | undefined) => (g !== undefined && known.has(g) ? g : undefined);
    const groups = this.#groups.value
      .filter((g) => parentOf(g.parent) === id)
      .sort((a, b) => (a.order ?? 0) - (b.order ?? 0));
    const commands = this.#commands.value.filter((c) => parentOf(c.group) === id);
    return { groups, commands };
  }

  #makeSoon(): void {
    if (this.#soon) return;
    this.#soon = true;
    queueMicrotask(() => {
      this.#soon = false;
      this.#make();
    });
  }

  /** Makes the stale lazy commands; each contribution's swap is one change of the lists. */
  #make(): void {
    if (this.#making || !this.#staleCount) return;
    this.#making = true;
    try {
      for (const lazy of this.#lazy) {
        if (!lazy.stale) continue;
        lazy.stale = false;
        this.#staleCount--;
        // A make that throws (or an id already taken) is said here, not to the read that asked: the lists stay as they
        // are, and the next change tries again.
        try {
          const made = lazy.source.make();
          if (!made) continue;
          batch(() => {
            lazy.drop?.();
            lazy.drop = undefined;
            lazy.ids = new Set();
            const dropGroups = this.registerGroups(made.groups);
            let dropCommands: () => void;
            try {
              dropCommands = this.register(made.commands);
            } catch (e) {
              dropGroups();
              throw e;
            }
            lazy.drop = () => batch(() => (dropCommands(), dropGroups()));
            lazy.ids = new Set(made.commands.map((c) => c.id));
          });
        } catch (e) {
          console.error('Lazy commands could not be made:', e);
        }
      }
    } finally {
      this.#making = false;
    }
  }

  /** Whether a stale lazy contribution may hold `id`: one it made, or an id not registered now. */
  #mayMake(id: string): boolean {
    if (this.#wanted || this.#making) return false;
    let unknown: boolean | undefined;
    for (const lazy of this.#lazy) {
      if (!lazy.stale) continue;
      if (lazy.ids.has(id)) return true;
      unknown ??= !this.byId().has(id);
      if (unknown) return true;
    }
    return false;
  }
}
