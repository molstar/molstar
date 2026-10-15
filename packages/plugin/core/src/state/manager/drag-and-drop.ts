/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommands } from '@molstar/plugin/commands';
import type { PluginContext } from '@molstar/plugin/context';

export type PluginDragAndDropHandler = (files: File[], plugin: PluginContext) => Promise<boolean> | boolean;

export interface PluginDragAndDropEntry {
  readonly name: string;
  readonly handle: PluginDragAndDropHandler;
  /** Fallback handlers run after the built-in session handling, so they only get files no other handler took. */
  readonly fallback?: boolean;
}

/**
 * Handlers are keyed by `name`; a provider's identity is its `handle` function.
 * Registering the same name and handle again increments a count, a different handle
 * under an existing name throws.
 */
export class DragAndDropManager {
  private entries: { entry: PluginDragAndDropEntry; count: number }[] = [];

  /** Returns the message `addEntry` would throw for `entry`, without changing anything. */
  findConflict(entry: PluginDragAndDropEntry): string | undefined {
    const existing = this.entries.find((e) => e.entry.name === entry.name);
    if (existing && existing.entry.handle !== entry.handle) {
      return `DragAndDropManager: a different handler is already registered under the name '${entry.name}'.`;
    }
    return undefined;
  }

  addEntry(entry: PluginDragAndDropEntry) {
    const conflict = this.findConflict(entry);
    if (conflict) throw new Error(conflict);

    const existing = this.entries.find((e) => e.entry.name === entry.name);
    if (existing) existing.count++;
    else this.entries.push({ entry, count: 1 });
  }

  /** Decrements the provider registered with the same name and handle, no-op when there is none. */
  removeEntry(entry: PluginDragAndDropEntry) {
    const index = this.entries.findIndex((e) => e.entry.name === entry.name && e.entry.handle === entry.handle);
    if (index >= 0) this.decrement(index);
  }

  addHandler(name: string, handler: PluginDragAndDropHandler, options?: { fallback?: boolean }) {
    this.addEntry({ name, handle: handler, ...options });
  }

  /** Decrements the provider registered under `name`, no-op for an unknown name. */
  removeHandler(name: string) {
    const index = this.entries.findIndex((e) => e.entry.name === name);
    if (index >= 0) this.decrement(index);
  }

  /** The registered entries in registration order. */
  list(): readonly PluginDragAndDropEntry[] {
    return this.entries.map((e) => e.entry);
  }

  private decrement(index: number) {
    if (--this.entries[index].count <= 0) this.entries.splice(index, 1);
  }

  /**
   * Tries non-fallback handlers, then the built-in session handling, then fallback handlers;
   * each group in reverse registration order.
   */
  async handle(files: File[]) {
    if (await this.tryHandlers(files, false)) return;
    if (openSession(this.plugin, files)) return;
    await this.tryHandlers(files, true);
  }

  private async tryHandlers(files: File[], fallback: boolean) {
    // snapshot, handlers may (un)register while running
    const entries = this.entries.filter((e) => !!e.entry.fallback === fallback).map((e) => e.entry);
    for (let i = entries.length - 1; i >= 0; i--) {
      const handled = await entries[i].handle(files, this.plugin);
      if (handled) return true;
    }
    return false;
  }

  dispose() {
    this.entries.length = 0;
  }

  constructor(public plugin: PluginContext) {}
}

function openSession(plugin: PluginContext, files: File[]) {
  const sessions = files.filter((f) => {
    const fn = f.name.toLowerCase();
    return fn.endsWith('.molx') || fn.endsWith('.molj');
  });
  if (sessions.length === 0) return false;

  PluginCommands.State.Snapshots.OpenFile(plugin, { file: sessions[0] });
  return true;
}
