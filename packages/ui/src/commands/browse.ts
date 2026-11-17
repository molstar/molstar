/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The menu's highlighted row (commands §4.5): previewed and counted once the highlight rests, for that row only.

import { type ReadonlyState, State } from '@molstar/core/reactive';
import type { CommandDefinition, PreviewSink } from './types.js';

/** The menu's highlighted row: its command, and how much it would do once counted. */
export interface CommandBrowse {
  id: string;
  count?: string;
}

export interface BrowseHost<Ctx> {
  context: ReadonlyState<Ctx>;
  preview: PreviewSink | undefined;
  /** Ms the highlight rests before the row is previewed and counted. */
  delay: number;
  blocked(command: CommandDefinition<Ctx>): string | undefined;
  /** What a command without params runs with. */
  args(command: CommandDefinition<Ctx>, ctx: Ctx): Record<string, unknown>;
}

interface Browse<Ctx> {
  command: CommandDefinition<Ctx>;
  abort: AbortController;
  timer?: ReturnType<typeof setTimeout>;
}

export class CommandBrowser<Ctx> {
  readonly #host: BrowseHost<Ctx>;
  readonly #browsing = new State<CommandBrowse | undefined>(undefined);
  #at: Browse<Ctx> | undefined;
  #previewed = false;

  constructor(host: BrowseHost<Ctx>) {
    this.#host = host;
  }

  get browsing(): ReadonlyState<CommandBrowse | undefined> {
    return this.#browsing.asReadonly;
  }

  /** The command highlighted now, if any. */
  get command(): CommandDefinition<Ctx> | undefined {
    return this.#at?.command;
  }

  /**
   * Highlights a command: one without params that can run is previewed and counted after the delay; a move of the
   * highlight (or none) aborts the count's signal and clears the preview first. Why the delay: moving through the list
   * computes only where it rests.
   */
  browse(command: CommandDefinition<Ctx> | undefined): void {
    if (command && command === this.#at?.command) return;
    this.end();
    if (!command) return;
    const at: Browse<Ctx> = { command, abort: new AbortController() };
    this.#at = at;
    this.#browsing.set({ id: command.id });
    const host = this.#host;
    if (command.params || (!command.preview && !command.count) || host.blocked(command)) return;
    if (host.delay > 0) at.timer = setTimeout(() => void this.#start(at), host.delay);
    else void this.#start(at);
  }

  /** No row is highlighted any more: its count is aborted, its preview cleared. */
  end(): void {
    const at = this.#at;
    if (!at) return;
    this.#at = undefined;
    if (at.timer !== undefined) clearTimeout(at.timer);
    at.abort.abort();
    if (this.#previewed) {
      this.#previewed = false;
      this.#host.preview?.clear();
    }
    this.#browsing.set(undefined);
  }

  async #start(at: Browse<Ctx>): Promise<void> {
    at.timer = undefined;
    const { command, abort } = at;
    const host = this.#host;
    const ctx = host.context.value;
    let args: Record<string, unknown>;
    try {
      args = host.args(command, ctx);
    } catch {
      return;
    }
    if (command.preview && host.preview) {
      try {
        const preview = command.preview(args, ctx);
        if (preview !== undefined) {
          host.preview.set(preview);
          this.#previewed = true;
        }
      } catch {
        // A preview that fails shows nothing.
      }
    }
    if (!command.count) return;
    try {
      const count = await command.count(args, ctx, abort.signal);
      if (!abort.signal.aborted && count) this.#browsing.set({ id: command.id, count });
    } catch {
      // A count that fails says nothing.
    }
  }
}
