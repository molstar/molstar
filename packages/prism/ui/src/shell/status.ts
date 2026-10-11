/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observableState } from '../reactive/observable-state.js';
import { ReactiveModel } from '../reactive/reactive-model.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

export type MessageKind = 'info' | 'success' | 'warning' | 'error';

export interface StatusMessage {
  id: number;
  kind: MessageKind;
  text: string;
  time: Date;
}

export interface StatusTask {
  id: string;
  label: string;
  /** 0..1, or undefined when indeterminate. */
  progress?: number;
  /**
   * It holds the whole plugin while it runs (rendering an animation, exporting): shown over everything, with its
   * progress and Cancel, instead of in the status line. Mol*'s task "overlay" mode.
   */
  overlay?: boolean;
  abort?(): void;
}

export interface StatusChip {
  id: string;
  label: string;
  /** Clicking the chip (e.g. focus the selection). */
  activate?(): void;
  /** The chip's ✕ (e.g. clear the selection). */
  dismiss?(): void;
}

/**
 * What the status line on the canvas shows: the latest message (it fades after a while; errors stay until dismissed),
 * running tasks, and chips (selection, focus). Keeps a bounded history for the log overlay.
 */
export class StatusModel extends ReactiveModel {
  readonly #current = observableState<StatusMessage | undefined>(undefined);
  readonly #history = observableState<readonly StatusMessage[]>([]);
  readonly #tasks = observableState<readonly StatusTask[]>([]);
  readonly #chips = observableState<readonly StatusChip[]>([]);
  #serial = 0;
  #fade: ReturnType<typeof setTimeout> | undefined;

  constructor(readonly options: { fadeMs?: number; historySize?: number } = {}) {
    super();
  }

  /** The message on show, if any. */
  get current(): ReadonlyState<StatusMessage | undefined> {
    return this.#current.asReadonly;
  }

  /** Recent messages, newest first. */
  get history(): ReadonlyState<readonly StatusMessage[]> {
    return this.#history.asReadonly;
  }

  get tasks(): ReadonlyState<readonly StatusTask[]> {
    return this.#tasks.asReadonly;
  }

  get chips(): ReadonlyState<readonly StatusChip[]> {
    return this.#chips.asReadonly;
  }

  post(text: string, kind: MessageKind = 'info'): void {
    const message: StatusMessage = { id: ++this.#serial, kind, text, time: new Date() };
    this.#history.set([message, ...this.#history.value].slice(0, this.options.historySize ?? 200));
    this.#current.set(message);
    clearTimeout(this.#fade);
    if (kind !== 'error') this.#fade = setTimeout(() => this.dismiss(message.id), this.options.fadeMs ?? 4000);
  }

  /** Hides the message on show (only that one, if `id` is given). */
  dismiss(id?: number): void {
    if (id !== undefined && this.#current.value?.id !== id) return;
    clearTimeout(this.#fade);
    this.#current.set(undefined);
  }

  setTasks(tasks: readonly StatusTask[]): void {
    this.#tasks.set(tasks);
  }

  setChips(chips: readonly StatusChip[]): void {
    this.#chips.set(chips);
  }

  protected onRetire(): void {
    clearTimeout(this.#fade);
  }
}
