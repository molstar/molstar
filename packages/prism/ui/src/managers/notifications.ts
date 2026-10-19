/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Observable, Subject } from 'rxjs';
import { observableState } from '../reactive/observable-state.js';
import { ReactiveModel } from '../reactive/reactive-model.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

export type MessageKind = 'info' | 'success' | 'warning' | 'error';

export interface StatusMessage {
  id: number;
  kind: MessageKind;
  text: string;
  time: Date;
  /** Said in the log only, beside the text (how long it took: "120 ms"); the status line never shows it. */
  detail?: string;
  /**
   * What wrote a log-only entry (`log`: the plugin's transforms, for one); absent: a message the status line showed. The
   * Log tab filters by it.
   */
  source?: string;
}

export interface LogOptions {
  readonly detail?: string;
  /** Entries kept of this source, the newest (default: the history's size). */
  readonly keep?: number;
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
 * running tasks, and chips (selection, focus). Keeps a bounded history for the log overlay. The viewer's
 * `managers.notifications` (api.md §11, prism-viewer README decision 22: `current`, `tasks`, `chips`; `posted`).
 */
export class NotificationsManager extends ReactiveModel {
  readonly #current = observableState<StatusMessage | undefined>(undefined);
  readonly #history = observableState<readonly StatusMessage[]>([]);
  readonly #tasks = observableState<readonly StatusTask[]>([]);
  /** The tasks put one by one (`putTask`), by id: `setTasks` leaves them. */
  readonly #put = new Set<string>();
  readonly #chips = observableState<readonly StatusChip[]>([]);
  readonly #posted = new Subject<StatusMessage>();
  readonly #spoken = observableState('');
  #serial = 0;
  #fade: ReturnType<typeof setTimeout> | undefined;
  /** Log-only entries waiting for the end of the task: a burst (a load's hundreds of transforms) is one change. */
  #pending: { message: StatusMessage; keep: number | undefined }[] = [];
  /** The sources with their own bound (`log`'s `keep`). */
  readonly #keep = new Map<string, number>();

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

  /**
   * What the shell's one polite live region says (accessibility.md): the latest message posted or announced. A repeat
   * of the same words differs by a trailing no-break space, so it is said again.
   */
  get spoken(): ReadonlyState<string> {
    return this.#spoken.asReadonly;
  }

  /**
   * Says `text` to screen readers only: the live region, with no message on show and nothing in the history (a
   * selection's or a focus's result, the assistant's finished reply).
   */
  announce(text: string): void {
    if (!text) return;
    this.#spoken.set(text === this.#spoken.value ? `${text}\u00a0` : text);
  }

  /** Each message as it is posted (not replayed). */
  get posted(): Observable<StatusMessage> {
    return this.#posted.asObservable();
  }

  post(text: string, kind: MessageKind = 'info', options: { readonly detail?: string } = {}): void {
    const message: StatusMessage = {
      id: ++this.#serial,
      kind,
      text,
      time: new Date(),
      ...(options.detail ? { detail: options.detail } : {}),
    };
    this.#flush();
    this.#history.set(this.#bounded([message, ...this.#history.value]));
    this.#current.set(message);
    this.announce(text);
    clearTimeout(this.#fade);
    if (kind !== 'error') this.#fade = setTimeout(() => this.dismiss(message.id), this.options.fadeMs ?? 4000);
    this.#posted.next(message);
  }

  /**
   * Writes `text` to the log only (the Log tab, filtered by `source`): never on show, never said. Entries of one source
   * past `keep` drop the oldest of that source. Written once a task, so a burst is one change of the history.
   */
  log(text: string, kind: MessageKind, source: string, options: LogOptions = {}): void {
    const message: StatusMessage = {
      id: ++this.#serial,
      kind,
      text,
      time: new Date(),
      source,
      ...(options.detail ? { detail: options.detail } : {}),
    };
    if (!this.#pending.length) queueMicrotask(() => this.#flush());
    this.#pending.push({ message, keep: options.keep });
  }

  #flush() {
    if (!this.#pending.length) return;
    const pending = this.#pending;
    this.#pending = [];
    for (const { message, keep } of pending) if (keep !== undefined) this.#keep.set(message.source!, keep);
    const next: StatusMessage[] = [];
    for (let i = pending.length - 1; i >= 0; i--) next.push(pending[i].message);
    for (const m of this.#history.value) next.push(m);
    this.#history.set(this.#bounded(next));
  }

  /**
   * The history within its bounds, newest kept: a source with its own (`keep`) to that many, the rest (the messages)
   * to the history's size, so a burst of transforms does not push the messages out.
   */
  #bounded(list: StatusMessage[]): StatusMessage[] {
    const size = this.options.historySize ?? 200;
    if (!this.#keep.size) return list.length > size ? list.slice(0, size) : list;
    const counts = new Map<string | undefined, number>();
    return list.filter((m) => {
      const own = m.source === undefined ? undefined : this.#keep.get(m.source);
      const key = own === undefined ? undefined : m.source;
      const n = (counts.get(key) ?? 0) + 1;
      counts.set(key, n);
      return n <= (own ?? size);
    });
  }

  /** Adds or replaces a task by id; it stays until `removeTask`, whatever `setTasks` sets. */
  putTask(task: StatusTask): void {
    this.#put.add(task.id);
    const tasks = this.#tasks.value;
    const i = tasks.findIndex((t) => t.id === task.id);
    this.#tasks.set(i < 0 ? [...tasks, task] : tasks.map((t, j) => (j === i ? task : t)));
  }

  removeTask(id: string): void {
    this.#put.delete(id);
    if (this.#tasks.value.some((t) => t.id === id)) this.#tasks.set(this.#tasks.value.filter((t) => t.id !== id));
  }

  /** Hides the message on show (only that one, if `id` is given). */
  dismiss(id?: number): void {
    if (id !== undefined && this.#current.value?.id !== id) return;
    clearTimeout(this.#fade);
    this.#current.set(undefined);
  }

  /**
   * Sets the tasks one source mirrors (the plugin's): every task but those put one by one (`putTask`), which stay until
   * they are removed.
   */
  setTasks(tasks: readonly StatusTask[]): void {
    const put = this.#tasks.value.filter((t) => this.#put.has(t.id) && !tasks.some((n) => n.id === t.id));
    this.#tasks.set(put.length ? [...tasks, ...put] : tasks);
  }

  setChips(chips: readonly StatusChip[]): void {
    this.#chips.set(chips);
  }

  protected onRetire(): void {
    clearTimeout(this.#fade);
    this.#posted.complete();
  }
}
