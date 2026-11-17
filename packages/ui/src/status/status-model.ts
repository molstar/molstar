/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { batch, KeyedState, type ReadonlyState, ReactiveModel, State } from '@molstar/core/reactive';
import { type LogEntry, LogModel } from '../log/log-model.js';

// What a status line shows (status §1): the latest message (it fades; an error stays until dismissed), the running
// tasks, and chips (what the person works on). DOM-free (status §1): its timer is injected, its messages go into a
// `LogModel`, so a headless host posts, reads and saves the same.

export type StatusKind = 'info' | 'success' | 'warning' | 'error';

/** A message posted: the log's entry for it (no source: a message). */
/**
 * A message on show: its log entry, and its scope (status §1: what it is about, a plain string the app chooses, e.g.
 * a pane's id; the line of that scope shows it; none: the line that takes the rest).
 */
export type StatusMessage = LogEntry & { readonly kind: StatusKind; readonly scope?: string };

export interface StatusTask {
  readonly id: string;
  readonly label: string;
  /** 0..1, or undefined when indeterminate. */
  readonly progress?: number;
  /** It holds the whole app while it runs (an export): the app shows it over everything, not in the status line. */
  readonly overlay?: boolean;
  abort?(): void;
}

export interface StatusChip {
  readonly id: string;
  readonly label: string;
  /** Clicking the chip (focus the selection). */
  activate?(): void;
  /** The chip's ✕ (clear the selection). */
  dismiss?(): void;
}

/** The fade's timer (the globals by default; a test or a host gives its own). */
export interface StatusTimers {
  setTimeout(fn: () => void, ms: number): unknown;
  clearTimeout(id: unknown): void;
}

export interface StatusModelOptions {
  /** How long a message stays on show (4000 ms); errors stay until dismissed. */
  readonly fadeMs?: number;
  /** The log its messages are written into (one is made when absent). */
  readonly log?: LogModel;
  readonly timers?: StatusTimers;
}

const GlobalTimers: StatusTimers = {
  setTimeout: (fn, ms) => setTimeout(fn, ms),
  clearTimeout: (id) => clearTimeout(id as ReturnType<typeof setTimeout>),
};

/**
 * The status (status §1): `post` a message (on show, said, in the log), `announce` (said only), `dismiss`; tasks from
 * two writers (`putTask`/`removeTask` one by one, `setTasks` mirroring a source whole, never dropping put ones), kept
 * by id (`taskIds`, `task(id)`) so a progress event tells one pill; `chips`.
 */
export class StatusModel extends ReactiveModel {
  readonly fadeMs: number;
  /** The log its messages go into (the status line's history; a Log view shows it). */
  readonly log: LogModel;
  readonly #timers: StatusTimers;
  readonly #ownsLog: boolean;
  readonly #current = new State<StatusMessage | undefined>(undefined);
  readonly #spoken = new State('');
  readonly #tasks = new State<readonly StatusTask[]>([]);
  readonly #lineTasks = new KeyedState<string, StatusTask>();
  readonly #taskIds = new State<readonly string[]>([], sameIds);
  /** The tasks put one by one (`putTask`): `setTasks` leaves them. */
  readonly #put = new Set<string>();
  readonly #chips = new State<readonly StatusChip[]>([]);
  #fade: unknown;

  constructor(options: StatusModelOptions = {}) {
    super();
    this.fadeMs = options.fadeMs ?? 4000;
    this.#ownsLog = !options.log;
    this.log = options.log ?? new LogModel();
    this.#timers = options.timers ?? GlobalTimers;
  }

  /** The message on show, if any. */
  get current(): ReadonlyState<StatusMessage | undefined> {
    return this.#current.asReadonly;
  }

  /**
   * What the one polite live region says: the latest message posted or announced. A repeat of the same words differs
   * by a trailing no-break space, so it is said again.
   */
  get spoken(): ReadonlyState<string> {
    return this.#spoken.asReadonly;
  }

  /** Every task, overlay ones too. */
  get tasks(): ReadonlyState<readonly StatusTask[]> {
    return this.#tasks.asReadonly;
  }

  /** The ids of the tasks the status line shows (all but `overlay`), in order; told when they change, not on progress. */
  get taskIds(): ReadonlyState<readonly string[]> {
    return this.#taskIds.asReadonly;
  }

  /** A status-line task by id: told when its record changes (its progress), undefined once it is gone. */
  task(id: string): ReadonlyState<StatusTask | undefined> {
    return this.#lineTasks.key(id);
  }

  get chips(): ReadonlyState<readonly StatusChip[]> {
    return this.#chips.asReadonly;
  }

  /** Says `text` to screen readers only: nothing on show, nothing in the log (a selection's result). */
  announce(text: string): void {
    if (!text) return;
    this.#spoken.set(text === this.#spoken.value ? `${text} ` : text);
  }

  /**
   * Shows a message (it fades after `fadeMs` unless an error), says it, and writes it into the log. `scope`: what it
   * is about (status §1: the line of that scope shows it; the app maps scopes to its parts, a pane each).
   */
  post(
    text: string,
    kind: StatusKind = 'info',
    options: { readonly detail?: string; readonly scope?: string } = {},
  ): StatusMessage {
    const entry = this.log.write({ text, kind, detail: options.detail }) as StatusMessage;
    const message: StatusMessage = options.scope === undefined ? entry : { ...entry, scope: options.scope };
    batch(() => {
      this.#current.set(message);
      this.announce(text);
    });
    this.#timers.clearTimeout(this.#fade);
    this.#fade = kind === 'error' ? undefined : this.#timers.setTimeout(() => this.dismiss(message.id), this.fadeMs);
    return message;
  }

  /** Hides the message on show (only that one, if `id` is given). */
  dismiss(id?: number): void {
    if (id !== undefined && this.#current.value?.id !== id) return;
    this.#timers.clearTimeout(this.#fade);
    this.#fade = undefined;
    this.#current.set(undefined);
  }

  /** Adds or replaces a task by id; it stays until `removeTask`, whatever `setTasks` sets. */
  putTask(task: StatusTask): void {
    this.#put.add(task.id);
    const tasks = this.#tasks.value;
    const i = tasks.findIndex((t) => t.id === task.id);
    this.#setTasks(i < 0 ? [...tasks, task] : tasks.map((t, j) => (j === i ? task : t)));
  }

  removeTask(id: string): void {
    this.#put.delete(id);
    if (this.#tasks.value.some((t) => t.id === id)) this.#setTasks(this.#tasks.value.filter((t) => t.id !== id));
  }

  /** Sets the tasks one source mirrors: every task but those put one by one, which stay until removed. */
  setTasks(tasks: readonly StatusTask[]): void {
    if (!this.#put.size) {
      this.#setTasks(tasks);
      return;
    }
    const given = new Set<string>();
    for (const t of tasks) given.add(t.id);
    const put = this.#tasks.value.filter((t) => this.#put.has(t.id) && !given.has(t.id));
    this.#setTasks(put.length ? [...tasks, ...put] : tasks);
  }

  setChips(chips: readonly StatusChip[]): void {
    this.#chips.set(chips);
  }

  /** The tasks, and the line's by id: the records first, so a pill made for a new id reads its task. */
  #setTasks(tasks: readonly StatusTask[]): void {
    const byId = new Map<string, StatusTask>();
    const ids: string[] = [];
    for (const t of tasks) {
      if (t.overlay || byId.has(t.id)) continue;
      byId.set(t.id, t);
      ids.push(t.id);
    }
    batch(() => {
      this.#lineTasks.replace(byId);
      this.#taskIds.set(ids);
      this.#tasks.set(tasks);
    });
  }

  protected override onRetire(): void {
    this.#timers.clearTimeout(this.#fade);
    this.#fade = undefined;
    if (this.#ownsLog) this.log.dispose();
  }
}

function sameIds(a: readonly string[], b: readonly string[]): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}
