/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Observable, Subject } from 'rxjs';
import { keyedState } from '../reactive/keyed-state.js';
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

/** What the Log shows: everything, the messages the status line showed, one log-only source, or the problems. */
export type LogFilter = 'all' | 'messages' | 'problems' | `source:${string}`;

export interface LogOptions {
  readonly detail?: string;
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

/** A warning or an error: what the Log's Problems filter shows, and what its bound drops last. */
function isProblem(m: StatusMessage): boolean {
  return m.kind === 'warning' || m.kind === 'error';
}

/**
 * What the status line on the canvas shows: the latest message (it fades after a while; errors stay until dismissed),
 * running tasks, and chips (selection, focus). Keeps a bounded history for the log overlay. The viewer's
 * `managers.notifications` (api.md §11, prism-viewer README decision 22: `current`, `tasks`, `chips`; `posted`).
 */
export class NotificationsManager extends ReactiveModel {
  readonly #current = observableState<StatusMessage | undefined>(undefined);
  readonly #history = observableState<readonly StatusMessage[]>([]);
  readonly #log = observableState<readonly StatusMessage[]>([]);
  readonly #tasks = observableState<readonly StatusTask[]>([]);
  /** The status line's tasks (not `overlay`) by id, and their ids: a pill follows its task (perf audit 2, 11/B3). */
  readonly #lineTasks = keyedState<string, StatusTask>();
  readonly #lineTaskIds = observableState<readonly string[]>([]);
  /** The tasks put one by one (`putTask`), by id: `setTasks` leaves them. */
  readonly #put = new Set<string>();
  readonly #chips = observableState<readonly StatusChip[]>([]);
  readonly #posted = new Subject<StatusMessage>();
  readonly #spoken = observableState('');
  readonly #logFilter = observableState<LogFilter>('all');
  #serial = 0;
  #fade: ReturnType<typeof setTimeout> | undefined;
  /** Log-only entries waiting for the end of the task: a burst (a load's hundreds of transforms) is one change. */
  #pending: StatusMessage[] = [];
  /** The log-only sources written, sorted: the Log's filters (told only when one is new). */
  readonly #sources = observableState<readonly string[]>([]);

  constructor(
    readonly options: {
      fadeMs?: number;
      /** Messages kept in `history` (200). */
      historySize?: number;
      /** Entries kept in the Log (`logged`), of every source together (300). */
      logSize?: number;
    } = {},
  ) {
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

  /**
   * What the Log shows, newest first: the messages (`history`) and what was written to the log only (`log`: the
   * plugin's transforms, for one), within one bound (`logSize`) where warnings and errors go last (`#bounded`).
   */
  get logged(): ReadonlyState<readonly StatusMessage[]> {
    return this.#log.asReadonly;
  }

  /** The log-only sources written so far, sorted (the Log's filter per source); told only when one is new. */
  get logSources(): ReadonlyState<readonly string[]> {
    return this.#sources.asReadonly;
  }

  /** The Log's filter, one for every Log view: the status line sets it when it opens the Log (an error: the problems). */
  get logFilter(): ReadonlyState<LogFilter> {
    return this.#logFilter.asReadonly;
  }

  setLogFilter(filter: LogFilter): void {
    this.#logFilter.set(filter);
  }

  get tasks(): ReadonlyState<readonly StatusTask[]> {
    return this.#tasks.asReadonly;
  }

  /**
   * The ids of the tasks the status line shows (all but `overlay` ones), in order; told only when they change, not on a
   * task's progress. Each task by `lineTask(id)`.
   */
  get lineTaskIds(): ReadonlyState<readonly string[]> {
    return this.#lineTaskIds.asReadonly;
  }

  /** A status-line task by id: told when that task's record changes (its progress), undefined once it is gone. */
  lineTask(id: string): ReadonlyState<StatusTask | undefined> {
    return this.#lineTasks.key(id);
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
    this.#history.set([message, ...this.#history.value].slice(0, this.options.historySize ?? 200));
    this.#log.set(this.#bounded([message, ...this.#log.value]));
    this.#current.set(message);
    this.announce(text);
    clearTimeout(this.#fade);
    if (kind !== 'error') this.#fade = setTimeout(() => this.dismiss(message.id), this.options.fadeMs ?? 4000);
    this.#posted.next(message);
  }

  /**
   * Writes `text` to the log only (`logged`, the Log tab, filtered by `source`): never on show, never said, not in
   * `history`. Written once a task, so a burst is one change.
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
    this.#pending.push(message);
  }

  #flush() {
    if (!this.#pending.length) return;
    const pending = this.#pending;
    this.#pending = [];
    const next: StatusMessage[] = [];
    let sources: string[] | undefined;
    for (let i = pending.length - 1; i >= 0; i--) {
      const m = pending[i];
      next.push(m);
      if (!(sources ?? this.#sources.value).includes(m.source!)) (sources ??= [...this.#sources.value]).push(m.source!);
    }
    for (const m of this.#log.value) next.push(m);
    this.#log.set(this.#bounded(next));
    if (sources) this.#sources.set(sources.sort());
  }

  /**
   * The log within its one bound, `logSize` (perf audit 2, 11/B5): past it, the oldest entries go, warnings and errors
   * last, so a burst (a load's transforms, a run of plugin warnings) never pushes a problem out while anything else is
   * left to drop. Only when the problems alone fill the bound do the oldest problems go.
   */
  #bounded(list: StatusMessage[]): StatusMessage[] {
    const size = this.options.logSize ?? 300;
    if (list.length <= size) return list;
    let problems = 0;
    for (const m of list) if (isProblem(m)) problems++;
    // Newest first: the first `left` of each kind are kept.
    if (problems >= size) {
      let left = size;
      return list.filter((m) => isProblem(m) && left-- > 0);
    }
    let left = size - problems;
    return list.filter((m) => isProblem(m) || left-- > 0);
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
    this.#setTasks(put.length ? [...tasks, ...put] : tasks);
  }

  /** The tasks, and the status line's by id: the records first, so a pill made for a new id reads its task. */
  #setTasks(tasks: readonly StatusTask[]): void {
    const byId = new Map<string, StatusTask>();
    const ids: string[] = [];
    for (const t of tasks) {
      if (t.overlay || byId.has(t.id)) continue;
      byId.set(t.id, t);
      ids.push(t.id);
    }
    this.#lineTasks.replace(byId);
    // Compared item by item: the same ids tell nothing.
    this.#lineTaskIds.set(ids);
    this.#tasks.set(tasks);
  }

  setChips(chips: readonly StatusChip[]): void {
    this.#chips.set(chips);
  }

  protected onRetire(): void {
    clearTimeout(this.#fade);
    this.#posted.complete();
  }
}
