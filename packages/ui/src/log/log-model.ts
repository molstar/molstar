/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { batch, derive, type ReadonlyState, State } from '@molstar/core/reactive';

// A generic log (log §1): what an app, a view or an agent wrote, oldest first, within one bound. DOM-free
// (headless): a script, a test or a headless host writes and reads it with no document; `LogView` shows it.

export type LogKind = 'info' | 'success' | 'warning' | 'error' | 'debug';

export interface LogEntry {
  readonly id: number;
  readonly kind: LogKind;
  readonly text: string;
  /** When it was written (ms since the epoch, the model's `now`). */
  readonly time: number;
  /** What wrote it (a filter each); absent: a message (the status line's). */
  readonly source?: string;
  /** Said beside the text (how long it took: "12 ms"). */
  readonly detail?: string;
}

/** What `write` takes: the entry without its id (and its time, unless given). */
export interface LogInput {
  readonly text: string;
  readonly kind?: LogKind;
  readonly source?: string;
  readonly detail?: string;
  readonly time?: number;
}

/**
 * What the log shows (log §1): the kinds (any of them; none: every kind) and the source (a source's entries, `null`
 * for the messages: the entries without one; absent: every source). The two combine.
 */
export interface LogFilter {
  readonly kinds?: readonly LogKind[];
  readonly source?: string | null;
}

/** Every entry. */
export const LogAll: LogFilter = {};
/** The warnings and errors (a status line opens the log at these from an error). */
export const LogProblems: LogFilter = { kinds: ['warning', 'error'] };

/** The kinds a filter shows, in a fixed order (none: every kind). */
export const LogKinds: readonly LogKind[] = ['info', 'success', 'warning', 'error', 'debug'];

/** Whether a filter shows everything (no kind, no source). */
export function isLogAll(filter: LogFilter): boolean {
  return !filter.kinds?.length && filter.source === undefined;
}

export interface LogModelOptions {
  /** Entries kept, of every source together (300); warnings and errors are dropped last (log §2). */
  readonly size?: number;
  /** The clock for an entry's `time` (ms since the epoch). */
  readonly now?: () => number;
  /** Runs the batched write (a microtask by default; a test or a host may run it otherwise). */
  readonly schedule?: (fn: () => void) => void;
}

/** A warning or an error: what `LogProblems` shows and what the bound drops last. */
export function isProblem(e: Pick<LogEntry, 'kind'>): boolean {
  return e.kind === 'warning' || e.kind === 'error';
}

/** Whether an entry passes a filter and a query (lower case, trimmed; empty: any). */
export function passesLog(e: LogEntry, filter: LogFilter, query: string): boolean {
  if (filter.kinds?.length && !filter.kinds.includes(e.kind)) return false;
  if (filter.source !== undefined && (e.source ?? null) !== filter.source) return false;
  if (!query) return true;
  return (
    e.text.toLowerCase().includes(query) ||
    !!e.detail?.toLowerCase().includes(query) ||
    !!e.source?.toLowerCase().includes(query)
  );
}

/**
 * The list within `size`, oldest first (log §2): past it, the oldest entries go, warnings and errors last, so a burst
 * (a load's hundreds of entries) never pushes a problem out while anything else is left to drop. Only when the
 * problems alone fill the bound do the oldest problems go. One pass, one array.
 */
export function boundLog(list: readonly LogEntry[], size: number): readonly LogEntry[] {
  if (list.length <= size) return list;
  let problems = 0;
  for (const e of list) if (isProblem(e)) problems++;
  const out: LogEntry[] = [];
  if (problems >= size) {
    let skip = problems - size;
    for (const e of list) if (isProblem(e) && skip-- <= 0) out.push(e);
    return out;
  }
  let skip = list.length - size;
  for (const e of list) {
    if (!isProblem(e) && skip > 0) {
      skip--;
      continue;
    }
    out.push(e);
  }
  return out;
}

/**
 * A log (log §1): `write` an entry (a burst within a task is one change: written a microtask later, or at `flush`),
 * one bound (`size`) with the problems kept longest, the sources as written (forgotten with their last entry), a
 * filter and a text query, the entries
 * they pass (`filtered`), `clear`. `save()` gives the entries as data (headless).
 */
export class LogModel {
  readonly size: number;
  readonly #now: () => number;
  readonly #schedule: (fn: () => void) => void;
  readonly #entries = new State<readonly LogEntry[]>([]);
  readonly #sources = new State<readonly string[]>([]);
  readonly #filter = new State<LogFilter>(LogAll);
  readonly #query = new State('');
  readonly #known = new Set<string>();
  #pending: LogEntry[] = [];
  #serial = 0;
  #disposed = false;
  readonly filtered: ReadonlyState<readonly LogEntry[]>;

  constructor(options: LogModelOptions = {}) {
    this.size = Math.max(1, options.size ?? 300);
    this.#now = options.now ?? Date.now;
    // Called as a function, not a method: a browser's `queueMicrotask` throws when called on another `this`.
    this.#schedule = options.schedule ?? ((fn) => queueMicrotask(fn));
    this.filtered = derive([this.#entries, this.#filter, this.#query], (list, filter, q): readonly LogEntry[] => {
      const query = q.trim().toLowerCase();
      return isLogAll(filter) && !query ? list : list.filter((e) => passesLog(e, filter, query));
    });
  }

  /** The entries, oldest first (newest last), within `size`. */
  get entries(): ReadonlyState<readonly LogEntry[]> {
    return this.#entries.asReadonly;
  }

  /**
   * The sources written, in the order first written (a filter each); told only when one is new, or when the bound
   * dropped a source's last entry (it is forgotten).
   */
  get sources(): ReadonlyState<readonly string[]> {
    return this.#sources.asReadonly;
  }

  get filter(): ReadonlyState<LogFilter> {
    return this.#filter.asReadonly;
  }

  /** The text filter as typed. */
  get query(): ReadonlyState<string> {
    return this.#query.asReadonly;
  }

  /** Shows what `filter` names (`LogAll`, `LogProblems`, `{ kinds, source }`). */
  setFilter(filter: LogFilter): void {
    this.#filter.set(filter);
  }

  /**
   * A kind's toggle (log §3): from every kind, shows that kind alone; otherwise adds it or takes it out (the last one out:
   * every kind again). The source stays.
   */
  toggleKind(kind: LogKind): void {
    const { kinds = [], source } = this.#filter.value;
    const next = kinds.includes(kind)
      ? kinds.filter((k) => k !== kind)
      : LogKinds.filter((k) => k === kind || kinds.includes(k));
    this.#filter.set({ ...(next.length ? { kinds: next } : {}), ...(source !== undefined ? { source } : {}) });
  }

  /** Shows one source's entries (`null`: the messages; `undefined`: every source). The kinds stay. */
  setSource(source: string | null | undefined): void {
    const { kinds } = this.#filter.value;
    this.#filter.set({ ...(kinds?.length ? { kinds } : {}), ...(source !== undefined ? { source } : {}) });
  }

  setQuery(query: string): void {
    this.#query.set(query);
  }

  /**
   * Writes an entry: given its id and time now, in `entries` at the end of the task (with the rest of a burst, one
   * change). Returns the entry.
   */
  write(input: LogInput): LogEntry {
    const entry: LogEntry = {
      id: ++this.#serial,
      kind: input.kind ?? 'info',
      text: input.text,
      time: input.time ?? this.#now(),
      ...(input.source !== undefined ? { source: input.source } : {}),
      ...(input.detail ? { detail: input.detail } : {}),
    };
    if (this.#disposed) return entry;
    if (!this.#pending.length) this.#schedule(() => this.flush());
    this.#pending.push(entry);
    return entry;
  }

  /** Writes what waits now (a reader that cannot wait for the task's end). */
  flush(): void {
    const pending = this.#pending;
    if (!pending.length) return;
    this.#pending = [];
    let sources: string[] | undefined;
    for (const e of pending) {
      if (e.source === undefined || this.#known.has(e.source)) continue;
      this.#known.add(e.source);
      (sources ??= [...this.#sources.value]).push(e.source);
    }
    const list = this.#entries.value;
    const next = new Array<LogEntry>(list.length + pending.length);
    for (let i = 0; i < list.length; i++) next[i] = list[i];
    for (let i = 0; i < pending.length; i++) next[list.length + i] = pending[i];
    const bounded = boundLog(next, this.size);
    // A source whose entries the bound dropped is forgotten (log §2): ids in source names would grow the list (and the
    // filters) without limit. One pass over the bounded list, only when the bound dropped something.
    let forgotten = false;
    if (bounded.length < next.length && this.#known.size) {
      const present = new Set<string>();
      for (const e of bounded) if (e.source !== undefined) present.add(e.source);
      if (present.size < this.#known.size) {
        for (const s of this.#known) if (!present.has(s)) this.#known.delete(s);
        sources = (sources ?? this.#sources.value).filter((s) => present.has(s));
        forgotten = true;
      }
    }
    batch(() => {
      this.#entries.set(bounded);
      if (sources) this.#sources.set(sources);
      // The filter's source is gone: every source again (the filter stays valid).
      const source = this.#filter.value.source;
      if (forgotten && typeof source === 'string' && !this.#known.has(source)) this.setSource(undefined);
    });
  }

  /** Forgets every entry and source (the filter stays, unless it named a source: then every source). */
  clear(): void {
    this.#pending = [];
    this.#known.clear();
    batch(() => {
      this.#entries.set([]);
      this.#sources.set([]);
      if (typeof this.#filter.value.source === 'string') this.setSource(undefined);
    });
  }

  /** The entries as data (what waits written first): a headless host's log. */
  save(): readonly LogEntry[] {
    this.flush();
    return this.#entries.value;
  }

  /** Writes nothing more. */
  dispose(): void {
    this.#disposed = true;
    this.#pending = [];
  }
}
