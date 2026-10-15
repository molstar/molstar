/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Observable, Subject } from 'rxjs';
import { observableState } from '../reactive/observable-state.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

// One undo for everything (spec sequence.md §7.2–7.4, README §3.3): every undoable step, a person's call or an agent's
// whole turn, is one entry with a checkpoint taken before it; Undo restores the latest. Generic and DOM-free: the app
// gives the checkpoints through a `HistoryPort` (the viewer's: the plugin's state and the sequence's).
//
// - A step is one undoable run, or one group. A group is atomic (`group`: a script's or a command's several calls, all
//   or nothing) or open (`beginGroup`: an agent's turn, whose calls join one entry while it is the latest).
// - Undoable runs are serialized: one at a time, each with its checkpoint, so a failed run can be restored without
//   touching another's work.
// - A failed run restores its checkpoint; a failed atomic group restores the group's.
// - A change outside the history (reported by the app with `invalidate`) makes every earlier step unsafe.
// - No redo yet: entries hold only the checkpoint before them.

/** Who ran it: `person`, `agent` (with its turn), `script`, `url`, `link`; the app's kinds. */
export interface HistoryOrigin {
  readonly kind: string;
  /** An agent's turn: its calls join one entry. */
  readonly turn?: string;
}

export interface HistoryPort<C> {
  /** Synchronous and cheap; measured against the budget (§7.4). */
  checkpoint(): C;
  restore(checkpoint: C): Promise<void>;
  /** The checkpoint's size in bytes, without data payloads; absent: not measured. */
  size?(checkpoint: C): number;
}

export interface HistoryStep {
  readonly id: number;
  /** What it did, for History and the Undo button ("Undo Change background"). */
  readonly label: string;
  readonly origin: HistoryOrigin;
  /** When it started (ms since the epoch). */
  readonly time: number;
  /** The runs inside a group, in order (their labels); a single run lists itself. */
  readonly runs: readonly string[];
  /** Why Undo cannot restore it (a change outside the API after it, a large scene); undefined: it can. */
  readonly unsafe?: string;
}

export type UndoResult = { ok: true; step: HistoryStep } | { ok: false; reason: string; laterSteps?: number };

export interface HistoryOptions {
  /** Most entries kept; the oldest go first. */
  readonly capacity?: number;
  /** Over either, per-step checkpoints turn off for the scene (sequence.md §7.4; a first guess). */
  readonly budget?: { readonly ms: number; readonly bytes: number };
  /** Over either, even a turn's one checkpoint is not taken. */
  readonly hardLimit?: { readonly ms: number; readonly bytes: number };
  /** The clock (tests). */
  readonly now?: () => number;
}

/** What `run` and `group` give back: the work's value, and whether a failure restored the checkpoint. */
export interface StepOutcome<T> {
  readonly value: T;
  readonly restored: boolean;
  /** The step recorded; undefined when nothing was recorded (a failure, or `record: false`). */
  readonly step?: HistoryStep;
}

export interface RunOptions<T> {
  /** Whether the value is a failure (restore the checkpoint, record nothing). Thrown errors are failures too. */
  readonly failed?: (value: T) => boolean;
  /** Not undoable (an upload): runs in order with the others, recorded as a step that cannot be undone. */
  readonly undoable?: boolean;
}

interface Entry<C> {
  step: HistoryStep;
  checkpoint: C | undefined;
  /** An open group's key while the group can still take runs. */
  openKey?: string;
}

interface OpenGroup {
  readonly key: string;
  readonly label: string;
  readonly origin: HistoryOrigin;
}

interface AtomicGroup<C> {
  readonly token: object;
  readonly label: string;
  readonly origin: HistoryOrigin;
  readonly checkpoint: C | undefined;
  readonly runs: string[];
}

/** The key that groups an origin's calls: an agent's turn. */
export function historyGroupKey(origin: HistoryOrigin): string | undefined {
  return origin.turn ? `${origin.kind}:${origin.turn}` : undefined;
}

export class HistoryManager<C = unknown> {
  readonly #port: HistoryPort<C>;
  readonly #options: Required<Omit<HistoryOptions, 'now'>> & { now: () => number };
  #entries: Entry<C>[] = [];
  readonly #steps = observableState<readonly HistoryStep[]>([]);
  readonly #running = observableState(0);
  readonly #restoring = observableState(false);
  readonly #checkpointing = observableState(true);
  readonly #undoBlocked = observableState<string | undefined>('Nothing to undo');
  readonly #undone = new Subject<HistoryStep>();
  readonly #checkpointingOff = new Subject<{ ms: number; bytes: number }>();
  readonly #open = new Map<string, OpenGroup>();
  #atomic: AtomicGroup<C> | undefined;
  #serial = 0;
  #lock: Promise<void> = Promise.resolve();

  constructor(port: HistoryPort<C>, options: HistoryOptions = {}) {
    this.#port = port;
    this.#options = {
      capacity: options.capacity ?? 100,
      budget: options.budget ?? { ms: 50, bytes: 5_000_000 },
      hardLimit: options.hardLimit ?? { ms: 500, bytes: 50_000_000 },
      now: options.now ?? (() => Date.now()),
    };
  }

  /** The undoable steps, newest first. */
  get steps(): ReadonlyState<readonly HistoryStep[]> {
    return this.#steps.asReadonly;
  }

  /** Why Undo cannot run now (a run in flight, nothing to undo, an unsafe step); undefined: it can. */
  get undoBlocked(): ReadonlyState<string | undefined> {
    return this.#undoBlocked.asReadonly;
  }

  /** Whether Undo can run now. */
  get canUndo(): boolean {
    return this.#undoBlocked.value === undefined;
  }

  /** Undoable runs in flight. */
  get running(): ReadonlyState<number> {
    return this.#running.asReadonly;
  }

  /** Per-step checkpoints are on; off for a large scene (§7.4). */
  get checkpointing(): ReadonlyState<boolean> {
    return this.#checkpointing.asReadonly;
  }

  /** A step undone. */
  get undone(): Observable<HistoryStep> {
    return this.#undone.asObservable();
  }

  /** Checkpoints just turned off for a large scene, with the measurement that did it (say so once). */
  get checkpointingOff(): Observable<{ ms: number; bytes: number }> {
    return this.#checkpointingOff.asObservable();
  }

  /**
   * Whether the history is changing the state itself (a run or a restore): the app's detector of changes outside the
   * API ignores what happens meanwhile.
   */
  get applying(): boolean {
    return this.#running.value > 0 || this.#restoring.value || this.#atomic !== undefined;
  }

  /** Turns per-step checkpoints on or off (Settings: "Undo for large scenes"); on again after a clear. */
  setCheckpointing(on: boolean): void {
    this.#checkpointing.set(on);
  }

  /**
   * Runs one undoable piece of work as a step: in order with the others, a checkpoint before it, recorded when it
   * succeeds, its checkpoint restored when it fails. Inside an atomic group it is part of the group; an open group of
   * its origin (an agent's turn) gets it as one more run.
   */
  async run<T>(label: string, origin: HistoryOrigin, work: () => Promise<T>, options: RunOptions<T> = {}) {
    const atomic = this.#atomic;
    if (atomic) return this.#runInAtomic(atomic, label, work, options);
    // In flight from the moment it is asked for: Undo waits for queued runs too.
    this.#running.set(this.#running.value + 1);
    this.#updateBlocked();
    const release = await this.#acquire();
    try {
      const key = historyGroupKey(origin);
      const group = key ? this.#open.get(key) : undefined;
      const top = this.#entries[0];
      // An entry takes more runs only while it is its group's and the latest; anything else closes it.
      if (top?.openKey && top.openKey !== key) top.openKey = undefined;
      const joining = !!group && top?.openKey === key;
      const checkpoint = this.#take(!!group);
      let value: T;
      try {
        value = await work();
      } catch (error) {
        await this.#restore(checkpoint.value);
        throw error;
      }
      if (options.failed?.(value)) {
        const restored = await this.#restore(checkpoint.value);
        return { value, restored } satisfies StepOutcome<T>;
      }
      if (joining && top) {
        top.step = { ...top.step, runs: [...top.step.runs, label] };
        this.#publish();
        return { value, restored: false, step: top.step } satisfies StepOutcome<T>;
      }
      const step = this.#push(
        group ? group.label : label,
        origin,
        [label],
        options.undoable === false ? undefined : checkpoint.value,
        options.undoable === false ? 'It cannot be undone' : checkpoint.unsafe,
        group?.key,
      );
      return { value, restored: false, step } satisfies StepOutcome<T>;
    } finally {
      this.#running.set(this.#running.value - 1);
      release();
      this.#updateBlocked();
    }
  }

  /**
   * Runs `work` as one atomic step: one checkpoint, every run inside it part of the step, and a failure (a throw, or a
   * value `failed` says is one) restores the checkpoint, so it is all or nothing.
   */
  async group<T>(label: string, origin: HistoryOrigin, work: () => Promise<T>, options: RunOptions<T> = {}) {
    if (this.#atomic) {
      // Nested: part of the outer group.
      return this.#runInAtomic(this.#atomic, label, work, options);
    }
    const release = await this.#acquire();
    const checkpoint = this.#take(false);
    const group: AtomicGroup<C> = { token: {}, label, origin, checkpoint: checkpoint.value, runs: [] };
    this.#atomic = group;
    this.#updateBlocked();
    try {
      let value: T;
      try {
        value = await work();
      } catch (error) {
        this.#atomic = undefined;
        await this.#restore(checkpoint.value);
        throw error;
      }
      this.#atomic = undefined;
      if (options.failed?.(value)) {
        const restored = await this.#restore(checkpoint.value);
        return { value, restored } satisfies StepOutcome<T>;
      }
      if (group.runs.length === 0) return { value, restored: false } satisfies StepOutcome<T>;
      const step = this.#push(label, origin, group.runs, checkpoint.value, checkpoint.unsafe);
      return { value, restored: false, step } satisfies StepOutcome<T>;
    } finally {
      this.#atomic = undefined;
      release();
      this.#updateBlocked();
    }
  }

  /**
   * Opens a group for an origin's calls (an agent's turn): they join one entry while it is the latest; a call of
   * another origin closes the entry, and the group's next call opens a new one under the same label. Returns `end`.
   */
  beginGroup(label: string, origin: HistoryOrigin): () => void {
    const key = historyGroupKey(origin);
    if (!key) throw new Error('beginGroup needs an origin with a turn');
    this.#open.set(key, { key, label, origin });
    return () => this.endGroup(origin);
  }

  endGroup(origin: HistoryOrigin): void {
    const key = historyGroupKey(origin);
    if (!key) return;
    this.#open.delete(key);
    for (const e of this.#entries) if (e.openKey === key) e.openKey = undefined;
  }

  /** Restores the latest step; blocked while a run is in flight, or when the step is unsafe. */
  async undo(): Promise<UndoResult> {
    const blocked = this.#undoBlocked.value;
    if (blocked) return { ok: false, reason: blocked };
    const release = await this.#acquire();
    try {
      const top = this.#entries[0];
      if (!top) return { ok: false, reason: 'Nothing to undo' };
      if (top.step.unsafe || top.checkpoint === undefined)
        return { ok: false, reason: top.step.unsafe ?? 'No checkpoint' };
      await this.#restore(top.checkpoint);
      this.#entries = this.#entries.slice(1);
      this.#publish();
      this.#undone.next(top.step);
      return { ok: true, step: top.step };
    } finally {
      release();
      this.#updateBlocked();
    }
  }

  /**
   * Undoes a turn's entries. While they are the latest, it restores the checkpoint before the turn's first one; with
   * later steps above them it asks first (`laterSteps`), and with `discardLater` restores past those too.
   */
  async undoTurn(turn: string, options: { discardLater?: boolean } = {}): Promise<UndoResult> {
    const blocked = this.#undoBlocked.value;
    if (blocked && blocked !== 'Nothing to undo') return { ok: false, reason: blocked };
    const indices = this.#entries.flatMap((e, i) => (e.step.origin.turn === turn ? [i] : []));
    if (!indices.length) return { ok: false, reason: 'That turn changed nothing that can be undone' };
    const last = indices[indices.length - 1];
    const later = this.#entries.slice(0, last + 1).filter((e) => e.step.origin.turn !== turn).length;
    if (later > 0 && !options.discardLater) {
      return {
        ok: false,
        reason: `${later} later step${later === 1 ? '' : 's'} would be undone too`,
        laterSteps: later,
      };
    }
    const release = await this.#acquire();
    try {
      const target = this.#entries[last];
      const unsafe = this.#entries.slice(0, last + 1).find((e) => e.step.unsafe || e.checkpoint === undefined);
      if (unsafe) return { ok: false, reason: unsafe.step.unsafe ?? 'No checkpoint' };
      await this.#restore(target.checkpoint);
      const undone = this.#entries.slice(0, last + 1);
      this.#entries = this.#entries.slice(last + 1);
      this.#publish();
      for (const e of undone) this.#undone.next(e.step);
      return { ok: true, step: target.step };
    } finally {
      release();
      this.#updateBlocked();
    }
  }

  /**
   * A change outside the history: every step taken before it can no longer be restored safely (its checkpoint would
   * drop the change, or fail on it). History says so on them.
   */
  invalidate(reason: string): void {
    let changed = false;
    for (const e of this.#entries) {
      if (e.step.unsafe) continue;
      e.step = { ...e.step, unsafe: reason };
      changed = true;
    }
    for (const e of this.#entries) e.openKey = undefined;
    if (changed) this.#publish();
    this.#updateBlocked();
  }

  /** Forgets every step (a cleared scene); checkpoints turn back on. */
  clear(): void {
    this.#entries = [];
    this.#checkpointing.set(true);
    this.#publish();
    this.#updateBlocked();
  }

  dispose(): void {
    this.#entries = [];
    this.#open.clear();
    this.#undone.complete();
    this.#checkpointingOff.complete();
  }

  async #runInAtomic<T>(
    group: AtomicGroup<C>,
    label: string,
    work: () => Promise<T>,
    options: RunOptions<T>,
  ): Promise<StepOutcome<T>> {
    this.#running.set(this.#running.value + 1);
    try {
      const value = await work();
      if (!options.failed?.(value)) group.runs.push(label);
      return { value, restored: false };
    } finally {
      this.#running.set(this.#running.value - 1);
    }
  }

  /** Takes a checkpoint unless checkpoints are off (a group's is taken below the hard limit). */
  #take(group: boolean): { value: C | undefined; unsafe?: string } {
    if (!this.#checkpointing.value && !group) {
      return { value: undefined, unsafe: 'Undo is off for this large scene' };
    }
    const t0 = performance.now();
    const value = this.#port.checkpoint();
    const ms = performance.now() - t0;
    const bytes = this.#port.size?.(value) ?? 0;
    const { budget, hardLimit } = this.#options;
    if (ms > hardLimit.ms || bytes > hardLimit.bytes) {
      this.#turnOff(ms, bytes);
      return { value: undefined, unsafe: 'This scene is too large to undo' };
    }
    if (ms > budget.ms || bytes > budget.bytes) this.#turnOff(ms, bytes);
    return { value };
  }

  #turnOff(ms: number, bytes: number) {
    if (!this.#checkpointing.value) return;
    this.#checkpointing.set(false);
    this.#checkpointingOff.next({ ms, bytes });
  }

  async #restore(checkpoint: C | undefined): Promise<boolean> {
    if (checkpoint === undefined) return false;
    this.#restoring.set(true);
    try {
      await this.#port.restore(checkpoint);
      return true;
    } finally {
      this.#restoring.set(false);
    }
  }

  #push(
    label: string,
    origin: HistoryOrigin,
    runs: readonly string[],
    checkpoint: C | undefined,
    unsafe: string | undefined,
    openKey?: string,
  ): HistoryStep {
    const step: HistoryStep = {
      id: ++this.#serial,
      label,
      origin,
      time: this.#options.now(),
      runs,
      ...(unsafe ? { unsafe } : {}),
    };
    this.#entries = [{ step, checkpoint, openKey }, ...this.#entries].slice(0, this.#options.capacity);
    this.#publish();
    return step;
  }

  #publish() {
    this.#steps.set(this.#entries.map((e) => e.step));
    this.#updateBlocked();
  }

  #updateBlocked() {
    const top = this.#entries[0];
    this.#undoBlocked.set(
      this.#running.value > 0 || this.#atomic
        ? 'A change is running; stop it first'
        : !top
          ? 'Nothing to undo'
          : top.step.unsafe
            ? top.step.unsafe
            : undefined,
    );
  }

  #acquire(): Promise<() => void> {
    let release!: () => void;
    const next = new Promise<void>((resolve) => (release = resolve));
    const previous = this.#lock;
    this.#lock = previous.then(() => next);
    return previous.then(() => release);
  }
}
