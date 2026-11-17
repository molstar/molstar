/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { batch, type ReadonlyState, State } from '@molstar/core/reactive';

// One undo for everything (history §1): every undoable step, a person's call or an agent's whole turn, is one entry
// with a checkpoint taken before it; Undo restores the latest. Generic and DOM-free (history §1): the app gives the
// checkpoints through a `HistoryPort`; a view (`HistoryView`), a command or a script call the same methods.
//
// - A step is one undoable run, or one group. A group is atomic (`group`: several calls, all or nothing) or open
//   (`beginGroup`: an agent's turn, whose runs and groups join one entry while it is the latest).
// - Top-level runs and groups are serialized: one at a time, each with its checkpoint, so a failed one is restored
//   without touching another's work. Each gets a context (`HistoryContext`); work that calls back in with it is part of
//   that step: no lock, no step of its own. Calls without it queue.
// - A failed run restores its checkpoint. A group fails as a whole when its work throws, returns a failure, or any run
//   inside it failed (unless it says `continue`, and then each run inside has its own checkpoint).
// - A change outside the history (`invalidate`) makes every earlier step unsafe.
// - A run that changes nothing undoable (`undoable: false`) is no step; one that changed nothing (`record` false)
//   neither.
// - Redo: Undo checkpoints the state it leaves, and Redo restores it, until a new step clears the list.

/** Who ran it: `person`, `agent` (with its turn), `script`, `url`, `link`; the app's kinds. */
export interface HistoryOrigin {
  readonly kind: string;
  /** An agent's turn: its calls join one entry. */
  readonly turn?: string;
}

export interface HistoryPort<C> {
  /** Synchronous and cheap; measured against the budget (history §3). */
  checkpoint(): C;
  restore(checkpoint: C): Promise<void>;
  /** The checkpoint's size in bytes; absent: not measured. */
  size?(checkpoint: C): number;
}

export interface HistoryStep {
  readonly id: number;
  /** What it did ("Change background"): History's row and the Undo button's tooltip. */
  readonly label: string;
  readonly origin: HistoryOrigin;
  /** When it started (ms, the model's `now`). */
  readonly time: number;
  /** The runs inside a group, in order (their labels); a single run lists itself. */
  readonly runs: readonly string[];
  /** Why Undo cannot restore it (a change outside the history after it, a large scene); undefined: it can. */
  readonly unsafe?: string;
}

/** A run or group in progress: work that passes it back is part of that step. Opaque. */
export interface HistoryContext {
  readonly id: number;
}

export type UndoResult = { ok: true; step: HistoryStep } | { ok: false; reason: string; laterSteps?: number };

/** What Undo or Redo last did: a new object each time (a step undone again is told again). */
export interface HistoryMove {
  readonly kind: 'undo' | 'redo';
  /** The steps it took back or put back (an undone turn: several, newest first). */
  readonly steps: readonly HistoryStep[];
}

export interface HistoryOptions {
  /** Most entries listed (100); the oldest go first. */
  readonly capacity?: number;
  /** The latest entries that keep their checkpoint (10); older ones stay listed, unsafe ("too old"). */
  readonly checkpoints?: number;
  /** Over either, per-step checkpoints turn off (history §3). */
  readonly budget?: { readonly ms: number; readonly bytes: number };
  /** Over either, even a group's one checkpoint is not taken. */
  readonly hardLimit?: { readonly ms: number; readonly bytes: number };
  /** The clock of a step's `time` and of the hard limit's pause (ms; `Date.now`). */
  readonly now?: () => number;
  /** The clock that times a checkpoint (ms; `performance.now`). */
  readonly measure?: () => number;
}

/** What `run` and `group` give back: the work's value, and whether a failure restored the checkpoint. */
export interface StepOutcome<T> {
  readonly value: T;
  /** The work failed (a throw is rethrown instead). */
  readonly failed: boolean;
  /** A checkpoint was restored because of the failure. */
  readonly restored: boolean;
  /** The step recorded or joined; undefined when nothing was recorded. */
  readonly step?: HistoryStep;
}

export interface GroupOutcome<T> extends StepOutcome<T> {
  /** The runs inside that changed something, in order. */
  readonly runs: readonly string[];
  /** The runs inside that failed. */
  readonly failures: number;
}

export interface RunOptions<T> {
  /** Whether the value is a failure (restore the checkpoint, record nothing). Thrown errors are failures too. */
  readonly failed?: (value: T) => boolean;
  /** Whether it changed anything worth a step; false: "it is already so", no step. Default: it did. */
  readonly record?: (value: T) => boolean;
  /** Changes nothing the history could restore (an upload): runs in order, no checkpoint, no step. */
  readonly undoable?: boolean;
  /** Part of this run or group in progress: no lock, no step of its own. */
  readonly within?: HistoryContext;
}

export interface GroupOptions<T> {
  readonly failed?: (value: T) => boolean;
  /**
   * A run inside that fails: `rollback` (default) fails the whole group at its end; `continue` restores just that run
   * (each run inside then takes its own checkpoint) and the group goes on.
   */
  readonly onFailure?: 'rollback' | 'continue';
  readonly within?: HistoryContext;
}

interface Entry<C> {
  step: HistoryStep;
  checkpoint: C | undefined;
  /** An open group's key while the group can still take runs. */
  openKey?: string;
  /** The step's runs, appended in place while its group takes more: each published step reads its own length. */
  log: string[];
}

interface OpenGroup {
  readonly key: string;
  readonly label: string;
  readonly origin: HistoryOrigin;
}

interface Active {
  readonly id: number;
  readonly kind: 'run' | 'group';
  readonly parent?: Active;
  readonly onFailure: 'rollback' | 'continue';
  readonly runs: string[];
  failures: number;
}

/** The key that groups an origin's calls: an agent's turn. */
export function historyGroupKey(origin: HistoryOrigin): string | undefined {
  return origin.turn ? `${origin.kind}:${origin.turn}` : undefined;
}

const labelOf = <T>(label: string | ((value: T) => string), value: T) =>
  typeof label === 'function' ? label(value) : label;

type StepFields = Omit<HistoryStep, 'runs'>;

/** A step's fields but its runs (read without copying them). */
function fieldsOf({ id, label, origin, time, unsafe }: HistoryStep): StepFields {
  return unsafe === undefined ? { id, label, origin, time } : { id, label, origin, time, unsafe };
}

/**
 * A step whose `runs` are the first `length` of `log`, an append-only list: copied out the first time they are read (and
 * kept), so a step published per run costs nothing until a reader asks, and an older step never sees later runs (a
 * copy per run was R² for a long turn; history §2).
 */
function stepWith(fields: StepFields, log: readonly string[], length: number): HistoryStep {
  let runs: readonly string[] | undefined;
  return {
    ...fields,
    get runs() {
      runs ??= length === 0 ? [] : log.slice(0, length);
      return runs;
    },
  };
}

/** Checkpoints over the time budget in a row that turn per-step checkpoints off (history §3). */
const SlowCheckpointsToTurnOff = 3;

/** After a checkpoint over the hard limit, how long none is even taken (history §3): each try is a snapshot thrown away. */
export const OverHardLimitPauseMs = 30_000;

const Running = 'A change is running; stop it first';
const NothingToUndo = 'Nothing to undo';
const NothingToRedo = 'Nothing to redo';

/**
 * Undo and redo over an app's checkpoints (history §1). States: `steps` (newest first), `redoable`, `canUndo`,
 * `undoBlocked`, `canRedo`, `running`, `checkpointing`, `lastMove`, `checkpointsOff`. Methods: `run`, `group`,
 * `beginGroup`/`endGroup`, `undo`, `redo`, `undoTo`, `redoTo`, `moveTo`, `undoTurn`, `invalidate`, `clear`, `setCheckpointing`.
 */
export class HistoryModel<C = unknown> {
  readonly #port: HistoryPort<C>;
  readonly #options: Required<Omit<HistoryOptions, 'now' | 'measure'>> & {
    now: () => number;
    measure: () => number;
  };
  #entries: Entry<C>[] = [];
  readonly #steps = new State<readonly HistoryStep[]>([]);
  readonly #running = new State(0);
  readonly #restoring = new State(false);
  readonly #checkpointing = new State(true);
  /** Checkpoints in a row over the time budget. */
  #slow = 0;
  /** When the last checkpoint over the hard limit was taken; none since a clear. */
  #overHardLimitAt: number | undefined;
  readonly #undoBlocked = new State<string | undefined>(NothingToUndo);
  readonly #canUndo = new State(false);
  readonly #lastMove = new State<HistoryMove | undefined>(undefined);
  /** What Undo took back, newest first, each with the state it left: Redo puts it back. */
  #redo: { entry: Entry<C>; after: C }[] = [];
  readonly #redoable = new State<readonly HistoryStep[]>([]);
  readonly #checkpointsOff = new State<{ readonly ms: number; readonly bytes: number } | undefined>(undefined);
  readonly #open = new Map<string, OpenGroup>();
  readonly #active = new Map<number, Active>();
  #serial = 0;
  #contexts = 0;
  #lock: Promise<void> = Promise.resolve();
  readonly #canRedo = new State(false);

  constructor(port: HistoryPort<C>, options: HistoryOptions = {}) {
    this.#port = port;
    this.#options = {
      capacity: options.capacity ?? 100,
      checkpoints: Math.max(1, Math.min(options.checkpoints ?? 10, options.capacity ?? 100)),
      budget: options.budget ?? { ms: 50, bytes: 5_000_000 },
      hardLimit: options.hardLimit ?? { ms: 500, bytes: 50_000_000 },
      now: options.now ?? Date.now,
      measure: options.measure ?? (() => performance.now()),
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

  get canUndo(): ReadonlyState<boolean> {
    return this.#canUndo.asReadonly;
  }

  /** Runs and groups in flight, queued ones included. */
  get running(): ReadonlyState<number> {
    return this.#running.asReadonly;
  }

  /** Per-step checkpoints are on; off for a large scene (history §3). */
  get checkpointing(): ReadonlyState<boolean> {
    return this.#checkpointing.asReadonly;
  }

  get canRedo(): ReadonlyState<boolean> {
    return this.#canRedo.asReadonly;
  }

  /** The steps Redo would put back, the next first. */
  get redoable(): ReadonlyState<readonly HistoryStep[]> {
    return this.#redoable.asReadonly;
  }

  /** What Undo or Redo last did (a new value each time). */
  get lastMove(): ReadonlyState<HistoryMove | undefined> {
    return this.#lastMove.asReadonly;
  }

  /** The measurement that turned checkpoints off, the last time they turned off (say so once). */
  get checkpointsOff(): ReadonlyState<{ readonly ms: number; readonly bytes: number } | undefined> {
    return this.#checkpointsOff.asReadonly;
  }

  /** Whether the history is changing the state itself (a run or a restore): an app's change detector ignores it. */
  get applying(): boolean {
    return this.#active.size > 0 || this.#restoring.value;
  }

  /** Whether a context is still in progress (work holding a stale one runs on its own). */
  isActive(context: HistoryContext | undefined): boolean {
    return !!context && this.#active.has(context.id);
  }

  /** Turns per-step checkpoints on or off (a setting); on again after a clear. */
  setCheckpointing(on: boolean): void {
    // Turned on by the person: try again at once, past the pause after a scene over the hard limit.
    if (on) this.#overHardLimitAt = undefined;
    this.#checkpointing.set(on);
  }

  /**
   * Runs one piece of work as a step: in order with the other top-level work, a checkpoint before it, recorded when it
   * succeeds and changed something, its checkpoint restored when it fails. With `within`, it is part of that run or
   * group. An open group of its origin (an agent's turn) gets it as one more run.
   */
  async run<T>(
    label: string | ((value: T) => string),
    origin: HistoryOrigin,
    work: (context: HistoryContext) => Promise<T>,
    options: RunOptions<T> = {},
  ): Promise<StepOutcome<T>> {
    const parent = options.within && this.#active.get(options.within.id);
    if (parent) {
      // Inside a step, work that changes nothing undoable is no run of it either.
      const record = options.undoable === false ? () => false : options.record;
      return this.#nested(parent, 'run', label, work, options.failed, record, 'rollback');
    }

    this.#running.set(this.#running.value + 1);
    this.#updateBlocked();
    const release = await this.#acquire();
    const active = this.#begin('run', undefined, 'rollback');
    try {
      const key = historyGroupKey(origin);
      const group = key ? this.#open.get(key) : undefined;
      this.#closeOthers(key);
      const undoable = options.undoable !== false;
      const checkpoint = undoable ? this.#take(!!group) : { value: undefined };
      let value: T;
      try {
        value = await work(active);
      } catch (error) {
        await this.#restore(checkpoint.value);
        throw error;
      }
      if (options.failed?.(value)) {
        const restored = await this.#restore(checkpoint.value);
        return { value, failed: true, restored };
      }
      if (!undoable || options.record?.(value) === false) return { value, failed: false, restored: false };
      const runs = [labelOf(label, value), ...active.runs];
      const step = this.#record(group ? group.label : labelOf(label, value), origin, runs, checkpoint, key);
      return { value, failed: false, restored: false, step };
    } finally {
      this.#end(active);
      this.#running.set(this.#running.value - 1);
      release();
      this.#updateBlocked();
    }
  }

  /**
   * Runs `work` as one step: one checkpoint (below the hard limit, even when per-step checkpoints are off), every run
   * inside it part of the step. It fails as a whole (its checkpoint restored) when the work throws, returns a value
   * `failed` says is one, or a run inside failed (unless `onFailure: 'continue'`). In an agent's turn it joins the
   * turn's entry.
   */
  async group<T>(
    label: string,
    origin: HistoryOrigin,
    work: (context: HistoryContext) => Promise<T>,
    options: GroupOptions<T> = {},
  ): Promise<GroupOutcome<T>> {
    const onFailure = options.onFailure ?? 'rollback';
    const parent = options.within && this.#active.get(options.within.id);
    if (parent) return this.#nested(parent, 'group', label, work, options.failed, undefined, onFailure);

    this.#running.set(this.#running.value + 1);
    this.#updateBlocked();
    const release = await this.#acquire();
    const active = this.#begin('group', undefined, onFailure);
    try {
      const key = historyGroupKey(origin);
      this.#closeOthers(key);
      const checkpoint = this.#take(true);
      let value: T;
      try {
        value = await work(active);
      } catch (error) {
        await this.#restore(checkpoint.value);
        throw error;
      }
      const failed = !!options.failed?.(value) || (onFailure === 'rollback' && active.failures > 0);
      if (failed) {
        const restored = await this.#restore(checkpoint.value);
        return { value, failed: true, restored, runs: active.runs, failures: active.failures };
      }
      if (!active.runs.length) return { value, failed: false, restored: false, runs: [], failures: active.failures };
      const turn = key ? this.#open.get(key) : undefined;
      const step = this.#record(turn ? turn.label : label, origin, active.runs, checkpoint, key);
      return { value, failed: false, restored: false, step, runs: active.runs, failures: active.failures };
    } finally {
      this.#end(active);
      this.#running.set(this.#running.value - 1);
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
    // Disabled while anything runs or waits to: it never queues behind a run.
    const early = this.#blockedReason();
    if (early) return { ok: false, reason: early };
    const release = await this.#acquire();
    try {
      return await this.#undoTop();
    } finally {
      release();
      this.#updateBlocked();
    }
  }

  /** Puts back the step the last Undo took back; blocked while a run is in flight. */
  async redo(): Promise<UndoResult> {
    if (this.#busy()) return { ok: false, reason: Running };
    if (!this.#redo.length) return { ok: false, reason: NothingToRedo };
    const release = await this.#acquire();
    try {
      return await this.#redoNext();
    } finally {
      release();
      this.#updateBlocked();
    }
  }

  /**
   * Undoes the step `id` and every later one, newest first ("Undo to here"), each as an Undo (so Redo puts them back);
   * stops at the first that cannot be undone and says why. Returns the last step undone.
   */
  async undoTo(id: number): Promise<UndoResult> {
    if (!this.#entries.some((e) => e.step.id === id)) return { ok: false, reason: 'That step is not in the history' };
    const early = this.#blockedReason();
    if (early) return { ok: false, reason: early };
    const release = await this.#acquire();
    try {
      let last: UndoResult = { ok: false, reason: NothingToUndo };
      while (this.#entries.some((e) => e.step.id === id)) {
        const blocked = this.#blockedReason();
        if (blocked) return { ok: false, reason: blocked };
        last = await this.#undoTop();
        if (!last.ok) return last;
      }
      return last;
    } finally {
      release();
      this.#updateBlocked();
    }
  }

  /** Redoes up to and including the undone step `id` ("Redo to here"), each as a Redo. Returns the last step redone. */
  async redoTo(id: number): Promise<UndoResult> {
    if (this.#busy()) return { ok: false, reason: Running };
    if (!this.#redo.some((r) => r.entry.step.id === id)) return { ok: false, reason: NothingToRedo };
    const release = await this.#acquire();
    try {
      let last: UndoResult = { ok: false, reason: NothingToRedo };
      while (this.#redo.some((r) => r.entry.step.id === id)) {
        last = await this.#redoNext();
        if (!last.ok) return last;
      }
      return last;
    } finally {
      release();
      this.#updateBlocked();
    }
  }

  /**
   * Moves the position to just after the step `id` (undefined: before every step), undoing the later steps or redoing
   * the undone ones up to it: History's "Undo to here" and "Redo to here".
   */
  async moveTo(id: number | undefined): Promise<UndoResult> {
    const steps = this.#steps.value;
    if (id === undefined) {
      const oldest = steps[steps.length - 1];
      return oldest ? this.undoTo(oldest.id) : { ok: false, reason: NothingToUndo };
    }
    const i = steps.findIndex((s) => s.id === id);
    if (i === 0) return { ok: false, reason: 'That is the latest step' };
    if (i > 0) return this.undoTo(steps[i - 1].id);
    return this.redoTo(id);
  }

  /**
   * Undoes a turn's entries. While they are the latest, it restores the checkpoint before the turn's first one; with
   * later steps above them it asks first (`laterSteps`), and with `discardLater` restores past those too.
   */
  async undoTurn(turn: string, options: { discardLater?: boolean } = {}): Promise<UndoResult> {
    if (this.#busy()) return { ok: false, reason: Running };
    const release = await this.#acquire();
    try {
      if (this.#active.size > 0) return { ok: false, reason: Running };
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
      const unsafe = this.#entries.slice(0, last + 1).find((e) => e.step.unsafe || e.checkpoint === undefined);
      if (unsafe) return { ok: false, reason: unsafe.step.unsafe ?? 'No checkpoint' };
      const target = this.#entries[last];
      await this.#restore(target.checkpoint);
      const undone = this.#entries.slice(0, last + 1);
      this.#entries = this.#entries.slice(last + 1);
      batch(() => {
        // A turn undone at once is not redone step by step.
        this.#setRedo([]);
        this.#publish();
        this.#lastMove.set({ kind: 'undo', steps: undone.map((e) => e.step) });
      });
      return { ok: true, step: target.step };
    } finally {
      release();
      this.#updateBlocked();
    }
  }

  /** A change outside the history: every step taken before it can no longer be restored safely. History says so. */
  invalidate(reason: string): void {
    batch(() => {
      // Redo would put back a state from before the change: no longer.
      this.#setRedo([]);
      let changed = false;
      for (const e of this.#entries) {
        e.openKey = undefined;
        if (e.step.unsafe) continue;
        e.step = stepWith({ ...fieldsOf(e.step), unsafe: reason }, e.log, e.log.length);
        changed = true;
      }
      if (changed) this.#publish();
      this.#updateBlocked();
    });
  }

  /** Forgets every step (a cleared scene); checkpoints turn back on. */
  clear(): void {
    batch(() => {
      this.#setRedo([]);
      this.#entries = [];
      this.#overHardLimitAt = undefined;
      this.#slow = 0;
      this.#checkpointing.set(true);
      this.#publish();
    });
  }

  dispose(): void {
    this.#entries = [];
    this.#redo = [];
    this.#open.clear();
  }

  /** The latest step undone (the lock held). */
  async #undoTop(): Promise<UndoResult> {
    const blocked = this.#blockedReason();
    if (blocked) return { ok: false, reason: blocked };
    const top = this.#entries[0];
    // The state Undo leaves, for Redo (a checkpoint like any other; none when they are off).
    const after = this.#checkpointing.value ? this.#port.checkpoint() : undefined;
    await this.#restore(top.checkpoint);
    this.#entries = this.#entries.slice(1);
    batch(() => {
      this.#setRedo(after === undefined ? [] : [{ entry: { ...top, openKey: undefined }, after }, ...this.#redo]);
      this.#publish();
      this.#lastMove.set({ kind: 'undo', steps: [top.step] });
    });
    return { ok: true, step: top.step };
  }

  /** The next step redone (the lock held). */
  async #redoNext(): Promise<UndoResult> {
    const next = this.#redo[0];
    if (!next) return { ok: false, reason: NothingToRedo };
    await this.#restore(next.after);
    batch(() => {
      this.#setEntries([next.entry, ...this.#entries]);
      this.#setRedo(this.#redo.slice(1));
      this.#publish();
      this.#lastMove.set({ kind: 'redo', steps: [next.entry.step] });
    });
    return { ok: true, step: next.entry.step };
  }

  #busy(): boolean {
    return this.#running.value > 0 || this.#active.size > 0;
  }

  #setRedo(list: { entry: Entry<C>; after: C }[]) {
    this.#redo = list;
    if (this.#redoable.value.length || list.length) this.#redoable.set(list.map((r) => r.entry.step));
    this.#canRedo.set(list.length > 0);
  }

  /** Work inside a run or group in progress: no lock, no step; recorded in the outer one. */
  async #nested<T>(
    parent: Active,
    kind: 'run' | 'group',
    label: string | ((value: T) => string),
    work: (context: HistoryContext) => Promise<T>,
    failed: ((value: T) => boolean) | undefined,
    record: ((value: T) => boolean) | undefined,
    onFailure: 'rollback' | 'continue',
  ): Promise<GroupOutcome<T>> {
    const active = this.#begin(kind, parent, onFailure);
    try {
      // Inside a group that goes on after a failure, each run keeps its own checkpoint to undo just itself.
      const own: { value: C | undefined } = parent.onFailure === 'continue' ? this.#take(true) : { value: undefined };
      let value: T;
      try {
        value = await work(active);
      } catch (error) {
        parent.failures++;
        await this.#restore(own.value);
        throw error;
      }
      const didFail = !!failed?.(value) || (kind === 'group' && onFailure === 'rollback' && active.failures > 0);
      if (didFail) {
        parent.failures++;
        const restored = await this.#restore(own.value);
        return { value, failed: true, restored, runs: active.runs, failures: active.failures };
      }
      const runs =
        kind === 'run'
          ? record?.(value) === false
            ? active.runs
            : [labelOf(label, value), ...active.runs]
          : active.runs;
      for (const r of runs) parent.runs.push(r);
      return { value, failed: false, restored: false, runs, failures: active.failures };
    } finally {
      this.#end(active);
    }
  }

  #begin(kind: 'run' | 'group', parent: Active | undefined, onFailure: 'rollback' | 'continue'): Active {
    const active: Active = { id: ++this.#contexts, kind, parent, onFailure, runs: [], failures: 0 };
    this.#active.set(active.id, active);
    return active;
  }

  #end(active: Active) {
    this.#active.delete(active.id);
  }

  /** An entry takes more runs only while it is its group's and the latest; anything else closes it. */
  #closeOthers(key: string | undefined) {
    const top = this.#entries[0];
    if (top?.openKey && top.openKey !== key) top.openKey = undefined;
  }

  /** Records a step, or adds the runs to the open entry of its turn when that is the latest. */
  #record(
    label: string,
    origin: HistoryOrigin,
    runs: readonly string[],
    checkpoint: { value: C | undefined; unsafe?: string },
    key: string | undefined,
  ): HistoryStep {
    const top = this.#entries[0];
    const group = key ? this.#open.get(key) : undefined;
    if (group && top?.openKey === key) {
      for (const run of runs) top.log.push(run);
      top.step = stepWith(fieldsOf(top.step), top.log, top.log.length);
      this.#publish();
      return top.step;
    }
    // Its own list of runs (the caller's may be a group's outcome): appended to while its group goes on.
    const log = [...runs];
    const step = stepWith(
      {
        id: ++this.#serial,
        label,
        origin,
        time: this.#options.now(),
        ...(checkpoint.unsafe ? { unsafe: checkpoint.unsafe } : {}),
      },
      log,
      log.length,
    );
    batch(() => {
      // A new step: what Undo took back can no longer be put back.
      this.#setRedo([]);
      this.#setEntries([
        { step, checkpoint: checkpoint.value, openKey: group ? key : undefined, log },
        ...this.#entries,
      ]);
      this.#publish();
    });
    return step;
  }

  /** The entries, the oldest past the capacity dropped and those past the checkpoints' count without one. */
  #setEntries(entries: Entry<C>[]) {
    const { capacity, checkpoints } = this.#options;
    if (entries.length > capacity) entries.length = capacity;
    for (let i = checkpoints; i < entries.length; i++) {
      const e = entries[i];
      if (e.checkpoint === undefined) continue;
      e.checkpoint = undefined;
      e.openKey = undefined;
      if (!e.step.unsafe) {
        const unsafe = `Too old to undo: the last ${checkpoints} step${checkpoints === 1 ? '' : 's'} can be undone`;
        e.step = stepWith({ ...fieldsOf(e.step), unsafe }, e.log, e.log.length);
      }
    }
    this.#entries = entries;
  }

  /** Takes a checkpoint unless checkpoints are off (a group's, and a turn's run's, are taken below the hard limit). */
  #take(group: boolean): { value: C | undefined; unsafe?: string } {
    if (!this.#checkpointing.value && !group) {
      return { value: undefined, unsafe: 'Undo is off for this large scene' };
    }
    const { budget, hardLimit, now, measure } = this.#options;
    // Over the hard limit a moment ago: the snapshot would be thrown away again (a turn's every call paid for it).
    if (this.#overHardLimitAt !== undefined && now() - this.#overHardLimitAt < OverHardLimitPauseMs) {
      return { value: undefined, unsafe: 'This scene is too large to undo' };
    }
    // The size pass is part of the checkpoint's cost: the time budget counts it.
    const t0 = measure();
    const value = this.#port.checkpoint();
    const bytes = this.#port.size?.(value) ?? 0;
    const ms = measure() - t0;
    if (ms > hardLimit.ms || bytes > hardLimit.bytes) {
      this.#overHardLimitAt = now();
      this.#turnOff(ms, bytes);
      return { value: undefined, unsafe: 'This scene is too large to undo' };
    }
    this.#overHardLimitAt = undefined;
    // Size is the scene's; time can be a one-off (a cold start, a pause): only a run of slow checkpoints counts.
    this.#slow = ms > budget.ms ? this.#slow + 1 : 0;
    if (bytes > budget.bytes || this.#slow >= SlowCheckpointsToTurnOff) this.#turnOff(ms, bytes);
    return { value };
  }

  #turnOff(ms: number, bytes: number) {
    if (!this.#checkpointing.value) return;
    batch(() => {
      this.#checkpointing.set(false);
      this.#checkpointsOff.set({ ms, bytes });
    });
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

  #publish() {
    this.#steps.set(this.#entries.map((e) => e.step));
    this.#updateBlocked();
  }

  #blockedReason(): string | undefined {
    const top = this.#entries[0];
    if (this.#busy()) return Running;
    if (!top) return NothingToUndo;
    if (top.step.unsafe) return top.step.unsafe;
    if (top.checkpoint === undefined) return 'No checkpoint';
    return undefined;
  }

  #updateBlocked() {
    const reason = this.#blockedReason();
    batch(() => {
      this.#undoBlocked.set(reason);
      this.#canUndo.set(reason === undefined);
    });
  }

  #acquire(): Promise<() => void> {
    let release!: () => void;
    const next = new Promise<void>((resolve) => (release = resolve));
    const previous = this.#lock;
    this.#lock = previous.then(() => next);
    return previous.then(() => release);
  }
}
