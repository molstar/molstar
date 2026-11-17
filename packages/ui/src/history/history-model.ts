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
//   neither. One that is not recorded (`recording` off, `recordIf` false) is a change outside the history: no
//   checkpoint, no step, and every earlier step unsafe (history §1).
// - Redo: Undo checkpoints the state it leaves, and Redo restores it, until a new step clears the list.
// - Limits (history §5): steps kept, checkpoints kept, redo steps kept, runs listed in a step; past one, the oldest go
//   and their checkpoints are released.

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

/**
 * Why a step cannot be undone, as a code an app can say its own way (its text is `unsafe`): `changed` outside the
 * history (`invalidate`), `not-recorded` (a run while recording was off), `too-old` (past the checkpoints kept),
 * `checkpoints-off` (per-step checkpoints off), `too-large` (over the hard limit).
 */
export type HistoryUnsafeCode = 'changed' | 'not-recorded' | 'too-old' | 'checkpoints-off' | 'too-large';

export interface HistoryStep {
  readonly id: number;
  /** What it did ("Change background"): History's row and the Undo button's tooltip. */
  readonly label: string;
  readonly origin: HistoryOrigin;
  /** When it started (ms, the model's `now`). */
  readonly time: number;
  /**
   * The runs listed, in order (their labels): the latest `maxRuns` of a group's (history §5); a single run lists itself.
   * Copied out when first read.
   */
  readonly runs: readonly string[];
  /** How many runs the step had, those past `maxRuns` (no longer listed) included. */
  readonly runCount: number;
  /** The runs past `maxRuns`, no longer listed (`runCount` less those listed). */
  readonly droppedRuns: number;
  /** The latest `count` runs listed, oldest first, without copying the rest (a view's capped rows). */
  lastRuns(count: number): readonly string[];
  /** Why Undo cannot restore it (a change outside the history after it, too large); undefined: it can. */
  readonly unsafe?: string;
  readonly unsafeCode?: HistoryUnsafeCode;
}

/** A run or group in progress: work that passes it back is part of that step. Opaque. */
export interface HistoryContext {
  readonly id: number;
}

/**
 * Why an Undo or a Redo was refused, as a code an app can map to its own words (`reason` is the labels' text):
 * `running` (a run in flight: transient), `nothing-to-undo`, `nothing-to-redo`, `not-in-history` (an unknown step),
 * `latest` (moving to the position), `turn-empty`, `later-steps` (`undoTurn` asks first), `unsafe` (the step's reason),
 * `no-checkpoint`.
 */
export type HistoryRefusalCode =
  | 'running'
  | 'nothing-to-undo'
  | 'nothing-to-redo'
  | 'not-in-history'
  | 'latest'
  | 'turn-empty'
  | 'later-steps'
  | 'unsafe'
  | 'no-checkpoint';

export type UndoResult =
  | { ok: true; step: HistoryStep }
  | { ok: false; code: HistoryRefusalCode; reason: string; laterSteps?: number };

/** The history's words (history §6): the refusals' and the unsafe steps' text. Domain-neutral by default. */
export interface HistoryLabels {
  readonly running: string;
  readonly nothingToUndo: string;
  readonly nothingToRedo: string;
  readonly notInHistory: string;
  readonly latest: string;
  readonly turnEmpty: string;
  readonly laterSteps: (count: number) => string;
  readonly noCheckpoint: string;
  /** A step past the checkpoints kept. */
  readonly tooOld: (kept: number) => string;
  /** A step recorded while per-step checkpoints were off. */
  readonly checkpointsOff: string;
  /** A checkpoint over the hard limit. */
  readonly tooLarge: string;
  /** A step before a run that was not recorded. */
  readonly notRecorded: string;
}

export const DefaultHistoryLabels: HistoryLabels = {
  running: 'A change is running',
  nothingToUndo: 'Nothing to undo',
  nothingToRedo: 'Nothing to redo',
  notInHistory: 'That step is not in the history',
  latest: 'That is the latest step',
  turnEmpty: 'That turn changed nothing that can be undone',
  laterSteps: (n) => `${n} later step${n === 1 ? '' : 's'} would be undone too`,
  noCheckpoint: 'No checkpoint to restore',
  tooOld: (kept) => `Too old to undo: the last ${kept} step${kept === 1 ? '' : 's'} can be undone`,
  checkpointsOff: 'Undo is off: checkpoints are too large or too slow',
  tooLarge: 'Too large to undo',
  notRecorded: 'Changed while the history was not recording',
};

/** What Undo or Redo last did: a new object each time (a step undone again is told again). */
export interface HistoryMove {
  readonly kind: 'undo' | 'redo';
  /** The steps it took back or put back (an undone turn: several, newest first). */
  readonly steps: readonly HistoryStep[];
}

/** What the history keeps (history §5); past a limit the oldest go, their checkpoints released. */
export interface HistoryLimits {
  /** Undo steps kept (100); the oldest go first. */
  readonly capacity: number;
  /** The latest steps that keep their checkpoint (10, at most `capacity`); older ones stay listed, unsafe ("too old"). */
  readonly checkpoints: number;
  /** Redo steps kept (100; 0: no Redo, and Undo takes no checkpoint of the state it leaves). */
  readonly maxRedo: number;
  /** Runs listed in one step (1000): past it, a group's or a turn's oldest runs leave its list (`runCount` keeps the count). */
  readonly maxRuns: number;
}

export const DefaultHistoryLimits: HistoryLimits = { capacity: 100, checkpoints: 10, maxRedo: 100, maxRuns: 1000 };

export interface HistoryOptions extends Partial<HistoryLimits> {
  /** Over either, per-step checkpoints turn off (history §3). */
  readonly budget?: { readonly ms: number; readonly bytes: number };
  /** Over either, even a group's one checkpoint is not taken. */
  readonly hardLimit?: { readonly ms: number; readonly bytes: number };
  /** The clock of a step's `time` and of the hard limit's pause (ms; `Date.now`). */
  readonly now?: () => number;
  /** The clock that times a checkpoint (ms; `performance.now`). */
  readonly measure?: () => number;
  /** Whether runs become steps (true; `setRecording`). Off, they run in order as changes outside the history (§1). */
  readonly recording?: boolean;
  /**
   * Whether a top-level run or group is recorded, by its origin and label (undefined for a label made from the value):
   * false is as if recording were off for it. Asked before it runs.
   */
  readonly recordIf?: (origin: HistoryOrigin, label: string | undefined) => boolean;
  /** The refusals' and unsafe steps' words (`DefaultHistoryLabels`). */
  readonly labels?: Partial<HistoryLabels>;
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

interface Unsafe {
  readonly text: string;
  readonly code: HistoryUnsafeCode;
}

interface Checkpoint<C> {
  readonly value: C | undefined;
  readonly unsafe?: Unsafe;
  /** The history's epoch when it was taken (history §1): a clear or an invalidate since makes it stale. */
  readonly epoch: number;
}

interface Entry<C> {
  step: HistoryStep;
  checkpoint: C | undefined;
  /** An open group's key while the group can still take runs. */
  openKey?: string;
  /** The step's runs, appended in place while its group takes more: each published step reads its own length. */
  log: string[];
  /** Runs dropped before `log[0]` (past `maxRuns`, history §5). */
  base: number;
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

type StepFields = Omit<HistoryStep, 'runs' | 'runCount' | 'droppedRuns' | 'lastRuns'>;

/** A step's fields but its runs (read without copying them). */
function fieldsOf({ id, label, origin, time, unsafe, unsafeCode }: HistoryStep): StepFields {
  return unsafe === undefined ? { id, label, origin, time } : { id, label, origin, time, unsafe, unsafeCode };
}

/**
 * A step whose `runs` are the latest `maxRuns` of the first `length` of `log`, an append-only list: copied out the first
 * time they are read (and kept), so a step published per run costs nothing until a reader asks, and an older step never
 * sees later runs (a copy per run was R² for a long turn; history §2).
 */
function stepWith(
  fields: StepFields,
  log: readonly string[],
  base: number,
  length: number,
  maxRuns: number,
): HistoryStep {
  const from = Math.max(0, length - maxRuns);
  let runs: readonly string[] | undefined;
  const all = () => (runs ??= from === length ? [] : log.slice(from, length));
  return {
    ...fields,
    runCount: base + length,
    droppedRuns: base + from,
    get runs() {
      return all();
    },
    lastRuns: (count) => (count >= length - from ? all() : log.slice(Math.max(from, length - count), length)),
  };
}

/** Checkpoints over the time budget in a row that turn per-step checkpoints off (history §3). */
const SlowCheckpointsToTurnOff = 3;

/** After a checkpoint over the hard limit, how long none is even taken (history §3): each try is a snapshot thrown away. */
export const OverHardLimitPauseMs = 30_000;

function limitsOf(l: Partial<HistoryLimits>, base: HistoryLimits): HistoryLimits {
  const capacity = Math.max(1, Math.floor(l.capacity ?? base.capacity));
  return {
    capacity,
    checkpoints: Math.max(1, Math.min(Math.floor(l.checkpoints ?? base.checkpoints), capacity)),
    maxRedo: Math.max(0, Math.floor(l.maxRedo ?? base.maxRedo)),
    maxRuns: Math.max(1, Math.floor(l.maxRuns ?? base.maxRuns)),
  };
}

/**
 * Undo and redo over an app's checkpoints (history §1). States: `steps` (newest first), `redoable`, `canUndo`,
 * `undoBlocked` (and its code), `canRedo`, `running`, `checkpointing`, `recording`, `lastMove`, `checkpointsOff`.
 * Methods: `run`, `group`, `beginGroup`/`endGroup`, `undo`, `redo`, `undoTo`, `redoTo`, `moveTo`, `undoTurn`,
 * `invalidate`, `clear`, `setCheckpointing`, `setRecording`, `setLimits`.
 */
export class HistoryModel<C = unknown> {
  readonly #port: HistoryPort<C>;
  readonly #options: {
    budget: { readonly ms: number; readonly bytes: number };
    hardLimit: { readonly ms: number; readonly bytes: number };
    now: () => number;
    measure: () => number;
    recordIf?: (origin: HistoryOrigin, label: string | undefined) => boolean;
  };
  #limits: HistoryLimits;
  readonly #labels: HistoryLabels;
  #entries: Entry<C>[] = [];
  readonly #steps = new State<readonly HistoryStep[]>([]);
  readonly #running = new State(0);
  readonly #restoring = new State(false);
  readonly #checkpointing = new State(true);
  readonly #recording: State<boolean>;
  /** Checkpoints in a row over the time budget. */
  #slow = 0;
  /** When the last checkpoint over the hard limit was taken; none since a clear. */
  #overHardLimitAt: number | undefined;
  readonly #undoBlocked = new State<string | undefined>(undefined);
  readonly #undoBlockedCode = new State<HistoryRefusalCode | undefined>(undefined);
  readonly #canUndo = new State(false);
  readonly #lastMove = new State<HistoryMove | undefined>(undefined);
  /** What Undo took back, newest first, each with the state it left: Redo puts it back. */
  #redo: { entry: Entry<C>; after: C }[] = [];
  readonly #redoable = new State<readonly HistoryStep[]>([]);
  readonly #checkpointsOff = new State<{ readonly ms: number; readonly bytes: number } | undefined>(undefined);
  readonly #open = new Map<string, OpenGroup>();
  readonly #active = new Map<number, Active>();
  /** Bumped by `clear` and `invalidate` (history §1): a checkpoint taken before is stale. */
  #epoch = 0;
  /** The epoch the last `clear` started: a run whose checkpoint is older records nothing. */
  #clearedAt = 0;
  /** Why the last `invalidate` made the steps unsafe (a stale run's step says it too). */
  #invalid: Unsafe | undefined;
  #serial = 0;
  #contexts = 0;
  #lock: Promise<void> = Promise.resolve();
  readonly #canRedo = new State(false);

  constructor(port: HistoryPort<C>, options: HistoryOptions = {}) {
    this.#port = port;
    this.#options = {
      budget: options.budget ?? { ms: 50, bytes: 5_000_000 },
      hardLimit: options.hardLimit ?? { ms: 500, bytes: 50_000_000 },
      now: options.now ?? Date.now,
      measure: options.measure ?? (() => performance.now()),
      recordIf: options.recordIf,
    };
    this.#limits = limitsOf(options, DefaultHistoryLimits);
    this.#labels = { ...DefaultHistoryLabels, ...options.labels };
    this.#recording = new State(options.recording ?? true);
    this.#updateBlocked();
  }

  /** The undoable steps, newest first. */
  get steps(): ReadonlyState<readonly HistoryStep[]> {
    return this.#steps.asReadonly;
  }

  /** Why Undo cannot run now (a run in flight, nothing to undo, an unsafe step); undefined: it can. */
  get undoBlocked(): ReadonlyState<string | undefined> {
    return this.#undoBlocked.asReadonly;
  }

  /** `undoBlocked` as a code: `running` is transient (a view greys its buttons and says nothing more). */
  get undoBlockedCode(): ReadonlyState<HistoryRefusalCode | undefined> {
    return this.#undoBlockedCode.asReadonly;
  }

  get canUndo(): ReadonlyState<boolean> {
    return this.#canUndo.asReadonly;
  }

  /** Runs and groups in flight, queued ones included. */
  get running(): ReadonlyState<number> {
    return this.#running.asReadonly;
  }

  /** Per-step checkpoints are on; off for a large state (history §3). */
  get checkpointing(): ReadonlyState<boolean> {
    return this.#checkpointing.asReadonly;
  }

  /** Runs become steps; off, they are changes outside the history (history §1). */
  get recording(): ReadonlyState<boolean> {
    return this.#recording.asReadonly;
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

  /** The measurement that turned checkpoints off, the last time they turned off (say so once); cleared by `clear`. */
  get checkpointsOff(): ReadonlyState<{ readonly ms: number; readonly bytes: number } | undefined> {
    return this.#checkpointsOff.asReadonly;
  }

  /** The words the history says (its options' over `DefaultHistoryLabels`): a command or a view says the same. */
  get labels(): HistoryLabels {
    return this.#labels;
  }

  /** What the history keeps now (history §5). */
  get limits(): HistoryLimits {
    return this.#limits;
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
    // Turned on by the person: try again at once, past the pause after a state over the hard limit.
    if (on) this.#overHardLimitAt = undefined;
    this.#checkpointing.set(on);
  }

  /**
   * Starts or stops recording (history §1). Off, a run or group still runs in order, but takes no checkpoint and
   * records no step: a change outside the history, so the steps before it turn unsafe once it changes something.
   * Unlike `undoable: false` (a run that changes nothing the history restores), which leaves them as they are.
   */
  setRecording(on: boolean): void {
    this.#recording.set(on);
  }

  /**
   * Changes what the history keeps (history §5), at once: past a lowered limit the oldest steps, checkpoints, redo
   * steps and runs go, their checkpoints released.
   */
  setLimits(limits: Partial<HistoryLimits>): void {
    this.#limits = limitsOf(limits, this.#limits);
    batch(() => {
      for (const e of this.#entries) this.#republish(e, true);
      for (const r of this.#redo) this.#republish(r.entry, true);
      this.#setEntries([...this.#entries]);
      this.#setRedo(this.#redo);
      this.#publish();
    });
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
      const recorded = undoable && this.#records(origin, label);
      const checkpoint = recorded ? this.#take(!!group) : this.#none();
      let value: T;
      try {
        value = await work(active);
      } catch (error) {
        await this.#restore(checkpoint);
        throw error;
      }
      if (options.failed?.(value)) {
        const restored = await this.#restore(checkpoint);
        return { value, failed: true, restored };
      }
      if (!undoable || options.record?.(value) === false) return { value, failed: false, restored: false };
      if (!recorded) {
        this.#notRecorded();
        return { value, failed: false, restored: false };
      }
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
      const recorded = this.#records(origin, label);
      const checkpoint = recorded ? this.#take(true) : this.#none();
      let value: T;
      try {
        value = await work(active);
      } catch (error) {
        await this.#restore(checkpoint);
        throw error;
      }
      const failed = !!options.failed?.(value) || (onFailure === 'rollback' && active.failures > 0);
      if (failed) {
        const restored = await this.#restore(checkpoint);
        return { value, failed: true, restored, runs: active.runs, failures: active.failures };
      }
      if (!active.runs.length) return { value, failed: false, restored: false, runs: [], failures: active.failures };
      if (!recorded) {
        this.#notRecorded();
        return { value, failed: false, restored: false, runs: active.runs, failures: active.failures };
      }
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
    const early = this.#blocked();
    if (early) return early;
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
    if (this.#busy()) return this.#refuse('running');
    if (!this.#redo.length) return this.#refuse('nothing-to-redo');
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
    if (!this.#entries.some((e) => e.step.id === id)) return this.#refuse('not-in-history');
    const early = this.#blocked();
    if (early) return early;
    const release = await this.#acquire();
    try {
      let last: UndoResult = this.#refuse('nothing-to-undo');
      while (this.#entries.some((e) => e.step.id === id)) {
        const blocked = this.#blocked();
        if (blocked) return blocked;
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
    if (this.#busy()) return this.#refuse('running');
    if (!this.#redo.some((r) => r.entry.step.id === id)) return this.#refuse('nothing-to-redo');
    const release = await this.#acquire();
    try {
      let last: UndoResult = this.#refuse('nothing-to-redo');
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
      return oldest ? this.undoTo(oldest.id) : this.#refuse('nothing-to-undo');
    }
    const i = steps.findIndex((s) => s.id === id);
    if (i === 0) return this.#refuse('latest');
    if (i > 0) return this.undoTo(steps[i - 1].id);
    return this.redoTo(id);
  }

  /**
   * Undoes a turn's entries: the origin's (kind and turn, as they group: `historyGroupKey`), or a turn id alone for an
   * `agent`'s. While they are the latest, it restores the checkpoint before the turn's first one; with later steps
   * above them it asks first (`laterSteps`), and with `discardLater` restores past those too.
   */
  async undoTurn(turn: HistoryOrigin | string, options: { discardLater?: boolean } = {}): Promise<UndoResult> {
    const key = historyGroupKey(typeof turn === 'string' ? { kind: 'agent', turn } : turn);
    if (!key) return this.#refuse('turn-empty');
    if (this.#busy()) return this.#refuse('running');
    const release = await this.#acquire();
    try {
      if (this.#active.size > 0) return this.#refuse('running');
      const ofTurn = (e: Entry<C>) => historyGroupKey(e.step.origin) === key;
      const indices = this.#entries.flatMap((e, i) => (ofTurn(e) ? [i] : []));
      if (!indices.length) return this.#refuse('turn-empty');
      const last = indices[indices.length - 1];
      const later = this.#entries.slice(0, last + 1).filter((e) => !ofTurn(e)).length;
      if (later > 0 && !options.discardLater) {
        return { ok: false, code: 'later-steps', reason: this.#labels.laterSteps(later), laterSteps: later };
      }
      const unsafe = this.#entries.slice(0, last + 1).find((e) => e.step.unsafe || e.checkpoint === undefined);
      if (unsafe)
        return unsafe.step.unsafe ? this.#refuse('unsafe', unsafe.step.unsafe) : this.#refuse('no-checkpoint');
      const target = this.#entries[last];
      const epoch = this.#epoch;
      await this.#restore({ value: target.checkpoint, epoch });
      // Cleared or changed while it restored: the lists are no longer the ones it undid.
      if (epoch !== this.#epoch) return { ok: true, step: target.step };
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

  /**
   * A change outside the history: every step taken before it can no longer be restored safely. History says so. A run
   * in flight whose checkpoint predates it records its step unsafe too, and does not restore it on failure (§1).
   */
  invalidate(reason: string): void {
    this.#invalidate({ text: reason, code: 'changed' });
  }

  /**
   * Forgets every step and what Redo would put back (a cleared scene); checkpoints turn back on and `checkpointsOff`
   * clears. A run in flight records nothing (its checkpoint is from before the clear) and restores nothing on failure;
   * an open turn's next run starts a new entry.
   */
  clear(): void {
    this.#epoch++;
    this.#clearedAt = this.#epoch;
    batch(() => {
      this.#setRedo([]);
      this.#entries = [];
      this.#overHardLimitAt = undefined;
      this.#slow = 0;
      this.#checkpointing.set(true);
      this.#checkpointsOff.set(undefined);
      this.#publish();
    });
  }

  dispose(): void {
    this.#entries = [];
    this.#redo = [];
    this.#open.clear();
  }

  #invalidate(unsafe: Unsafe): void {
    this.#epoch++;
    this.#invalid = unsafe;
    batch(() => {
      // Redo would put back a state from before the change: no longer.
      this.#setRedo([]);
      let changed = false;
      for (const e of this.#entries) {
        e.openKey = undefined;
        if (e.step.unsafe) continue;
        this.#markUnsafe(e, unsafe);
        changed = true;
      }
      if (changed) this.#publish();
      this.#updateBlocked();
    });
  }

  /** A run that changed something but was not recorded: a change outside the history (§1). */
  #notRecorded() {
    this.#invalidate({ text: this.#labels.notRecorded, code: 'not-recorded' });
  }

  #records(origin: HistoryOrigin, label: string | ((value: never) => string)): boolean {
    if (!this.#recording.value) return false;
    const recordIf = this.#options.recordIf;
    return !recordIf || recordIf(origin, typeof label === 'string' ? label : undefined);
  }

  /** The latest step undone (the lock held). */
  async #undoTop(): Promise<UndoResult> {
    const blocked = this.#blocked();
    if (blocked) return blocked;
    const top = this.#entries[0];
    // The state Undo leaves, for Redo (a checkpoint like any other; none when they are off or Redo keeps none).
    const after = this.#checkpointing.value && this.#limits.maxRedo > 0 ? this.#port.checkpoint() : undefined;
    const epoch = this.#epoch;
    await this.#restore({ value: top.checkpoint, epoch });
    if (epoch !== this.#epoch) return { ok: true, step: top.step };
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
    if (!next) return this.#refuse('nothing-to-redo');
    const epoch = this.#epoch;
    await this.#restore({ value: next.after, epoch });
    if (epoch !== this.#epoch) return { ok: true, step: next.entry.step };
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

  /** The redo list, past `maxRedo` the farthest dropped (their checkpoints released). */
  #setRedo(list: { entry: Entry<C>; after: C }[]) {
    if (list.length > this.#limits.maxRedo) list = list.slice(0, this.#limits.maxRedo);
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
      const own = parent.onFailure === 'continue' ? this.#take(true) : this.#none();
      let value: T;
      try {
        value = await work(active);
      } catch (error) {
        parent.failures++;
        await this.#restore(own);
        throw error;
      }
      const didFail = !!failed?.(value) || (kind === 'group' && onFailure === 'rollback' && active.failures > 0);
      if (didFail) {
        parent.failures++;
        const restored = await this.#restore(own);
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

  /**
   * Records a step, or adds the runs to the open entry of its turn when that is the latest. A checkpoint from before a
   * clear records nothing; one from before an invalidate records the step unsafe (§1).
   */
  #record(
    label: string,
    origin: HistoryOrigin,
    runs: readonly string[],
    checkpoint: Checkpoint<C>,
    key: string | undefined,
  ): HistoryStep | undefined {
    if (checkpoint.epoch < this.#clearedAt) return undefined;
    if (checkpoint.epoch !== this.#epoch) checkpoint = { value: undefined, unsafe: this.#invalid, epoch: this.#epoch };
    const top = this.#entries[0];
    const group = key ? this.#open.get(key) : undefined;
    if (group && top?.openKey === key) {
      for (const run of runs) top.log.push(run);
      this.#republish(top);
      this.#publish();
      return top.step;
    }
    // Its own list of runs (the caller's may be a group's outcome): appended to while its group goes on.
    const { maxRuns } = this.#limits;
    const log = runs.length > maxRuns ? runs.slice(runs.length - maxRuns) : [...runs];
    const base = runs.length - log.length;
    const unsafe = checkpoint.unsafe;
    const step = stepWith(
      {
        id: ++this.#serial,
        label,
        origin,
        time: this.#options.now(),
        ...(unsafe ? { unsafe: unsafe.text, unsafeCode: unsafe.code } : {}),
      },
      log,
      base,
      log.length,
      maxRuns,
    );
    batch(() => {
      // A new step: what Undo took back can no longer be put back.
      this.#setRedo([]);
      this.#setEntries([
        { step, checkpoint: checkpoint.value, openKey: group ? key : undefined, log, base },
        ...this.#entries,
      ]);
      this.#publish();
    });
    return step;
  }

  /**
   * The entry's step again over its runs: past twice `maxRuns`, the log keeps its latest `maxRuns` in a new list (the
   * steps published before read the old one), so a long turn holds O(maxRuns) runs and copies each once (history §5);
   * `now`: a lowered limit, at once.
   */
  #republish(e: Entry<C>, now = false) {
    const { maxRuns } = this.#limits;
    if (e.log.length > maxRuns && (now || e.log.length >= 2 * maxRuns)) {
      const drop = e.log.length - maxRuns;
      e.log = e.log.slice(drop);
      e.base += drop;
    }
    e.step = stepWith(fieldsOf(e.step), e.log, e.base, e.log.length, maxRuns);
  }

  #markUnsafe(e: Entry<C>, unsafe: Unsafe) {
    e.step = stepWith(
      { ...fieldsOf(e.step), unsafe: unsafe.text, unsafeCode: unsafe.code },
      e.log,
      e.base,
      e.log.length,
      this.#limits.maxRuns,
    );
  }

  /** The entries, the oldest past the capacity dropped and those past the checkpoints' count without one. */
  #setEntries(entries: Entry<C>[]) {
    const { capacity, checkpoints } = this.#limits;
    if (entries.length > capacity) entries.length = capacity;
    for (let i = checkpoints; i < entries.length; i++) {
      const e = entries[i];
      if (e.checkpoint === undefined) continue;
      e.checkpoint = undefined;
      e.openKey = undefined;
      if (!e.step.unsafe) this.#markUnsafe(e, { text: this.#labels.tooOld(checkpoints), code: 'too-old' });
    }
    this.#entries = entries;
  }

  #none(): Checkpoint<C> {
    return { value: undefined, epoch: this.#epoch };
  }

  /** Takes a checkpoint unless checkpoints are off (a group's, and a turn's run's, are taken below the hard limit). */
  #take(group: boolean): Checkpoint<C> {
    const epoch = this.#epoch;
    if (!this.#checkpointing.value && !group) {
      return { value: undefined, unsafe: { text: this.#labels.checkpointsOff, code: 'checkpoints-off' }, epoch };
    }
    const tooLarge: Unsafe = { text: this.#labels.tooLarge, code: 'too-large' };
    const { budget, hardLimit, now, measure } = this.#options;
    // Over the hard limit a moment ago: the snapshot would be thrown away again (a turn's every call paid for it).
    if (this.#overHardLimitAt !== undefined && now() - this.#overHardLimitAt < OverHardLimitPauseMs) {
      return { value: undefined, unsafe: tooLarge, epoch };
    }
    // The size pass is part of the checkpoint's cost: the time budget counts it.
    const t0 = measure();
    const value = this.#port.checkpoint();
    const bytes = this.#port.size?.(value) ?? 0;
    const ms = measure() - t0;
    if (ms > hardLimit.ms || bytes > hardLimit.bytes) {
      this.#overHardLimitAt = now();
      this.#turnOff(ms, bytes);
      return { value: undefined, unsafe: tooLarge, epoch };
    }
    this.#overHardLimitAt = undefined;
    // Size is the state's; time can be a one-off (a cold start, a pause): only a run of slow checkpoints counts.
    this.#slow = ms > budget.ms ? this.#slow + 1 : 0;
    if (bytes > budget.bytes || this.#slow >= SlowCheckpointsToTurnOff) this.#turnOff(ms, bytes);
    return { value, epoch };
  }

  #turnOff(ms: number, bytes: number) {
    if (!this.#checkpointing.value) return;
    batch(() => {
      this.#checkpointing.set(false);
      this.#checkpointsOff.set({ ms, bytes });
    });
  }

  /** Restores a checkpoint, unless a clear or an invalidate came after it (it would undo that change too, §1). */
  async #restore(checkpoint: Checkpoint<C>): Promise<boolean> {
    if (checkpoint.value === undefined || checkpoint.epoch !== this.#epoch) return false;
    this.#restoring.set(true);
    try {
      await this.#port.restore(checkpoint.value);
      return true;
    } finally {
      this.#restoring.set(false);
    }
  }

  #publish() {
    this.#steps.set(this.#entries.map((e) => e.step));
    this.#updateBlocked();
  }

  #refuse(code: Exclude<HistoryRefusalCode, 'later-steps'>, reason?: string): UndoResult {
    const l = this.#labels;
    const text =
      reason ??
      {
        running: l.running,
        'nothing-to-undo': l.nothingToUndo,
        'nothing-to-redo': l.nothingToRedo,
        'not-in-history': l.notInHistory,
        latest: l.latest,
        'turn-empty': l.turnEmpty,
        unsafe: l.noCheckpoint,
        'no-checkpoint': l.noCheckpoint,
      }[code];
    return { ok: false, code, reason: text };
  }

  #blocked(): UndoResult | undefined {
    const top = this.#entries[0];
    if (this.#busy()) return this.#refuse('running');
    if (!top) return this.#refuse('nothing-to-undo');
    if (top.step.unsafe) return this.#refuse('unsafe', top.step.unsafe);
    if (top.checkpoint === undefined) return this.#refuse('no-checkpoint');
    return undefined;
  }

  #updateBlocked() {
    const blocked = this.#blocked() as Extract<UndoResult, { ok: false }> | undefined;
    batch(() => {
      this.#undoBlocked.set(blocked?.reason);
      this.#undoBlockedCode.set(blocked?.code);
      this.#canUndo.set(blocked === undefined);
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
