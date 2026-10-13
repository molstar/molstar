/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Json } from '@molstar/ai-agent';
import { HarnessError } from '../contracts/ports.js';
import type {
  ActiveTurn,
  CallRecord,
  ControlEvent,
  Execution,
  HarnessEvent,
  HarnessId,
  LaunchKey,
  StepEvent,
  TurnEnd,
  TurnSummary,
  Usage,
} from '../contracts/records.js';
import type { StorageTx } from '../contracts/storage.js';

type Mutable<T> = { -readonly [K in keyof T]: T[K] };

export const ZERO: Usage = {
  modelCalls: 0,
  toolCalls: 0,
  uncachedInput: 0,
  cacheRead: 0,
  cacheWrite: 0,
  output: 0,
  reasoning: 0,
};

export function addUsage(a: Usage, b: Partial<Usage>): Usage {
  const out: Mutable<Usage> = { ...a };
  for (const key of [
    'modelCalls',
    'toolCalls',
    'uncachedInput',
    'cacheRead',
    'cacheWrite',
    'output',
    'reasoning',
  ] as const)
    out[key] = a[key] + (b[key] ?? 0);
  if (a.cost !== undefined || b.cost !== undefined) out.cost = (a.cost ?? 0) + (b.cost ?? 0);
  if (b.estimated === true) out.estimated = true;
  return out;
}

/** Usage of the calls' attempts: all of it for the turn's totals, only top-level calls for its budget. */
export function usageOf(calls: readonly CallRecord[]): { all: Usage; charged: Usage } {
  let all = ZERO;
  let charged = ZERO;
  for (const call of calls)
    for (const attempt of call.attempts) {
      all = addUsage(all, attempt.usage);
      if (call.parentCallId === undefined) charged = addUsage(charged, attempt.usage);
    }
  return { all, charged };
}

/** An error raised after the operation's writes commit (a claim that blocks a turn, for example). */
export class Raise {
  readonly error: HarnessError;
  constructor(error: HarnessError) {
    this.error = error;
  }
}

export interface HostSettings {
  readonly leaseMs: number;
  readonly launchMs: number;
  readonly schemaVersion: number;
}

/** One operation's view of a harness inside a transaction: the Execution, and events with seqs. */
export class Writer {
  readonly exec: Mutable<Execution>;
  readonly written: HarnessEvent[] = [];
  readonly tx: StorageTx;
  readonly now: number;
  readonly settings: HostSettings;
  readonly key: string | undefined;
  readonly inputDigest: string | undefined;
  changed = false;

  constructor(
    tx: StorageTx,
    execution: Execution,
    now: number,
    settings: HostSettings,
    key: string | undefined,
    inputDigest: string | undefined,
  ) {
    this.exec = { ...execution };
    this.tx = tx;
    this.now = now;
    this.settings = settings;
    this.key = key;
    this.inputDigest = inputDigest;
  }

  get harnessId(): HarnessId {
    return this.exec.harnessId;
  }

  set<K extends keyof Execution>(field: K, value: Execution[K] | undefined): void {
    if (value === undefined) delete this.exec[field];
    else this.exec[field] = value;
    this.changed = true;
  }

  private append<E extends HarnessEvent>(event: E): E {
    this.set('seq', event.seq);
    this.tx.events.append(event);
    this.written.push(event);
    return event;
  }

  control(
    control: ControlEvent['control'],
    fields: Partial<Omit<ControlEvent, 'kind' | 'seq' | 'at'>> = {},
  ): ControlEvent {
    return this.append({
      kind: 'control',
      seq: this.exec.seq + 1,
      at: this.now,
      control,
      detail: {},
      ...(this.key !== undefined ? { key: this.key } : {}),
      ...(this.inputDigest !== undefined ? { inputDigest: this.inputDigest } : {}),
      ...fields,
    });
  }

  step(fields: Omit<StepEvent, 'kind' | 'seq' | 'at' | 'commitKey' | 'inputDigest'>): StepEvent {
    return this.append({
      kind: 'step',
      seq: this.exec.seq + 1,
      at: this.now,
      commitKey: this.key ?? '',
      ...(this.inputDigest !== undefined ? { inputDigest: this.inputDigest } : {}),
      ...fields,
    });
  }

  /** Work exists (harness.md §4): an active turn neither parked nor blocked, or no turn and an item waiting. */
  workExists(): boolean {
    const turn = this.exec.turn;
    return turn !== undefined ? turn.park === undefined && turn.blocked === undefined : this.exec.inbox.length > 0;
  }

  /** No runner is on the work: `none`, or a `launching` or lease that has expired at host time. */
  runnerFree(): boolean {
    const runner = this.exec.runner;
    if (runner.phase === 'launching') return runner.until <= this.now;
    if (runner.phase === 'live') return runner.leaseUntil <= this.now;
    return true;
  }

  liveLease(): boolean {
    const runner = this.exec.runner;
    return runner.phase === 'live' && runner.leaseUntil > this.now;
  }

  /** Launches when work exists and no runner is on it (harness.md §4); writes `launched` last. */
  launchIfNeeded(): LaunchKey | undefined {
    if (!this.workExists() || !this.runnerFree()) return undefined;
    const launch = `${this.harnessId}/launch/${this.exec.seq + 1}`;
    this.set('runner', { phase: 'launching', key: launch, until: this.now + this.settings.launchMs });
    this.control('launched', { detail: { launch } });
    return launch;
  }

  /** The active turn, or `stale` if it is gone or another. */
  turn(turnId?: string): Mutable<ActiveTurn> {
    const turn = this.exec.turn;
    if (turn === undefined || (turnId !== undefined && turn.id !== turnId))
      throw new HarnessError('stale', `turn ${turnId ?? ''} is not the active turn`);
    return turn;
  }

  /** The summary of the active turn as it ends: what started it comes from its begin step event. */
  async summary(status: TurnEnd['status'], usage: Usage, lastStep?: string): Promise<TurnSummary> {
    const turn = this.turn();
    const [begin] = await this.tx.events.byTurn(turn.id, undefined, 1);
    const begun = begin?.kind === 'step' ? begin : undefined;
    return {
      id: turn.id,
      status,
      startedBy: begun?.begun?.items ?? [],
      firstStep: begun?.step ?? turn.pending.id,
      ...(lastStep !== undefined ? { lastStep } : {}),
      usage,
      endedAt: this.now,
    };
  }

  /** Ends the active turn by a control event (`stopped`, `reset`). */
  async endByControl(control: 'stopped' | 'reset', reason?: Json): Promise<ControlEvent> {
    const turn = this.turn();
    const summary = await this.summary(control, turn.usage);
    this.set('turn', undefined);
    return this.control(control, { turnId: turn.id, summary, detail: reason !== undefined ? { reason } : {} });
  }
}
