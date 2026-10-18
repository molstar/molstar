/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { digest, estimatedTokens, type Json, type Message } from '../../messages/index.js';
import type {
  ActiveTurn,
  AttemptRecord,
  BlobRef,
  CallId,
  CallRecord,
  Compaction,
  EntryId,
  HarnessEvent,
  HeadId,
  Limits,
  ScopeId,
  StepEvent,
  TurnEnd,
  Usage,
} from '../contracts/records.js';

/** How a step stops early (harness.md §6); the first of stop, fail and park wins. */
export class Exit extends Error {
  readonly kind: 'stop' | 'fail' | 'park';
  readonly reason?: Json;
  readonly path: readonly CallId[];
  readonly note?: Json;
  constructor(kind: Exit['kind'], options: { reason?: Json; path?: readonly CallId[]; note?: Json } = {}) {
    super(`exit: ${kind}`);
    this.kind = kind;
    this.path = options.path ?? [];
    if (options.reason !== undefined) this.reason = options.reason;
    if (options.note !== undefined) this.note = options.note;
  }
}
const RANK: Record<Exit['kind'], number> = { stop: 0, fail: 1, park: 2 };
export const firstExit = (exits: readonly Exit[]): Exit | undefined =>
  [...exits].sort((a, b) => RANK[a.kind] - RANK[b.kind])[0];

/** A step's work ends without a commit: a reset ended its turn (harness.md §5). */
export class Halt extends Error {
  constructor() {
    super('halt: reset');
  }
}

/** What one top-level step's runs share: its abort and its exits. */
export class StepState {
  readonly controller = new AbortController();
  readonly exits: Exit[] = [];
  /** The turn's limits when the step began; a commit records them again if an `extend` replaced them. */
  limits?: Limits;
  /** A reset ended the turn: nothing this step did is written. */
  halted = false;
  /** The keys of the calls that started in an earlier, parked run of this step: run again with `prior: "possible"`. */
  prior?: ReadonlySet<string>;
  /** Records an exit and aborts the step's other work (harness.md §6). */
  exit(exit: Exit): Exit {
    this.exits.push(exit);
    if (!this.controller.signal.aborted) this.controller.abort(exit);
    return exit;
  }
}

/** A call record while its step is open. */
export interface OpenCall {
  readonly id: CallId;
  readonly key: string;
  readonly parentCallId?: CallId;
  readonly kind: 'model' | 'tool';
  readonly request: BlobRef;
  readonly attempts: AttemptRecord[];
  result?: BlobRef;
}

/** A context's transient head within a step: a pointer's (`scope`), or a local one (a fork, or empty). */
export interface Staged {
  readonly scope?: ScopeId;
  /** The head it started from. */
  readonly from?: HeadId;
  entries: EntryId[];
  /** Hashes appended since the step started; undefined after an edit. */
  appended: EntryId[] | undefined;
}

export const local = (entries: readonly EntryId[] = []): Staged => ({
  entries: [...entries],
  appended: [],
});

/** What a top-level step writes at its end (harness.md §2): calls, new messages, pointer moves, compactions, blobs. */
export class Collector {
  readonly calls: OpenCall[] = [];
  readonly messages = new Map<EntryId, Message>();
  readonly pointers = new Map<ScopeId, Staged>();
  readonly compaction = new Map<ScopeId, Compaction>();
  readonly blobs = new Map<BlobRef, Json>();
  /** Who writes each pointer in this step (harness.md §8). */
  readonly holders = new Map<ScopeId, string>();

  stage(scope: ScopeId, from: HeadId | undefined, entries: readonly EntryId[]): Staged {
    let staged = this.pointers.get(scope);
    if (staged === undefined) {
      staged = { scope, ...(from !== undefined ? { from } : {}), ...local(entries) };
      this.pointers.set(scope, staged);
    }
    return staged;
  }

  append(staged: Staged, ids: readonly EntryId[], messages: readonly Message[]): void {
    staged.entries.push(...ids);
    staged.appended?.push(...ids);
    ids.forEach((id, i) => {
      const message = messages[i];
      if (message !== undefined) this.messages.set(id, message);
    });
  }

  edit(staged: Staged, entries: readonly EntryId[]): void {
    staged.entries = [...entries];
    staged.appended = undefined;
  }
}

export const ZERO: Usage = {
  modelCalls: 0,
  toolCalls: 0,
  uncachedInput: 0,
  cacheRead: 0,
  cacheWrite: 0,
  output: 0,
  reasoning: 0,
};

type Mutable<T> = { -readonly [K in keyof T]: T[K] };

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

/** A head's id: the digest of its entries (harness.md §2). */
export const headIdOf = (entries: readonly EntryId[]): HeadId => digest({ head: [...entries] });

/**
 * What a step's parked runs left (harness.md §5): the results of the top-level calls that succeeded and hold no
 * pointer (reused on resume, not run again), and the keys of every call that started.
 */
export interface Parked {
  readonly results: Map<CallId, BlobRef>;
  readonly started: Set<string>;
}

/**
 * The harness's state, one object in memory (harness.md §2–§4): entries by hash, heads as lists of entry ids,
 * pointers, the event log (`seq` is an event's position, from 1), blobs, the active turn. A step changes it only at
 * its end, in one synchronous assignment.
 */
export class Store {
  readonly entries = new Map<EntryId, Message>();
  readonly heads = new Map<HeadId, readonly EntryId[]>();
  /** Pointers other than `main`. */
  readonly pointers = new Map<ScopeId, HeadId>();
  readonly events: HarnessEvent[] = [];
  readonly blobs = new Map<BlobRef, Json>();
  /** How each ended turn ended. */
  readonly ends = new Map<string, TurnEnd>();
  private readonly tokens = new Map<EntryId, number>();
  main: HeadId;
  turn: ActiveTurn | undefined;
  /** The turn's begin step, while it is active. */
  firstStep: string | undefined;
  /** `main` at the end of the last completed turn (reset's default). */
  completedMain: HeadId | undefined;
  /** What the parked runs of the active turn's pending step left, for its re-run on resume (harness.md §5). */
  readonly parked = new Map<string, Parked>();

  constructor() {
    this.main = this.head([]);
  }

  get seq(): number {
    return this.events.length;
  }

  /** Records a head; its id. */
  head(entries: readonly EntryId[]): HeadId {
    const id = headIdOf(entries);
    if (!this.heads.has(id)) this.heads.set(id, [...entries]);
    return id;
  }

  /** The last compaction of a pointer, from the event log (the `boundary` and `turn-start` payloads, harness.md §7). */
  lastCompaction(pointer: ScopeId): { seq: number; sizeAfter: Compaction['sizeAfter'] } | undefined {
    for (let i = this.events.length - 1; i >= 0; i -= 1) {
      const event = this.events[i];
      if (event?.kind !== 'step') continue;
      const found = event.compaction?.findLast((c) => c.pointer === pointer);
      if (found !== undefined) return { seq: event.seq, sizeAfter: found.sizeAfter };
    }
    return undefined;
  }

  /**
   * An input-token estimate of a head (a quarter of its JSON length), cached per entry. A deferred tool declaration
   * counts its name only: a provider shows at most that. A label is never sent, so it never counts.
   */
  tokensOf(entries: readonly EntryId[], message: (id: EntryId) => Message | undefined): number {
    let total = 0;
    for (const id of entries) {
      let tokens = this.tokens.get(id);
      if (tokens === undefined) {
        const found = message(id);
        const shown =
          found?.role === 'system' && found.parts.some((p) => p.kind === 'tool-definition' && p.deferred)
            ? found.parts.map((p) => (p.kind === 'tool-definition' && p.deferred ? p.definition.name : p))
            : found && { role: found.role, parts: found.parts };
        tokens = estimatedTokens(JSON.stringify(shown ?? null));
        if (found !== undefined) this.tokens.set(id, tokens);
      }
      total += tokens;
    }
    return total;
  }

  steps(ids: readonly string[]): StepEvent[] {
    const wanted = new Set(ids);
    return this.events.filter((e): e is StepEvent => e.kind === 'step' && wanted.has(e.step));
  }
}
