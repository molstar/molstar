/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { batch, KeyedState, ReactiveModel, type ReadonlyState, State } from '@molstar/core/reactive';
import { type Issue, issueOf, type Node, nodesEqual } from '@molstar/core/repr/tree';
import { SourceIssues, StoreIssues } from './issues.js';
import { DefaultTileBins, dataKey, fitLevel, nodeKey, tileRegion, tileSpan } from './keys.js';
import {
  type Coverage,
  type DynamicNode,
  type FailureKind,
  isAbortError,
  type OpenedSource,
  type RegionAnswer,
  type RegionRequest,
  type Shape,
  type Tables,
} from './protocol.js';
import { answerBytes, checkTables, isEmptyAnswer } from './tables.js';

// The data store (_spec/README.md "Data", Store): entries by key in keyed states, versioned; requests deduped by
// tile, ordered and bounded in flight, cancelled after a grace, retried by policy; a byte-bounded cache that never
// evicts what is in view; nodes held by what carried them.

/** Where a wanted region stands: in the view, its margin, or prefetched (EC-172's order). */
export type Priority = 'view' | 'margin' | 'prefetch';

const Rank: Record<Priority, number> = { view: 0, margin: 1, prefetch: 2 };
const Unwanted = 3;

/** A region a visible layer needs: positions `[start, end)` of a system of an attached source, at a level. */
export interface Want {
  readonly source: string;
  readonly system: string;
  readonly start: number;
  readonly end: number;
  /** A bin width; an unlisted one falls back to the coarsest listed no wider (EC-170). */
  readonly level: number;
  /** `view` by default. */
  readonly priority?: Priority;
}

/** The tile an entry holds: its parts and its positions. */
export interface TileRegion {
  readonly source: string;
  /** The coverage's name of the system (an alias resolved). */
  readonly system: string;
  readonly level: number;
  readonly tile: number;
  readonly start: number;
  readonly end: number;
}

/** A data entry: pending is never drawn or hit as empty (EC-169); failed carries its issue. */
export type DataEntry =
  | { readonly status: 'pending'; readonly version: number; readonly region: TileRegion }
  | {
      readonly status: 'loaded';
      readonly version: number;
      readonly region: TileRegion;
      readonly tables: Tables;
      /** The node keys (`<source>#<key>`) this region carried. */
      readonly nodes: readonly string[];
      readonly bytes: number;
    }
  | { readonly status: 'empty'; readonly version: number; readonly region: TileRegion }
  | {
      readonly status: 'failed';
      readonly version: number;
      readonly region: TileRegion;
      readonly kind: FailureKind;
      readonly issue: Issue;
    };

export type EntryStatus = DataEntry['status'];

/** A node entry: a dynamic node of a source, held while something holds it. */
export interface NodeEntry {
  readonly source: string;
  readonly key: string;
  readonly node: Node;
  readonly version: number;
}

/** The store's time: real timers by default, a fake one in tests. */
export interface StoreClock {
  now(): number;
  setTimeout(fn: () => void, ms: number): unknown;
  clearTimeout(handle: unknown): void;
}

export const RealClock: StoreClock = {
  now: () => performance.now(),
  setTimeout: (fn, ms) => setTimeout(fn, ms),
  clearTimeout: (h) => clearTimeout(h as ReturnType<typeof setTimeout>),
};

export interface StoreSettings {
  /** Bins per tile. */
  readonly tileBins: number;
  /** The data entries' byte bound. */
  readonly maxBytes: number;
  /** Requests in flight at most (EC-172). */
  readonly maxInFlight: number;
  /** How long a request that left the wanted set stays in flight before it is cancelled (EC-173). */
  readonly graceMs: number;
  /** A `failed` (or malformed) answer: tried `attempts` times with backoff `baseMs · factor^n`, then a cooldown. */
  readonly retry: {
    readonly attempts: number;
    readonly baseMs: number;
    readonly factor: number;
    readonly cooldownMs: number;
  };
}

export const DefaultStoreSettings: StoreSettings = {
  tileBins: DefaultTileBins,
  maxBytes: 64 * 1024 * 1024,
  maxInFlight: 6,
  graceMs: 300,
  retry: { attempts: 3, baseMs: 500, factor: 2, cooldownMs: 30_000 },
};

type Phase = 'queued' | 'flight' | 'backoff' | 'cooldown' | 'done';

/** A tile's bookkeeping (one per entry key). */
interface Slot {
  readonly key: string;
  readonly region: TileRegion;
  readonly shapes: readonly Shape[];
  rank: number;
  /** The order it was first wanted in: ties within a rank go first come, first served. */
  order: number;
  phase: Phase;
  /** Bumped when a request is let go, so its late answer is dropped. */
  token: number;
  controller: AbortController | undefined;
  /** The grace or retry timer. */
  timer: unknown;
  /** Answers that failed in a row. */
  attempts: number;
  /** A failure past its cooldown: asked again when wanted. */
  retryable: boolean;
  bytes: number;
  nodes: readonly string[];
}

interface AttachedSource {
  readonly opened: OpenedSource;
  readonly names: ReadonlyMap<string, string>;
  /** The node keys its coverage brought, held while it is attached. */
  readonly nodes: readonly string[];
}

/**
 * One store per viewer: data entries by `<source>|<system>|<level>|<tile>`, node entries by `<source>#<key>`, each a
 * keyed state so only what reads an entry is told. `want` sets the regions needed now; the store requests the tiles
 * they cover from the attached sources.
 */
export class DataStore extends ReactiveModel {
  readonly settings: StoreSettings;
  private readonly _clock: StoreClock;
  private readonly _sources = new Map<string, AttachedSource>();
  /** In least recently used order (a touch moves a slot to the end). */
  private readonly _slots = new Map<string, Slot>();
  private readonly _entries = new KeyedState<string, DataEntry>();
  private readonly _nodes = new KeyedState<string, NodeEntry>();
  private readonly _holders = new Map<string, number>();
  private readonly _queue: Slot[] = [];
  private _inFlight = 0;
  private _bytes = 0;
  private _counter = 0;
  private _order = 0;
  private readonly _reported = new Set<string>();
  private readonly _version = new State(0);
  private readonly _settled = new State(true);
  private readonly _overBudget = new State(false);
  private readonly _issues = new State<readonly Issue[]>([]);

  constructor(options: { settings?: Partial<StoreSettings>; clock?: StoreClock } = {}) {
    super();
    this.settings = { ...DefaultStoreSettings, ...options.settings };
    this._clock = options.clock ?? RealClock;
  }

  /** Moves when any entry does. */
  get version(): ReadonlyState<number> {
    return this._version;
  }

  /** No wanted tile pending and nothing wanted in flight. */
  get settled(): ReadonlyState<boolean> {
    return this._settled;
  }

  /** What is in view alone is over the byte bound (reported once, `seq.store.over_budget`). */
  get overBudget(): ReadonlyState<boolean> {
    return this._overBudget;
  }

  /** What the sources and the store reported, each cause once per source. */
  get issues(): ReadonlyState<readonly Issue[]> {
    return this._issues;
  }

  /** The bytes the data entries hold. */
  get bytes(): number {
    return this._bytes;
  }

  /** Requests in flight now. */
  get inFlight(): number {
    return this._inFlight;
  }

  // Sources

  attach(address: string, opened: OpenedSource): void {
    if (this._sources.has(address)) this.detach(address);
    const names = new Map<string, string>();
    for (const s of opened.coverage.systems) {
      names.set(s.name, s.name);
      for (const a of s.aliases ?? []) if (!names.has(a)) names.set(a, s.name);
    }
    batch(() => {
      const nodes = this._intakeNodes(address, opened.coverage.nodes);
      this._sources.set(address, { opened, names, nodes });
    });
  }

  /** Lets a source go: its requests cancelled, its entries and the holds they and its coverage had released. */
  detach(address: string): void {
    const s = this._sources.get(address);
    if (!s) return;
    batch(() => {
      for (const slot of [...this._slots.values()]) if (slot.region.source === address) this._drop(slot);
      for (const k of s.nodes) this._release(k);
      this._sources.delete(address);
      this._pump();
      this._updateSettled();
    });
    s.opened.dispose();
  }

  coverage(address: string): Coverage | undefined {
    return this._sources.get(address)?.opened.coverage;
  }

  // Entries

  entry(key: string): ReadonlyState<DataEntry | undefined> {
    return this._entries.key(key);
  }

  get(key: string): DataEntry | undefined {
    return this._entries.get(key);
  }

  node(source: string, key: string): ReadonlyState<NodeEntry | undefined> {
    return this._nodes.key(nodeKey(source, key));
  }

  getNode(source: string, key: string): NodeEntry | undefined {
    return this._nodes.get(nodeKey(source, key));
  }

  /** Holds a node (a placement shown, the selection, a highlight) until the returned release is called. */
  hold(source: string, key: string): () => void {
    const k = nodeKey(source, key);
    this._holders.set(k, (this._holders.get(k) ?? 0) + 1);
    let held = true;
    return () => {
      if (!held) return;
      held = false;
      batch(() => this._release(k));
    };
  }

  /** The tiles a want covers: the keys of the entries it reads, in position order; none for an unknown source or system. */
  keysOf(want: Want): string[] {
    const t = this._tiles(want);
    if (!t) return [];
    const out: string[] = [];
    for (let i = t.first; i <= t.last; i++) out.push(dataKey(want.source, t.system, t.level, i));
    return out;
  }

  // Requests

  /**
   * The regions needed now, all of them: each tile they cover is requested once (deduped), in the order of its best
   * priority, then of first want. What left the set is cancelled after the grace if in flight, dropped if not yet
   * asked for, and otherwise kept in the cache.
   */
  want(wants: readonly Want[]): void {
    batch(() => {
      const next = new Map<string, { rank: number; region: TileRegion; shapes: readonly Shape[] }>();
      for (const w of wants) {
        const t = this._tiles(w);
        if (!t) continue;
        const rank = Rank[w.priority ?? 'view'];
        for (let tile = t.first; tile <= t.last; tile++) {
          const key = dataKey(w.source, t.system, t.level, tile);
          const at = next.get(key);
          if (at) at.rank = Math.min(at.rank, rank);
          else {
            const { start, end } = tileRegion(tile, t.level, this.settings.tileBins, t.length);
            const region = { source: w.source, system: t.system, level: t.level, tile, start, end };
            next.set(key, { rank, region, shapes: t.shapes });
          }
        }
      }
      for (const slot of [...this._slots.values()])
        if (slot.rank !== Unwanted && !next.has(slot.key)) this._unwant(slot);
      for (const [key, w] of next) {
        const slot = this._slots.get(key);
        if (!slot) this._create(key, w.region, w.shapes, w.rank);
        else this._rewant(slot, w.rank);
      }
      this._pump();
      this._evict();
      this._updateSettled();
    });
  }

  /** Resolves `true` when settled, `false` after `timeoutMs` (if given) without settling (EC-125). */
  whenSettled(timeoutMs?: number): Promise<boolean> {
    if (this._settled.value) return Promise.resolve(true);
    return new Promise((resolve) => {
      let timer: unknown;
      const stop = this._settled.listen((s) => {
        if (!s) return;
        stop();
        if (timer !== undefined) this._clock.clearTimeout(timer);
        resolve(true);
      });
      if (timeoutMs !== undefined)
        timer = this._clock.setTimeout(() => {
          stop();
          resolve(false);
        }, timeoutMs);
    });
  }

  protected override onRetire(): void {
    batch(() => {
      for (const address of [...this._sources.keys()]) this.detach(address);
    });
  }

  // Tiles and slots

  private _tiles(w: Want) {
    const s = this._sources.get(w.source);
    const system = s?.names.get(w.system);
    if (!s || system === undefined) return undefined;
    const cov = s.opened.coverage;
    const level = cov.levels.includes(w.level) ? w.level : fitLevel(cov.levels, w.level);
    const length = cov.systems.find((x) => x.name === system)?.length;
    const end = length === undefined ? w.end : Math.min(w.end, length);
    const [first, last] = tileSpan(w.start, end, level, this.settings.tileBins);
    const shapes = cov.shapes.filter((x) => x.system === system).map((x) => x.shape);
    return { system, level, length, first, last, shapes };
  }

  private _setEntry(slot: Slot, entry: DataEntry): void {
    this._entries.set(slot.key, entry);
    this._version.set(entry.version);
  }

  private _next(): number {
    return ++this._counter;
  }

  private _create(key: string, region: TileRegion, shapes: readonly Shape[], rank: number): void {
    const slot: Slot = {
      key,
      region,
      shapes,
      rank,
      order: this._order++,
      phase: 'queued',
      token: 0,
      controller: undefined,
      timer: undefined,
      attempts: 0,
      retryable: false,
      bytes: 0,
      nodes: [],
    };
    this._slots.set(key, slot);
    this._queue.push(slot);
    this._setEntry(slot, { status: 'pending', version: this._next(), region });
  }

  private _touch(slot: Slot): void {
    this._slots.delete(slot.key);
    this._slots.set(slot.key, slot);
  }

  private _clearTimer(slot: Slot): void {
    if (slot.timer !== undefined) this._clock.clearTimeout(slot.timer);
    slot.timer = undefined;
  }

  private _rewant(slot: Slot, rank: number): void {
    const was = slot.rank;
    slot.rank = rank;
    if (was === Unwanted) slot.order = this._order++;
    this._touch(slot);
    // Back within its grace: it stays in flight.
    if (slot.phase === 'flight') this._clearTimer(slot);
    else if (slot.phase === 'done' && slot.retryable) this._enqueue(slot);
  }

  private _unwant(slot: Slot): void {
    slot.rank = Unwanted;
    switch (slot.phase) {
      case 'flight':
        // EC-173: a grace before cancelling, so a jitter back keeps the request; the oldest stale goes first.
        slot.timer = this._clock.setTimeout(() => this._cancel(slot), this.settings.graceMs);
        break;
      case 'queued':
      case 'backoff':
        if (this.get(slot.key)?.status === 'pending') this._drop(slot);
        else {
          this._dequeue(slot);
          this._clearTimer(slot);
          slot.phase = 'done';
          slot.retryable = true;
        }
        break;
    }
  }

  private _cancel(slot: Slot): void {
    slot.timer = undefined;
    if (slot.phase !== 'flight' || slot.rank !== Unwanted) return;
    batch(() => {
      if (this.get(slot.key)?.status === 'pending') this._drop(slot);
      else {
        this._letGo(slot);
        slot.phase = 'done';
        slot.retryable = true;
      }
      this._pump();
      this._updateSettled();
    });
  }

  /** Stops a request in flight: aborted, its answer dropped when it comes. */
  private _letGo(slot: Slot): void {
    if (slot.phase !== 'flight') return;
    slot.token++;
    slot.controller?.abort();
    slot.controller = undefined;
    this._inFlight--;
  }

  private _dequeue(slot: Slot): void {
    const i = this._queue.indexOf(slot);
    if (i >= 0) this._queue.splice(i, 1);
  }

  /** Removes a slot and its entry, releasing what it held. */
  private _drop(slot: Slot): void {
    this._letGo(slot);
    this._dequeue(slot);
    this._clearTimer(slot);
    this._slots.delete(slot.key);
    this._bytes -= slot.bytes;
    slot.bytes = 0;
    for (const k of slot.nodes) this._release(k);
    slot.nodes = [];
    this._entries.delete(slot.key);
    this._version.set(this._next());
  }

  private _enqueue(slot: Slot): void {
    slot.phase = 'queued';
    slot.retryable = false;
    this._queue.push(slot);
  }

  // Requests in flight

  private _pump(): void {
    const max = this.settings.maxInFlight;
    if (this._inFlight >= max || this._queue.length === 0) return;
    this._queue.sort((a, b) => a.rank - b.rank || a.order - b.order);
    while (this._inFlight < max && this._queue.length) this._start(this._queue.shift()!);
  }

  private _start(slot: Slot): void {
    const source = this._sources.get(slot.region.source);
    if (!source) return;
    slot.phase = 'flight';
    this._inFlight++;
    const controller = new AbortController();
    slot.controller = controller;
    const token = ++slot.token;
    const { system, start, end, level } = slot.region;
    const request: RegionRequest = { system, start, end, level, shapes: slot.shapes };
    let reading: Promise<RegionAnswer>;
    try {
      reading = source.opened.read(request, controller.signal);
    } catch (e) {
      reading = Promise.reject(e);
    }
    reading.then(
      (answer) => this._onAnswer(slot, token, answer),
      (error) => this._onError(slot, token, error),
    );
  }

  private _current(slot: Slot, token: number): boolean {
    return slot.token === token && slot.phase === 'flight' && this._slots.get(slot.key) === slot;
  }

  private _onAnswer(slot: Slot, token: number, answer: RegionAnswer): void {
    if (!this._current(slot, token)) return;
    batch(() => {
      this._settle(slot);
      if (!answer.ok) this._failure(slot, answer.kind, answer.issue);
      else {
        const bad = checkTables(answer.tables);
        if (bad) {
          const issue = issueOf(SourceIssues, 'seq.source.malformed', `The source answered a region malformed: ${bad}`);
          this._failure(slot, 'malformed', issue);
        } else this._loaded(slot, answer.tables, answer.nodes);
      }
      this._pump();
      this._evict();
      this._updateSettled();
    });
  }

  private _onError(slot: Slot, token: number, error: unknown): void {
    // An abort the store asked for was let go already (`_current` is false): never an error state.
    if (!this._current(slot, token)) return;
    batch(() => {
      this._settle(slot);
      const message = isAbortError(error) ? 'The read was aborted by the source' : 'The read failed';
      this._failure(slot, 'failed', issueOf(SourceIssues, 'seq.source.failed', message));
      this._pump();
      this._updateSettled();
    });
  }

  /** An answer came: the request is no longer in flight. */
  private _settle(slot: Slot): void {
    this._clearTimer(slot);
    slot.controller = undefined;
    slot.phase = 'done';
    this._inFlight--;
  }

  private _loaded(slot: Slot, tables: Tables, nodes: readonly DynamicNode[] | undefined): void {
    slot.attempts = 0;
    for (const k of slot.nodes) this._release(k);
    slot.nodes = this._intakeNodes(slot.region.source, nodes);
    this._bytes -= slot.bytes;
    const empty = isEmptyAnswer(tables, nodes);
    slot.bytes = empty ? 0 : answerBytes(tables, nodes);
    this._bytes += slot.bytes;
    this._touch(slot);
    const version = this._next();
    const region = slot.region;
    if (empty) this._setEntry(slot, { status: 'empty', version, region });
    else this._setEntry(slot, { status: 'loaded', version, region, tables, nodes: slot.nodes, bytes: slot.bytes });
    if (slot.rank === Unwanted) this._clearTimer(slot);
  }

  /** EC-174: `failed` and `malformed` retried with backoff, then a cooldown; the others never; each reported once. */
  private _failure(slot: Slot, kind: FailureKind, issue: Issue): void {
    const source = slot.region.source;
    this._report(source, { ...issue, address: issue.address ?? source });
    const retry = this.settings.retry;
    const fail = () =>
      this._setEntry(slot, { status: 'failed', version: this._next(), region: slot.region, kind, issue });
    if (kind !== 'failed' && kind !== 'malformed') {
      fail();
      return;
    }
    slot.attempts++;
    if (slot.attempts < retry.attempts) {
      slot.phase = 'backoff';
      const delay = retry.baseMs * retry.factor ** (slot.attempts - 1);
      slot.timer = this._clock.setTimeout(() => this._retry(slot), delay);
      return;
    }
    slot.attempts = 0;
    fail();
    slot.phase = 'cooldown';
    slot.timer = this._clock.setTimeout(() => this._cooled(slot), retry.cooldownMs);
  }

  private _retry(slot: Slot): void {
    slot.timer = undefined;
    if (this._slots.get(slot.key) !== slot || slot.phase !== 'backoff') return;
    batch(() => {
      if (slot.rank !== Unwanted) this._enqueue(slot);
      else if (this.get(slot.key)?.status === 'pending') this._drop(slot);
      else {
        slot.phase = 'done';
        slot.retryable = true;
      }
      this._pump();
      this._updateSettled();
    });
  }

  private _cooled(slot: Slot): void {
    slot.timer = undefined;
    if (this._slots.get(slot.key) !== slot || slot.phase !== 'cooldown') return;
    slot.phase = 'done';
    slot.retryable = true;
    if (slot.rank === Unwanted) return;
    batch(() => {
      this._enqueue(slot);
      this._pump();
    });
  }

  private _report(source: string, issue: Issue): void {
    const id = `${source}\u0000${issue.code}`;
    if (this._reported.has(id)) return;
    this._reported.add(id);
    this._issues.set([...this._issues.value, issue]);
  }

  private _updateSettled(): void {
    let settled = true;
    for (const slot of this._slots.values()) {
      if (slot.rank === Unwanted) continue;
      if (slot.phase === 'queued' || slot.phase === 'flight' || slot.phase === 'backoff') {
        settled = false;
        break;
      }
    }
    this._settled.set(settled);
  }

  // Eviction

  /**
   * Least recently used first among what is not wanted, then prefetched, then in the margin; never what is in view or
   * a request on its way. When what is left is still over the bound, the store is over budget (reported once a time).
   */
  private _evict(): void {
    const max = this.settings.maxBytes;
    while (this._bytes > max) {
      let victim: Slot | undefined;
      for (const slot of this._slots.values()) {
        if (slot.rank === Rank.view || slot.phase === 'queued' || slot.phase === 'flight' || slot.phase === 'backoff')
          continue;
        if (slot.bytes === 0) continue;
        if (!victim || slot.rank > victim.rank) victim = slot;
        if (victim.rank === Unwanted) break;
      }
      if (!victim) break;
      this._drop(victim);
    }
    const over = this._bytes > max;
    if (over && !this._overBudget.value) {
      this._issues.set([
        ...this._issues.value,
        issueOf(StoreIssues, 'seq.store.over_budget', 'What is in view is over the store’s byte bound'),
      ]);
    }
    this._overBudget.set(over);
  }

  // Nodes

  /** Holds the nodes an answer (or a coverage) brought; the same key with other content keeps the held node. */
  private _intakeNodes(source: string, nodes: readonly DynamicNode[] | undefined): readonly string[] {
    if (!nodes?.length) return [];
    const keys: string[] = [];
    const seen = new Set<string>();
    for (const n of nodes) {
      const k = nodeKey(source, n.key);
      if (seen.has(k)) continue;
      seen.add(k);
      keys.push(k);
      this._holders.set(k, (this._holders.get(k) ?? 0) + 1);
      const held = this._nodes.get(k);
      if (!held) {
        const version = this._next();
        this._nodes.set(k, { source, key: n.key, node: n.node, version });
        this._version.set(version);
      } else if (!nodesEqual(held.node, n.node)) {
        const issue = issueOf(StoreIssues, 'seq.store.node_conflict', 'A node came again with other content: kept', {
          address: source,
        });
        this._report(source, issue);
      }
    }
    return keys;
  }

  /** A node with no holder left goes, quietly. */
  private _release(k: string): void {
    const n = (this._holders.get(k) ?? 0) - 1;
    if (n > 0) {
      this._holders.set(k, n);
      return;
    }
    this._holders.delete(k);
    if (this._nodes.has(k)) {
      this._nodes.delete(k);
      this._version.set(this._next());
    }
  }
}
