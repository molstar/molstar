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
import { checkCoverage } from './seq-format.js';
import { DefaultSummarySettings, type SummarySettings, withSummaries } from './summaries.js';
import { answerBytes, checkAnswer, DefaultMaxRows, isEmptyAnswer } from './tables.js';

// The data store (_spec/README.md "Data", Store): entries by key in keyed states, versioned; the want sets of its
// consumers merged; requests deduped by tile, ordered and bounded in flight, cancelled after a grace (at once when the
// view waits), retried by policy; a byte-bounded cache that never drops what is wanted; nodes held by what carried
// them.

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

/**
 * A data entry: pending is never drawn or hit as empty (EC-169); failed carries its issue; evicted is wanted but
 * dropped for the byte bound, not asked again while it stays wanted.
 */
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
    }
  | { readonly status: 'evicted'; readonly version: number; readonly region: TileRegion };

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
  /** Tiles one want may cover; more is an issue and the want is not requested (work bounded before it starts). */
  readonly maxTilesPerWant: number;
  /** Rows a list or density table may hold in one answer. */
  readonly maxRows: number;
  /** Issues kept; then one `seq.store.limit_exceeded`. */
  readonly maxIssues: number;
  /** How long a request that left the wanted set stays in flight before it is cancelled (EC-173). */
  readonly graceMs: number;
  /**
   * A `failed` answer: tried `attempts` times with backoff `baseMs · factor^n`, then a cooldown (a `malformed` one goes
   * to the cooldown at once).
   */
  readonly retry: {
    readonly attempts: number;
    readonly baseMs: number;
    readonly factor: number;
    readonly cooldownMs: number;
  };
  /** The value pyramid given a source that serves level 1 only (`withSummaries`); `false`: none. */
  readonly summaries: SummarySettings | false;
}

export const DefaultStoreSettings: StoreSettings = {
  tileBins: DefaultTileBins,
  maxBytes: 64 * 1024 * 1024,
  maxInFlight: 6,
  maxTilesPerWant: 256,
  maxRows: DefaultMaxRows,
  maxIssues: 100,
  graceMs: 300,
  retry: { attempts: 3, baseMs: 500, factor: 2, cooldownMs: 30_000 },
  summaries: DefaultSummarySettings,
};

/** A consumer's want set (the live view, an export): merged with the others by best priority. */
export interface Wanter {
  /** The regions this consumer needs now, all of them. */
  want(wants: readonly Want[]): void;
  /** `true` when none of this consumer's tiles is pending; `false` after `timeoutMs` (if given) without (EC-125). */
  whenSettled(timeoutMs?: number): Promise<boolean>;
  /** The keys of the entries this consumer wants. */
  readonly keys: ReadonlySet<string>;
  /** Withdraws its wants. */
  dispose(): void;
}

/**
 * `queued`, `flight`: being asked; `backoff`: waiting to be asked again after a failure; `cooldown`: failed, waiting
 * longer; `done`: answered (or let go: `retryable` says it is asked again when wanted); `evicted`: wanted, dropped for
 * the bound.
 */
type Phase = 'queued' | 'flight' | 'backoff' | 'cooldown' | 'done' | 'evicted';

interface PlannedTile {
  rank: number;
  readonly region: TileRegion;
  readonly shapes: readonly Shape[];
}

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
  request: RegionRequest | undefined;
  /** The grace, backoff or cooldown timer. */
  timer: unknown;
  /** When it left the wanted set while in flight (EC-173's order), else -1. */
  staleSince: number;
  /** `failed` answers in a row; kept while it is away and back. */
  attempts: number;
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

interface Reported {
  readonly source: string;
  readonly id: string;
  readonly issue: Issue;
}

/** The store-wide reports' source. */
const StoreWide = '';

/**
 * One store per viewer: data entries by `<source>|<system>|<level>|<tile>`, node entries by `<source>#<key>`, each a
 * keyed state so only what reads an entry is told. Each consumer's `Wanter` (the live view's is `want`) sets the
 * regions it needs now; the store requests the tiles they cover from the attached sources.
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
  private readonly _wanters = new Set<WantSet>();
  private readonly _live: WantSet;
  private _inFlight = 0;
  private _bytes = 0;
  private _counter = 0;
  private _order = 0;
  private readonly _reported: Reported[] = [];
  private _capped = false;
  private readonly _version = new State(0);
  private readonly _nodeVersion = new State(0);
  private readonly _settled = new State(true);
  private readonly _overBudget = new State(false);
  private readonly _issues = new State<readonly Issue[]>([]);

  /** The clock the store's timers run on (a caller sharing the store keeps time on it). */
  get clock(): StoreClock {
    return this._clock;
  }

  constructor(options: { settings?: Partial<StoreSettings>; clock?: StoreClock } = {}) {
    super();
    this.settings = { ...DefaultStoreSettings, ...options.settings };
    this._clock = options.clock ?? RealClock;
    this._live = this.wanter() as WantSet;
  }

  /** Moves when any entry does. */
  get version(): ReadonlyState<number> {
    return this._version;
  }

  /** Moves when a node entry comes or goes (the overlay's version, composition.md's tier 2). */
  get nodeVersion(): ReadonlyState<number> {
    return this._nodeVersion;
  }

  /** No wanted entry pending (a failed one waiting to be asked again is settled). */
  get settled(): ReadonlyState<boolean> {
    return this._settled;
  }

  /** What is in view alone is over the byte bound (reported once each time, `seq.store.over_budget`). */
  get overBudget(): ReadonlyState<boolean> {
    return this._overBudget;
  }

  /** What the sources and the store reported: each cause once per source, at most `maxIssues`. */
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

  /**
   * Attaches an opened source under its address. Its coverage is checked first: a malformed one is an issue and the
   * source stays unavailable (`false`); a retired store attaches nothing (`false`).
   */
  attach(address: string, source: OpenedSource): boolean {
    // a late opening after the store is retired is the caller's to dispose
    if (this.isRetired) return false;
    if (this._sources.has(address)) this.detach(address);
    const bad = checkCoverage(source.coverage);
    if (bad) {
      this._run(() => this._report(address, { ...bad, address }));
      return false;
    }
    // a level-1 source's values get a pyramid, so a values layer zoomed out reads a coarse level
    const summaries = this.settings.summaries;
    const opened = summaries ? withSummaries(source, summaries) : source;
    const names = new Map<string, string>();
    for (const s of opened.coverage.systems) {
      names.set(s.name, s.name);
      for (const a of s.aliases ?? []) if (!names.has(a)) names.set(a, s.name);
    }
    this._run(() => {
      const nodes = this._intakeNodes(address, opened.coverage.nodes);
      this._sources.set(address, { opened, names, nodes });
    });
    return true;
  }

  /**
   * Lets a source go: its requests cancelled, its entries and the holds they and its coverage had released, its
   * reports cleared. Wants of it stay and wait for a source attached again under the address.
   */
  detach(address: string): void {
    const s = this._sources.get(address);
    if (!s) return;
    this._run(() => {
      for (const slot of [...this._slots.values()]) if (slot.region.source === address) this._drop(slot);
      for (const k of s.nodes) this._release(k);
      this._sources.delete(address);
      this._clearReports(address);
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

  /** Every node entry held now, by `<source>#<key>` (live: read it, never keep it). */
  get nodes(): ReadonlyMap<string, NodeEntry> {
    return this._nodes.values;
  }

  /** Holds a node (a placement shown, the selection, a highlight) until the returned release is called. */
  hold(source: string, key: string): () => void {
    const k = nodeKey(source, key);
    this._holders.set(k, (this._holders.get(k) ?? 0) + 1);
    let held = true;
    return () => {
      if (!held) return;
      held = false;
      this._run(() => this._release(k));
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

  /** A new consumer's want set (an export's, beside the live view's), merged with the others by best priority. */
  wanter(): Wanter {
    const set = new WantSet(this);
    this._wanters.add(set);
    return set;
  }

  /** The live view's want set (`wanter()` for another consumer). */
  want(wants: readonly Want[]): void {
    this._live.want(wants);
  }

  /** Resolves `true` when settled, `false` after `timeoutMs` (if given) without settling (EC-125). */
  whenSettled(timeoutMs?: number): Promise<boolean> {
    return this._waitFor(() => this._settled.value, timeoutMs);
  }

  protected override onRetire(): void {
    batch(() => {
      for (const w of [...this._wanters]) w.dispose();
      for (const address of [...this._sources.keys()]) this.detach(address);
    });
  }

  /** @internal Plans a consumer's wants into tiles. */
  _plan(wants: readonly Want[]): Map<string, PlannedTile> {
    const next = new Map<string, PlannedTile>();
    for (const w of wants) {
      const t = this._tiles(w);
      if (!t || t.last < t.first) continue;
      if (t.last - t.first + 1 > this.settings.maxTilesPerWant) {
        const message = `A wanted region covers ${t.last - t.first + 1} tiles, more than ${this.settings.maxTilesPerWant}`;
        this._run(() =>
          this._report(w.source, issueOf(StoreIssues, 'seq.store.too_many_tiles', message, { address: w.source })),
        );
        continue;
      }
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
    return next;
  }

  /**
   * @internal The want sets merged and applied: each tile requested once, in the order of its best priority, then of
   * first want. What left every set is cancelled after the grace if in flight, dropped if not yet asked for or
   * evicted, and otherwise kept (in the cache, or in its backoff).
   */
  _apply(): void {
    this._run(() => {
      const merged = new Map<string, PlannedTile>();
      for (const set of this._wanters)
        for (const [key, t] of set.tiles) {
          const at = merged.get(key);
          if (!at) merged.set(key, { ...t });
          else at.rank = Math.min(at.rank, t.rank);
        }
      for (const slot of [...this._slots.values()])
        if (slot.rank !== Unwanted && !merged.has(slot.key)) this._unwant(slot);
      for (const [key, w] of merged) {
        const slot = this._slots.get(key);
        if (!slot) this._create(key, w.region, w.shapes, w.rank);
        else this._rewant(slot, w.rank);
      }
    });
  }

  /** @internal */
  _removeSet(set: WantSet): void {
    if (this._wanters.delete(set)) this._apply();
  }

  /** @internal Resolves `true` once `check` holds after a change, `false` after `timeoutMs`. */
  _waitFor(check: () => boolean, timeoutMs?: number): Promise<boolean> {
    if (check()) return Promise.resolve(true);
    return new Promise((resolve) => {
      let timer: unknown;
      const stop = this._version.listen(() => {
        if (!check()) return;
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

  /** Runs a change in one batch, then starts what may start, keeps the bound and the settled state. */
  private _run(fn: () => void): void {
    batch(() => {
      fn();
      this._pump();
      this._evict();
      this._updateSettled();
    });
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
      request: undefined,
      timer: undefined,
      staleSince: -1,
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
    if (slot.phase === 'flight' && slot.staleSince >= 0) {
      // Back within its grace: it stays in flight.
      this._clearTimer(slot);
      slot.staleSince = -1;
    } else if (slot.phase === 'done' && slot.retryable) this._enqueue(slot);
    // A slot in backoff or cooldown keeps its timer and its attempts; an evicted one waits to leave the set.
  }

  private _unwant(slot: Slot): void {
    slot.rank = Unwanted;
    switch (slot.phase) {
      case 'flight':
        // EC-173: a grace before cancelling, so a jitter back keeps the request; the oldest stale goes first.
        slot.staleSince = this._clock.now();
        slot.timer = this._clock.setTimeout(() => this._graceOver(slot), this.settings.graceMs);
        break;
      case 'queued':
        if (this.get(slot.key)?.status === 'pending') this._drop(slot);
        else {
          this._dequeue(slot);
          slot.phase = 'done';
          slot.retryable = true;
        }
        break;
      case 'evicted':
        this._drop(slot);
        break;
    }
  }

  private _graceOver(slot: Slot): void {
    slot.timer = undefined;
    if (slot.phase !== 'flight' || slot.rank !== Unwanted) return;
    this._run(() => this._cancel(slot));
  }

  /** A stale request let go: a pending entry goes with it, a failed one stays, asked again when wanted. */
  private _cancel(slot: Slot): void {
    if (this.get(slot.key)?.status === 'pending') this._drop(slot);
    else {
      this._clearTimer(slot);
      this._letGo(slot);
      slot.phase = 'done';
      slot.retryable = true;
    }
  }

  /** Stops a request in flight: aborted, its answer dropped when it comes. */
  private _letGo(slot: Slot): void {
    if (slot.phase !== 'flight') return;
    slot.token++;
    slot.staleSince = -1;
    const c = slot.controller;
    slot.controller = undefined;
    this._inFlight--;
    c?.abort();
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
    this._free(slot);
    this._entries.delete(slot.key);
    this._version.set(this._next());
  }

  private _free(slot: Slot): void {
    this._bytes -= slot.bytes;
    slot.bytes = 0;
    const nodes = slot.nodes;
    slot.nodes = [];
    for (const k of nodes) this._release(k);
  }

  private _enqueue(slot: Slot): void {
    slot.phase = 'queued';
    slot.retryable = false;
    this._queue.push(slot);
  }

  // Requests in flight

  /** The oldest stale request in flight. */
  private _oldestStale(): Slot | undefined {
    let out: Slot | undefined;
    for (const s of this._slots.values())
      if (s.phase === 'flight' && s.staleSince >= 0 && (!out || s.staleSince < out.staleSince)) out = s;
    return out;
  }

  private _pump(): void {
    const max = this.settings.maxInFlight;
    if (this._queue.length === 0) return;
    this._queue.sort((a, b) => a.rank - b.rank || a.order - b.order);
    while (this._queue.length) {
      if (this._inFlight >= max) {
        // EC-173: the view waits, so the oldest stale request goes now, not after its grace.
        const stale = this._queue[0].rank === Rank.view ? this._oldestStale() : undefined;
        if (!stale) break;
        this._cancel(stale);
        continue;
      }
      this._start(this._queue.shift()!);
    }
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
    slot.request = request;
    let reading: unknown;
    try {
      reading = source.opened.read(request, controller.signal);
    } catch (e) {
      reading = Promise.reject(e);
    }
    if (typeof (reading as Promise<unknown> | null)?.then !== 'function') {
      this._answered(slot, () => this._malformed(slot, 'the read did not return a promise'));
      return;
    }
    (reading as Promise<RegionAnswer>).then(
      (answer) => {
        if (this._current(slot, token)) this._run(() => this._answered(slot, () => this._onAnswer(slot, answer)));
      },
      (error) => {
        // An abort the store asked for was let go already (`_current` is false): never an error state.
        if (this._current(slot, token)) this._run(() => this._answered(slot, () => this._onError(slot, error)));
      },
    );
  }

  private _current(slot: Slot, token: number): boolean {
    return slot.token === token && slot.phase === 'flight' && this._slots.get(slot.key) === slot;
  }

  /** An answer came: the request is no longer in flight; what handling it throws makes it malformed. */
  private _answered(slot: Slot, handle: () => void): void {
    this._clearTimer(slot);
    slot.controller = undefined;
    slot.staleSince = -1;
    slot.phase = 'done';
    this._inFlight--;
    try {
      handle();
    } catch {
      this._malformed(slot, 'the answer could not be taken in');
    }
  }

  private _malformed(slot: Slot, why: string): void {
    const message = `The source answered a region malformed: ${why}`;
    this._failure(slot, 'malformed', issueOf(SourceIssues, 'seq.source.malformed', message));
  }

  private _onAnswer(slot: Slot, answer: RegionAnswer): void {
    const bad = checkAnswer(answer, slot.request!, this.settings.maxRows);
    if (bad) this._malformed(slot, bad);
    else if (!answer.ok) this._failure(slot, answer.kind, answer.issue);
    else this._loaded(slot, answer.tables, answer.nodes);
  }

  private _onError(slot: Slot, error: unknown): void {
    // An abort the store did not ask for is the source's failure.
    const message = isAbortError(error) ? 'The read was aborted by the source' : 'The read failed';
    this._failure(slot, 'failed', issueOf(SourceIssues, 'seq.source.failed', message));
  }

  private _loaded(slot: Slot, tables: Tables, nodes: readonly DynamicNode[] | undefined): void {
    slot.attempts = 0;
    // The new nodes first, then the old released: a node in both is never dropped between.
    const held = this._intakeNodes(slot.region.source, nodes);
    this._free(slot);
    slot.nodes = held;
    const empty = isEmptyAnswer(tables, nodes);
    slot.bytes = empty ? 0 : answerBytes(tables, nodes);
    this._bytes += slot.bytes;
    this._touch(slot);
    const version = this._next();
    const region = slot.region;
    if (empty) this._setEntry(slot, { status: 'empty', version, region });
    else this._setEntry(slot, { status: 'loaded', version, region, tables, nodes: held, bytes: slot.bytes });
  }

  /**
   * EC-174: `failed` is retried with backoff (the entry as it was), then fails with a cooldown; `malformed` fails with
   * the cooldown at once; `not-served` and `too-large` are never retried. Each cause is reported once.
   */
  private _failure(slot: Slot, kind: FailureKind, issue: Issue): void {
    const source = slot.region.source;
    this._report(source, { ...issue, address: issue.address ?? source });
    const retry = this.settings.retry;
    if (kind === 'failed' && slot.attempts + 1 < retry.attempts) {
      slot.attempts++;
      slot.phase = 'backoff';
      const delay = retry.baseMs * retry.factor ** (slot.attempts - 1);
      slot.timer = this._clock.setTimeout(() => this._backoffOver(slot), delay);
      return;
    }
    slot.attempts = 0;
    this._setEntry(slot, { status: 'failed', version: this._next(), region: slot.region, kind, issue });
    if (kind === 'failed' || kind === 'malformed') {
      slot.phase = 'cooldown';
      slot.timer = this._clock.setTimeout(() => this._backoffOver(slot), retry.cooldownMs);
    } else slot.phase = 'done';
  }

  /** A backoff or a cooldown is over: asked again if wanted, else when it is wanted again. */
  private _backoffOver(slot: Slot): void {
    slot.timer = undefined;
    if (this._slots.get(slot.key) !== slot || (slot.phase !== 'backoff' && slot.phase !== 'cooldown')) return;
    this._run(() => {
      slot.phase = 'done';
      slot.retryable = true;
      if (slot.rank !== Unwanted) this._enqueue(slot);
    });
  }

  // Reports

  private _report(source: string, issue: Issue): void {
    const id = `${source}\u0000${issue.code}`;
    if (this._reported.some((r) => r.id === id)) return;
    if (this._reported.length >= this.settings.maxIssues) {
      if (!this._capped) {
        this._capped = true;
        this._publishIssues();
      }
      return;
    }
    this._reported.push({ source, id, issue });
    this._publishIssues();
  }

  private _clearReports(source: string): void {
    const before = this._reported.length;
    for (let i = this._reported.length - 1; i >= 0; i--)
      if (this._reported[i].source === source) this._reported.splice(i, 1);
    if (this._reported.length === before) return;
    this._capped = false;
    this._publishIssues();
  }

  private _publishIssues(): void {
    const list = this._reported.map((r) => r.issue);
    if (this._capped)
      list.push(issueOf(StoreIssues, 'seq.store.limit_exceeded', `More than ${this.settings.maxIssues} issues`));
    this._issues.set(list);
  }

  private _updateSettled(): void {
    let settled = true;
    for (const slot of this._slots.values()) {
      if (slot.rank !== Unwanted && this._entries.get(slot.key)?.status === 'pending') {
        settled = false;
        break;
      }
    }
    this._settled.set(settled);
  }

  // Eviction

  /**
   * Least recently used first among what is not wanted; still over, the wanted prefetch then margin tiles become
   * `evicted`. Never the view, never a request on its way. Still over: over budget, reported once each time.
   */
  private _evict(): void {
    const max = this.settings.maxBytes;
    for (const rank of [Unwanted, Rank.prefetch, Rank.margin]) {
      if (this._bytes <= max) break;
      for (const slot of [...this._slots.values()]) {
        if (this._bytes <= max) break;
        if (slot.rank !== rank || slot.bytes === 0 || slot.phase !== 'done') continue;
        if (rank === Unwanted) this._drop(slot);
        else {
          this._free(slot);
          slot.phase = 'evicted';
          this._setEntry(slot, { status: 'evicted', version: this._next(), region: slot.region });
        }
      }
    }
    const over = this._bytes > max;
    if (over && !this._overBudget.value) {
      const message = 'What is in view is over the store’s byte bound';
      this._report(StoreWide, issueOf(StoreIssues, 'seq.store.over_budget', message));
    } else if (!over && this._overBudget.value) this._clearReports(StoreWide);
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
        this._nodeVersion.set(version);
      } else if (!nodesEqual(held.node, n.node)) {
        const message = 'A node came again with other content: kept';
        this._report(source, issueOf(StoreIssues, 'seq.store.node_conflict', message, { address: source }));
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
      const version = this._next();
      this._version.set(version);
      this._nodeVersion.set(version);
    }
  }
}

class WantSet implements Wanter {
  tiles = new Map<string, PlannedTile>();
  private _disposed = false;

  constructor(private readonly _store: DataStore) {}

  get keys(): ReadonlySet<string> {
    return new Set(this.tiles.keys());
  }

  want(wants: readonly Want[]): void {
    if (this._disposed) return;
    this.tiles = this._store._plan(wants);
    this._store._apply();
  }

  whenSettled(timeoutMs?: number): Promise<boolean> {
    return this._store._waitFor(() => {
      for (const k of this.tiles.keys()) if (this._store.get(k)?.status === 'pending') return false;
      return true;
    }, timeoutMs);
  }

  dispose(): void {
    if (this._disposed) return;
    this._disposed = true;
    this.tiles = new Map();
    this._store._removeSet(this);
  }
}
