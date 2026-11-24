/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { batch, KeyedState, ReactiveModel, type ReadonlyState, State } from '@molstar/core/reactive';
import { type Issue, issueOf } from '@molstar/core/repr/tree';
import { channelSource } from './channels.js';
import { SourceIssues, StoreIssues } from './issues.js';
import { dataKey, nodeKey } from './keys.js';
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
import { StoreNodes } from './store-nodes.js';
import { StoreReports, StoreWide } from './store-reports.js';
import { StoreSources } from './store-sources.js';
import {
  type DataEntry,
  DefaultStoreSettings,
  type NodeEntry,
  Rank,
  RealClock,
  type StoreClock,
  type StoreSettings,
  type TileRegion,
  Unwanted,
  type Want,
  type Wanter,
} from './store-types.js';
import {
  eachTile,
  type LayerLoad,
  layerLoad,
  knownEmpty,
  type PlannedTile,
  planWants,
  tileRequest,
  WantSet,
  type WantTile,
  wantTiles,
} from './store-wants.js';
import { answerBytes, checkAnswer, isEmptyAnswer } from './tables.js';

export {
  type DataEntry,
  DefaultStoreSettings,
  type EntryStatus,
  type NodeEntry,
  type Priority,
  RealClock,
  type StoreClock,
  type StoreSettings,
  type TileRegion,
  type Want,
  type Wanter,
} from './store-types.js';

// The data store (_spec/README.md "Data", Store): entries by key in keyed states, versioned; the want sets of its
// consumers merged; requests deduped by tile, ordered and bounded in flight, cancelled after a grace (at once when the
// view waits), retried by policy; a byte-bounded cache that never drops what is wanted; nodes held by what carried
// them.

/**
 * `queued`, `flight`: being asked; `backoff`: waiting to be asked again after a failure; `cooldown`: failed, waiting
 * longer; `done`: answered (or let go: `retryable` says it is asked again when wanted); `evicted`: wanted, dropped for
 * the bound.
 */
type Phase = 'queued' | 'flight' | 'backoff' | 'cooldown' | 'done' | 'evicted';

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

/**
 * One store per viewer: data entries by `<source>|<system>|<level>|<tile>`, node entries by `<source>#<key>`, each a
 * keyed state so only what reads an entry is told. Each consumer's `Wanter` (the live view's is `want`) sets the
 * regions it needs now; the store requests the tiles they cover from the attached sources.
 */
export class DataStore extends ReactiveModel {
  readonly settings: StoreSettings;
  private readonly _clock: StoreClock;
  private readonly _sources: StoreSources;
  /** In least recently used order (a touch moves a slot to the end). */
  private readonly _slots = new Map<string, Slot>();
  private readonly _entries = new KeyedState<string, DataEntry>();
  private readonly _nodes: StoreNodes;
  private readonly _queue: Slot[] = [];
  private readonly _wanters = new Set<WantSet>();
  private readonly _live: WantSet;
  private _inFlight = 0;
  private _bytes = 0;
  private _counter = 0;
  private _order = 0;
  private readonly _version = new State(0);
  private readonly _nodeVersion = new State(0);
  private readonly _settled = new State(true);
  /** The summaries' builds the settled state waits on. */
  private readonly _awaited = new WeakSet<Promise<unknown>>();
  private readonly _overBudget = new State(false);
  private readonly _issues = new State<readonly Issue[]>([]);
  private readonly _reports: StoreReports;
  /** The layers whose wants together were past `maxTilesPerLayer` (`tiles`) or `maxBytesPerLayer` (`cost`). */
  private readonly _refused = new Map<string, 'tiles' | 'cost'>();
  private _refusals = 0;

  /** The clock the store's timers run on (a caller sharing the store keeps time on it). */
  get clock(): StoreClock {
    return this._clock;
  }

  constructor(options: { settings?: Partial<StoreSettings>; clock?: StoreClock } = {}) {
    super();
    this.settings = { ...DefaultStoreSettings, ...options.settings };
    this._clock = options.clock ?? RealClock;
    this._sources = new StoreSources(this.settings);
    this._reports = new StoreReports(this.settings.maxIssues, (list) => this._issues.set(list));
    const conflict = (source: string) => {
      const message = 'A node came again with other content: kept';
      this._reports.report(source, issueOf(StoreIssues, 'seq.store.node_conflict', message, { address: source }));
    };
    const bump = (version: number) => {
      this._version.set(version);
      this._nodeVersion.set(version);
    };
    const takeStale = (k: string) => this._sources.takeStale(k);
    this._nodes = new StoreNodes({ bump, conflict, takeStale }, () => this._next());
    this._live = this.wanter() as WantSet;
  }

  /** Moves when any entry does. */
  get version(): ReadonlyState<number> {
    return this._version;
  }

  /** Moves when a node entry comes or goes (the overlay's version, _spec/README.md "Overlay"). */
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

  /** The bytes the attached sources' value pyramids hold: counted against `maxBytes` with the entries'. */
  get summaryBytes(): number {
    return this._sources.summaryBytes;
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
    if (this._sources.attached(address)) this.detach(address);
    this._unopened.delete(address);
    const prepared = this._prepare(address, source);
    if (!prepared) return false;
    this._run(() => {
      const nodes = this._nodes.intake(address, prepared.opened.coverage.nodes);
      this._sources.set(address, { ...prepared, nodes });
    });
    return true;
  }

  /**
   * The attached sources by their addresses, as attached (before their summaries): what `findNames` asks (a locus box's
   * fallback to a name, EC-228).
   */
  *named(): IterableIterator<{ readonly address: string; readonly source: OpenedSource }> {
    for (const [address, a] of this._sources.all()) yield { address, source: a.raw };
  }

  /**
   * A new reading of an attached source (a host's edit, _spec/README.md "Data", Store, Reattach): nodes of a key it
   * still brings keep their identity and take its content (version moved); keys it no longer brings go, whoever holds
   * them (the reading says they are gone); every tile of it and its channels is dropped and asked again while wanted.
   * A malformed coverage is an issue and the old reading stays (`false`: the caller disposes `source`). Not attached:
   * `attach`.
   */
  reattach(address: string, source: OpenedSource): boolean {
    if (this.isRetired) return false;
    const old = this._sources.attached(address);
    if (!old) return this.attach(address, source);
    const prepared = this._prepare(address, source);
    if (!prepared) return false;
    this._run(() => {
      // the new nodes first (replacing content by key), then the old holds let go, so a kept key is never dropped
      const nodes = this._nodes.intake(address, prepared.opened.coverage.nodes, true);
      const kept = new Set(nodes);
      // what the old reading's tiles brought keeps its identity while held: replaced when served again, gone when the
      // tile that brought it is read without it
      for (const slot of this._slots.values())
        if (channelSource(slot.region.source) === address)
          for (const k of slot.nodes) if (!kept.has(k)) this._sources.markStale(address, k, slot.key);
      this._dropTiles(address);
      for (const k of old.nodes) this._nodes.release(k);
      for (const k of old.nodes) if (!kept.has(k)) this._nodes.remove(k);
      this._sources.set(address, { ...prepared, nodes });
      this._reports.clear(address);
      for (const set of this._wanters) set.tiles = this._plan(set.wants);
    });
    old.opened.dispose();
    this._apply();
    return true;
  }

  /** A source checked and summarized; undefined (said) when its coverage is malformed. */
  private _prepare(address: string, source: OpenedSource) {
    const prepared = this._sources.prepare(source);
    if ('opened' in prepared) return prepared;
    this._run(() => this._reports.report(address, { ...prepared, address }));
    return undefined;
  }

  /** Drops the tiles of a source and of its channels, and the channels. */
  private _dropTiles(address: string): void {
    for (const slot of [...this._slots.values()]) if (channelSource(slot.region.source) === address) this._drop(slot);
    this._sources.dropChannels(address);
  }

  /**
   * Lets a source go: its requests cancelled, its entries and the holds they and its coverage had released, its
   * reports cleared. Nothing of it stays pending: its tiles a want set still lists get no entry until it is wanted again
   * with a source attached under the address.
   */
  detach(address: string): void {
    this._unopened.delete(address);
    const s = this._sources.attached(address);
    if (!s) return;
    this._run(() => {
      this._dropTiles(address);
      for (const k of s.nodes) this._nodes.release(k);
      this._sources.delete(address);
      this._reports.clear(address);
    });
    s.opened.dispose();
  }

  coverage(address: string): Coverage | undefined {
    return this._sources.get(address)?.opened.coverage;
  }

  /** The sources whose opening failed, by address (the opener's word): cleared when one attaches or is let go. */
  private readonly _unopened = new Set<string>();

  /**
   * The opener says a source's opening failed (its issue its own, said where it opens): until a source attaches under
   * the address, its tiles read as failed, not pending (`unopened`), so a layer waits for nothing (_spec/README.md
   * "Data", Store, A failed opening).
   */
  markUnopened(address: string): void {
    if (!this.isRetired && !this._sources.attached(address)) this._unopened.add(address);
  }

  /** Whether a source's opening failed and none is attached under its address. */
  unopened(address: string): boolean {
    return this._unopened.has(address);
  }

  // Entries

  entry(key: string): ReadonlyState<DataEntry | undefined> {
    return this._entries.key(key);
  }

  get(key: string): DataEntry | undefined {
    return this._entries.get(key);
  }

  node(source: string, key: string): ReadonlyState<NodeEntry | undefined> {
    return this._nodes.entries.key(nodeKey(source, key));
  }

  getNode(source: string, key: string): NodeEntry | undefined {
    return this._nodes.entries.get(nodeKey(source, key));
  }

  /** Every node entry held now, by `<source>#<key>` (live: read it, never keep it). */
  get nodes(): ReadonlyMap<string, NodeEntry> {
    return this._nodes.entries.values;
  }

  /** Holds a node (a placement shown, the selection, a highlight) until the returned release is called. */
  hold(source: string, key: string): () => void {
    const k = nodeKey(source, key);
    this._nodes.hold(k);
    let held = true;
    return () => {
      if (!held) return;
      held = false;
      this._run(() => this._nodes.release(k));
    };
  }

  /**
   * The tiles a want covers: the keys of the entries it reads, in position order; none for an unknown source or system,
   * and none past `maxTilesPerWant` tiles (such a want is not asked for either: `seq.store.too_many_tiles`), so a
   * reader of the whole of a long system at a fine level does not build a key a tile.
   */
  keysOf(want: Want): string[] {
    // a 1D want's keys straight from its tiles' range: a layout lists them every frame, no tile object made
    if (!want.y) {
      const t = this._tiles(want);
      if (!t || t.count > this.settings.maxTilesPerWant) return [];
      const parts = want.side !== undefined || want.items === true ? { side: want.side, items: want.items } : undefined;
      const out: string[] = [];
      for (let i = t.first; i <= t.last; i++) out.push(dataKey(want.source, t.system, t.level, i, parts));
      return out;
    }
    const tiles = this.tilesOf(want);
    // a symmetric grid's mirrored tiles share a key: each once
    return [...new Set(tiles.map((t) => t.key))];
  }

  /**
   * The tiles a want covers, x-major (a 2D want's x tiles × y tiles): each its key, its region (a symmetric grid's served
   * tile) and whether the asked tile is that one's mirror; none as `keysOf`.
   */
  tilesOf(want: Want): WantTile[] {
    const t = this._tiles(want);
    if (!t || t.count > this.settings.maxTilesPerWant) return [];
    return eachTile(want, t);
  }

  /**
   * What a layer's wants cover together and whether the store asks them (each within `maxTilesPerWant`, together
   * within `maxTilesPerLayer` tiles and `maxBytesPerLayer` of read cost): a layer reads coarser until they fit.
   */
  layerLoad(wants: readonly Want[]): LayerLoad & { fits: boolean } {
    return layerLoad(wants, this._sources, this.settings);
  }

  /**
   * Whether a want covers more than `maxTilesPerWant` tiles (not asked, `seq.store.too_many_tiles`; `keysOf` lists none
   * of them): what tells a reader's empty list past the bound (pending) from one past a source's end (nothing there).
   */
  tooMany(want: Want): boolean {
    const t = this._tiles(want);
    return t !== undefined && t.count > this.settings.maxTilesPerWant && !knownEmpty(want, t, this._sources);
  }

  // Requests

  /** A new consumer's want set (an export's, beside the live view's), merged with the others by best priority. */
  wanter(): Wanter {
    const set = new WantSet(this);
    this._wanters.add(set);
    return set;
  }

  /**
   * Whether a layer's wants together (`Want.layer`) were last refused, past `maxTilesPerLayer` (`tiles`) or
   * `maxBytesPerLayer` (`cost`): what draws its "zoom in" state.
   */
  layerRefused(layer: string): 'tiles' | 'cost' | undefined {
    return this._refused.get(layer);
  }

  /** Moves when a layer's refusal comes or goes: what a reader of `layerRefused` keeps to know it read it again. */
  get refusals(): number {
    return this._refusals;
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
      for (const address of this._sources.addresses()) this.detach(address);
    });
  }

  /** @internal Plans a consumer's wants into tiles. */
  _plan(wants: readonly Want[]): Map<string, PlannedTile> {
    const refused = new Map<string, 'tiles' | 'cost'>();
    const planned = planWants(wants, this._sources, this.settings, ({ want: w, bound, value, max }) => {
      if (w.layer !== undefined && bound !== 'want') refused.set(w.layer, bound === 'layer' ? 'tiles' : 'cost');
      const message =
        bound === 'want'
          ? `A wanted region covers ${value} tiles, more than ${max}`
          : bound === 'layer'
            ? `A layer's wanted regions cover ${value} tiles together, more than ${max}`
            : `A layer's wanted regions cost ${value} bytes to read together, more than ${max}`;
      const address = channelSource(w.source);
      this._run(() =>
        this._reports.report(address, issueOf(StoreIssues, 'seq.store.too_many_tiles', message, { address })),
      );
    });
    // a layer's refusal moves the version when it changes, so what draws the layer reads it again
    let changed = false;
    for (const w of wants) {
      if (w.layer === undefined) continue;
      const now = refused.get(w.layer);
      if (this._refused.get(w.layer) !== now) changed = true;
      if (now) this._refused.set(w.layer, now);
      else this._refused.delete(w.layer);
    }
    if (changed) {
      this._refusals++;
      this._version.set(this._next());
    }
    return planned;
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
          // a tile of a source detached since it was planned: no entry pending that nothing can answer
          if (!this._sources.get(t.region.source)) continue;
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
      const onChange = () => {
        if (!check()) return;
        stop();
        if (timer !== undefined) this._clock.clearTimeout(timer);
        resolve(true);
      };
      // the settled state moves alone too (a source's summaries built: no entry changes)
      const stops = [this._version.listen(onChange), this._settled.listen(onChange)];
      const stop = () => {
        for (const s of stops) s();
      };
      if (timeoutMs !== undefined)
        timer = this._clock.setTimeout(() => {
          stop();
          resolve(false);
        }, timeoutMs);
    });
  }

  // Tiles and slots

  private _tiles(w: Want) {
    return wantTiles(w, this._sources, this.settings);
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
    for (const k of nodes) this._nodes.release(k);
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
    const source = this._sources.get(slot.region.source, true);
    // never left queued out of the queue: a slot whose source is gone goes with its pending entry
    if (!source) return this._drop(slot);
    slot.phase = 'flight';
    this._inFlight++;
    const controller = new AbortController();
    slot.controller = controller;
    const token = ++slot.token;
    const request: RegionRequest = tileRequest(slot.region, slot.shapes);
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
    const bad = checkAnswer(answer, slot.request!, this.settings.maxRows, this.settings);
    if (bad) this._malformed(slot, bad);
    else if (!answer.ok) this._failure(slot, answer.kind, answer.issue);
    else if (answerBytes(answer.tables, answer.nodes) > this.settings.maxTileBytes) {
      // one answer never more than `maxTileBytes`, a quarter of the cache (_spec/README.md "Data", Store, Bytes)
      const message = `An answer holds more than ${this.settings.maxTileBytes} bytes`;
      const hint = SourceIssues.codes['seq.source.too_large'].hint;
      this._failure(slot, 'too-large', issueOf(SourceIssues, 'seq.source.too_large', message, { hint }));
    } else this._loaded(slot, answer.tables, answer.nodes);
  }

  private _onError(slot: Slot, error: unknown): void {
    // An abort the store did not ask for is the source's failure.
    const message = isAbortError(error) ? 'The read was aborted by the source' : 'The read failed';
    this._failure(slot, 'failed', issueOf(SourceIssues, 'seq.source.failed', message));
  }

  private _loaded(slot: Slot, tables: Tables, nodes: readonly DynamicNode[] | undefined): void {
    slot.attempts = 0;
    // The new nodes first, then the old released: a node in both is never dropped between.
    const held = this._nodes.intake(slot.region.source, nodes);
    this._free(slot);
    slot.nodes = held;
    for (const k of this._sources.goneWith(slot.key, new Set(held))) this._nodes.remove(k);
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
    const source = channelSource(slot.region.source);
    this._reports.report(source, { ...issue, address: issue.address ?? source });
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

  private _updateSettled(): void {
    let settled = true;
    for (const slot of this._slots.values()) {
      if (slot.rank !== Unwanted && this._entries.get(slot.key)?.status === 'pending') {
        settled = false;
        break;
      }
    }
    // a source's summaries building (a read let go may have started them) keep it unsettled: what they make is drawn
    // after, and an export waits for it (EC-125); their end runs the checks again (and the bound, their bytes in)
    const building = settled ? this._sources.building : undefined;
    if (building) {
      settled = false;
      if (!this._awaited.has(building)) {
        this._awaited.add(building);
        const again = () => this._run(() => {});
        building.then(again, again);
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
    // a pyramid is not evicted (it goes with its source), but it leaves the entries less room
    const max = this.settings.maxBytes - this.summaryBytes;
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
      this._reports.report(StoreWide, issueOf(StoreIssues, 'seq.store.over_budget', message));
    } else if (!over && this._overBudget.value) this._reports.clear(StoreWide);
    this._overBudget.set(over);
  }
}
