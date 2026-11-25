/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from '@molstar/core/task';
import type { ParticleList, ParticleTarget } from './particle-list.js';
import { ParticleTileGrid, type ParticleTileSet } from './particle-tiles.js';

// The loader of a tile set's tiles and targets: what is wanted is the caller's (by
// priority), what to let go the app's (it measures memory): no budget here, the counts and bytes of what is loaded.

/** A tile's or a target's state */
export type ParticleTileLoadState = 'none' | 'loading' | 'ready' | 'failed';

/** A tile at a frame */
export interface ParticleTileRequest {
  readonly id: number;
  readonly frame: number;
}

export interface ParticleTileLoader {
  readonly set: ParticleTileSet;
  /**
   * The tiles wanted, by priority (the first first; a trajectory's next frames after the current one's: wanted ahead):
   * their loads started in that order, at most `concurrency` at once, a tile's targets with it (a target once for all
   * the tiles wanting it); a load of a tile no longer wanted cancelled (its task aborted where it can be, its result
   * dropped); what is ready kept, wanted or not (`release`)
   */
  setWanted(requests: readonly ParticleTileRequest[]): void;
  /** Lets go of tiles: ready, loading (cancelled) or failed, back to none */
  release(requests: readonly ParticleTileRequest[]): void;
  /** Lets go of targets no ready or loading tile uses (`getUnusedTargets`); the others kept */
  releaseTargets(ids: Iterable<number>): void;
  getState(id: number, frame: number): ParticleTileLoadState;
  /** A tile's list where it is ready: its particles, its targets loaded (`getTarget`) */
  getTile(id: number, frame: number): ParticleList | undefined;
  getError(id: number, frame: number): unknown;
  getTargetState(id: number): ParticleTileLoadState;
  getTarget(id: number): ParticleTarget | undefined;
  /** The loaded targets, one object each for all tiles */
  readonly targets: ReadonlyMap<number, ParticleTarget>;
  /** advanced as `targets` changes (a target loaded or let go): a view of the same version its mapping kept */
  readonly targetsVersion: number;
  /** The ready tiles, of a frame or of all */
  getReady(frame?: number): ParticleTileRequest[];
  /** The loaded targets no ready or loading tile uses */
  getUnusedTargets(): number[];
  /** What is loaded: the ready tiles, their particles, the loaded targets, and the bytes of both (estimates) */
  getLoaded(): { readonly tiles: number; readonly particles: number; readonly targets: number; readonly bytes: number };
  /** Called where a state changed */
  subscribe(listener: () => void): () => void;
  /** Resolves where no load runs and every wanted tile is ready or failed */
  whenIdle(): Promise<void>;
  dispose(): void;
}

export interface ParticleTileLoaderOptions {
  /** the tiles loading at once, 4 by default */
  readonly concurrency?: number;
  /** a tile's load failed (a request dropped, a short answer): told each time, the tile failed until wanted again */
  readonly onError?: (id: number, frame: number, error: unknown) => void;
}

/** a failed tile's loads at most, while wanted (`waitToRetry`, `setWanted`): a dropped request not a hole for good */
export const ParticleTileLoadAttempts = 3;
/**
 * the wait before a failed tile loads again, by its failures so far (1: before the second load), in ms: attempts a
 * frame apart would all fall within a short drop of the network
 */
export const ParticleTileRetryDelaysMs: readonly number[] = [1000, 3000];

interface TileEntry {
  readonly id: number;
  readonly frame: number;
  state: ParticleTileLoadState;
  list?: ParticleList;
  error?: unknown;
  /** set where a load runs: no longer wanted, its result dropped */
  cancelled?: boolean;
  abort?: () => void;
  /** its load's reads (`ParticleTileProvider.getTile`'s signal): aborted as it is cancelled */
  controller?: AbortController;
  /** its loads failed since it was last ready or released */
  failures: number;
  /** failed, waiting before it may load again (`ParticleTileRetryDelaysMs`): the timer that ends the wait */
  retryTimer?: ReturnType<typeof setTimeout>;
}

interface TargetEntry {
  state: ParticleTileLoadState;
  target?: ParticleTarget;
  error?: unknown;
  promise?: Promise<ParticleTarget>;
}

/** `value` resolved: as it is, or its task run (its abort kept where its progress offers one) */
function resolve<T>(value: T | Task<T>, onAbort?: (abort: () => void) => void): Promise<T> {
  if (!Task.is(value)) return Promise.resolve(value);
  return value.run((p) => onAbort?.(() => p.requestAbort('no longer wanted')), 50);
}

export function createParticleTileLoader(
  set: ParticleTileSet,
  options: ParticleTileLoaderOptions = {},
): ParticleTileLoader {
  const concurrency = Math.max(1, options.concurrency ?? 4);
  // a frame's keys past every tile id: a coarse level's nodes past the grid's cells (`ParticleTileLevels`, `cells + k`)
  let ids = ParticleTileGrid.cellCount(set.grid);
  for (const id of set.tiles.keys()) if (id >= ids) ids = id + 1;
  const keyOf = (id: number, frame: number) => frame * ids + id;
  const tiles = new Map<number, TileEntry>();
  const targetEntries = new Map<number, TargetEntry>();
  const loadedTargets = new Map<number, ParticleTarget>();
  let targetsVersion = 0;
  let wanted: number[] = [];
  let wantedSet = new Set<number>();
  let running = 0;
  let disposed = false;
  const listeners = new Set<() => void>();
  let idleWaiters: (() => void)[] = [];

  const isIdle = () => {
    if (running > 0) return false;
    for (const key of wanted) {
      const state = tiles.get(key)?.state ?? 'none';
      if (state === 'none' || state === 'loading') return false;
    }
    return true;
  };
  const notify = () => {
    for (const l of listeners) l();
    settleIdle();
  };
  const settleIdle = () => {
    if (idleWaiters.length > 0 && isIdle()) {
      const waiters = idleWaiters;
      idleWaiters = [];
      for (const w of waiters) w();
    }
  };

  const getEntry = (id: number, frame: number) => {
    const key = keyOf(id, frame);
    let entry = tiles.get(key);
    if (!entry) {
      entry = { id, frame, state: 'none', failures: 0 };
      tiles.set(key, entry);
    }
    return entry;
  };

  function loadTarget(id: number): Promise<ParticleTarget> {
    let entry = targetEntries.get(id);
    if (entry?.state === 'ready') return Promise.resolve(entry.target!);
    if (entry?.promise) return entry.promise;
    entry = { state: 'loading' };
    targetEntries.set(id, entry);
    const e = entry;
    e.promise = (async () => {
      try {
        const target = await resolve(set.provider.getTarget(id));
        e.state = 'ready';
        e.target = target;
        loadedTargets.set(id, target);
        targetsVersion += 1;
        return target;
      } catch (error) {
        e.state = 'failed';
        e.error = error;
        throw error;
      } finally {
        e.promise = undefined;
      }
    })();
    return e.promise;
  }

  function start(entry: TileEntry) {
    entry.state = 'loading';
    entry.cancelled = false;
    entry.error = undefined;
    running += 1;
    const tile = set.tiles.get(entry.id);
    (async () => {
      try {
        if (!tile) throw new Error(`particle tiles: no tile ${entry.id}`);
        // its targets with it, each once for all the tiles wanting it; those loaded already not asked again (a tile:
        // hundreds of targets, a promise each otherwise)
        let pending: Promise<ParticleTarget>[] | undefined;
        for (let i = 0; i < tile.targets.length; ++i) {
          const t = tile.targets[i];
          if (targetEntries.get(t)?.state === 'ready') continue;
          (pending ??= []).push(loadTarget(t));
        }
        const targets = pending ? Promise.all(pending) : undefined;
        targets?.catch(() => {}); // awaited below, after the tile: not an unhandled rejection meanwhile
        entry.controller = new AbortController();
        const list = await resolve(set.provider.getTile(entry.id, entry.frame, entry.controller.signal), (abort) => {
          entry.abort = abort;
        });
        if (targets) await targets;
        if (entry.cancelled || disposed) {
          entry.state = 'none';
          return;
        }
        entry.state = 'ready';
        entry.list = list;
        entry.failures = 0;
      } catch (error) {
        if (entry.cancelled || disposed) {
          entry.state = 'none';
        } else {
          entry.state = 'failed';
          entry.error = error;
          entry.failures += 1;
          if (entry.failures < ParticleTileLoadAttempts) waitToRetry(entry);
          options.onError?.(entry.id, entry.frame, error);
        }
      } finally {
        entry.abort = undefined;
        entry.controller = undefined;
        running -= 1;
        if (!disposed) {
          pump();
          notify();
        }
      }
    })();
  }

  /**
   * a failed tile's wait before it may load again: as it ends, loaded again where still wanted, the listeners told. Why
   * not by the next `setWanted`: a view asks again only as its pick changes (a camera still, a view not following
   * never), the tile failed for good otherwise
   */
  function waitToRetry(entry: TileEntry) {
    const delays = ParticleTileRetryDelaysMs;
    const ms = delays[Math.min(entry.failures, delays.length) - 1] ?? 0;
    entry.retryTimer = setTimeout(() => {
      entry.retryTimer = undefined;
      if (disposed) return;
      // no longer wanted: failed until wanted again (`setWanted`)
      if (entry.state === 'failed' && wantedSet.has(keyOf(entry.id, entry.frame))) {
        entry.state = 'none';
        entry.error = undefined;
        pump();
      }
      notify();
    }, ms);
  }
  function clearRetry(entry: TileEntry) {
    if (entry.retryTimer === undefined) return;
    clearTimeout(entry.retryTimer);
    entry.retryTimer = undefined;
  }

  /** the next wanted tiles started, by priority, while fewer than `concurrency` load */
  function pump() {
    for (const key of wanted) {
      if (running >= concurrency) return;
      const entry = tiles.get(key)!;
      if (entry.state === 'none') start(entry);
    }
  }

  function cancel(entry: TileEntry) {
    if (entry.state !== 'loading') return;
    entry.cancelled = true;
    entry.abort?.();
    // its reads dropped too (a request not sent, one in flight aborted where nothing else wants its bytes)
    entry.controller?.abort();
  }

  return {
    set,
    setWanted(requests) {
      // disposed: none loads again (a view of its set let go after it, `ParticleTilesView`'s dispose)
      if (disposed) return;
      const keys: number[] = [];
      const next = new Set<number>();
      for (const { id, frame } of requests) {
        if (!set.tiles.has(id) || frame < 0 || frame >= set.frameCount) continue;
        const key = keyOf(id, frame);
        if (next.has(key)) continue;
        next.add(key);
        keys.push(key);
        // failed: loaded again where wanted once its wait is over (a wait over while not wanted: here), a few times (a
        // dropped request), then failed until released
        const entry = getEntry(id, frame);
        if (entry.state === 'failed' && entry.failures < ParticleTileLoadAttempts && entry.retryTimer === undefined) {
          entry.state = 'none';
          entry.error = undefined;
        }
      }
      for (const key of wantedSet) {
        if (!next.has(key)) cancel(tiles.get(key)!);
      }
      wanted = keys;
      wantedSet = next;
      // the listeners told where a load started (a state changed), not of the wants alone: a view's pass asks every
      // time the camera moves, and a listener's pass woken by it would ask again (one a frame while it moves). A load
      // cancelled tells them as it stops
      const before = running;
      pump();
      if (running !== before) notify();
      else settleIdle();
    },

    release(requests) {
      for (const { id, frame } of requests) {
        const entry = tiles.get(keyOf(id, frame));
        if (!entry) continue;
        cancel(entry);
        if (entry.state !== 'loading') entry.state = 'none';
        entry.list = undefined;
        entry.error = undefined;
        entry.failures = 0;
        clearRetry(entry);
      }

      notify();
    },
    releaseTargets(ids) {
      const unused = new Set(this.getUnusedTargets());
      for (const id of ids) {
        if (!unused.has(id)) continue;
        targetEntries.delete(id);
        if (loadedTargets.delete(id)) targetsVersion += 1;
      }
      notify();
    },
    getState: (id, frame) => tiles.get(keyOf(id, frame))?.state ?? 'none',
    getTile: (id, frame) => {
      const entry = tiles.get(keyOf(id, frame));
      return entry?.state === 'ready' ? entry.list : undefined;
    },
    getError: (id, frame) => tiles.get(keyOf(id, frame))?.error,
    getTargetState: (id) => targetEntries.get(id)?.state ?? 'none',
    getTarget: (id) => loadedTargets.get(id),
    targets: loadedTargets,
    get targetsVersion() {
      return targetsVersion;
    },
    getReady(frame) {
      const out: ParticleTileRequest[] = [];
      for (const entry of tiles.values()) {
        if (entry.state === 'ready' && (frame === undefined || entry.frame === frame)) {
          out.push({ id: entry.id, frame: entry.frame });
        }
      }
      return out.sort((a, b) => a.frame - b.frame || a.id - b.id);
    },
    getUnusedTargets() {
      const used = new Set<number>();
      for (const entry of tiles.values()) {
        if (entry.state !== 'ready' && entry.state !== 'loading') continue;
        for (const t of set.tiles.get(entry.id)?.targets ?? []) used.add(t);
      }
      return [...loadedTargets.keys()].filter((id) => !used.has(id)).sort((a, b) => a - b);
    },
    getLoaded() {
      let count = 0;
      let particles = 0;
      let bytes = 0;
      for (const entry of tiles.values()) {
        if (entry.state !== 'ready') continue;
        const tile = set.tiles.get(entry.id)!;
        count += 1;
        particles += tile.count;
        bytes += tile.byteCount ?? 0;
      }
      for (const id of loadedTargets.keys()) bytes += set.catalog.targets.get(id)?.byteCount ?? 0;
      return { tiles: count, particles, targets: loadedTargets.size, bytes };
    },
    subscribe(listener) {
      listeners.add(listener);
      return () => listeners.delete(listener);
    },
    whenIdle() {
      if (isIdle()) return Promise.resolve();
      return new Promise<void>((resolve) => idleWaiters.push(resolve));
    },
    dispose() {
      disposed = true;
      for (const entry of tiles.values()) {
        cancel(entry);
        clearRetry(entry);
      }
      tiles.clear();
      wanted = [];
      wantedSet = new Set();
      targetEntries.clear();
      loadedTargets.clear();
      listeners.clear();
      for (const w of idleWaiters) w();
      idleWaiters = [];
    },
  };
}
