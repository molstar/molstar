/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { Task } from '@molstar/core/task';
import type { ParticleList, ParticleTarget } from '../particle-list.js';
import { createParticleTileLoader } from '../particle-tile-loader.js';
import { type ParticleTileProvider, type ParticleTileSet, ParticleTiles } from '../particle-tiles.js';

const Count = 64;

/** 64 particles on a 4³ lattice over 0-30 Å, targets 1 and 2 alternating by x: 8 tiles of 2³, each using both */
function createList(): ParticleList {
  const coordinates = new Float32Array(Count * 3);
  const targets = new Int32Array(Count);
  for (let i = 0; i < Count; ++i) {
    const x = i % 4;
    const y = Math.floor(i / 4) % 4;
    const z = Math.floor(i / 16);
    coordinates.set([x * 10, y * 10, z * 10], i * 3);
    targets[i] = x % 2 ? 2 : 1;
  }
  const sites = (n: number): ParticleTarget => ({
    kind: 'sites',
    sites: { count: n, positions: new Float32Array(n * 3), radii: new Float32Array(n).fill(1) },
  });
  return {
    label: 'lattice',
    count: Count,
    targets,
    targetInfo: new Map([
      [1, { name: 'one' }],
      [2, { name: 'two' }],
    ]),
    targetMapping: new Map([
      [1, sites(2)],
      [2, sites(3)],
    ]),
    coordinates,
  };
}

const Delay = 10;

async function flush() {
  for (let i = 0; i < 8; ++i) await new Promise((r) => setImmediate(r));
}

/** the tasks run up to their wait, then `ms` passed, then they run on */
async function step(ms = Delay) {
  await flush();
  await vi.advanceTimersByTimeAsync(ms);
  await flush();
}

/** the set's provider behind delayed tasks that count their calls and aborts (a task updates before and after its wait) */
function delayed(set: ParticleTileSet) {
  const calls = { tiles: [] as string[], targets: [] as number[], aborted: [] as string[] };
  const fail = new Set<string>();
  const task = <T>(name: string, value: () => T) =>
    Task.create(name, async (ctx) => {
      await ctx.update(name);
      await new Promise((r) => setTimeout(r, Delay));
      try {
        await ctx.update(name);
      } catch (e) {
        calls.aborted.push(name);
        throw e;
      }
      if (fail.has(name)) throw new Error(`failed ${name}`);
      return value();
    });
  const provider: ParticleTileProvider = {
    getTile(id, frame) {
      calls.tiles.push(`${id}@${frame}`);
      return task(`${id}@${frame}`, () => set.provider.getTile(id, frame) as ParticleList);
    },
    getTarget(id) {
      calls.targets.push(id);
      return task(`t${id}`, () => set.provider.getTarget(id) as ParticleTarget);
    },
  };
  return { set: { ...set, provider } as ParticleTileSet, calls, fail };
}

describe('particle tile loader', () => {
  // only timeouts faked: the tasks' updates go through the scheduler's real setImmediate, flushed around each step
  beforeEach(() => vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout'] }));
  afterEach(() => vi.useRealTimers());

  const base = ParticleTiles.fromList(createList(), [2, 2, 2]);
  const ids = [...base.tiles.keys()].sort((a, b) => a - b);

  it('loads by priority, at most the concurrency at once, a target once for all tiles', async () => {
    const { set, calls } = delayed(base);
    const loader = createParticleTileLoader(set, { concurrency: 2 });
    const order = [ids[5], ids[2], ids[7], ids[0]];
    loader.setWanted(order.map((id) => ({ id, frame: 0 })));
    expect(calls.tiles).toEqual([`${ids[5]}@0`, `${ids[2]}@0`]);
    expect(loader.getState(ids[7], 0)).toBe('none');
    expect(loader.getState(ids[5], 0)).toBe('loading');

    await step();
    expect(loader.getState(ids[5], 0)).toBe('ready');
    expect(loader.getState(ids[2], 0)).toBe('ready');
    expect(calls.tiles).toEqual(order.slice(0, 4).map((id) => `${id}@0`));

    const idle = loader.whenIdle();
    await step();
    await idle;
    expect(order.map((id) => loader.getState(id, 0))).toEqual(['ready', 'ready', 'ready', 'ready']);
    // both targets, each loaded once, one object for all tiles
    expect(calls.targets.sort()).toEqual([1, 2]);
    expect(loader.getTarget(1)).toBe(base.provider.getTarget(1));
    expect(loader.targets.size).toBe(2);
    const tile = loader.getTile(ids[5], 0)!;
    expect(tile.count).toBe(base.tiles.get(ids[5])!.count);
    expect(loader.getLoaded()).toMatchObject({ tiles: 4, particles: 32, targets: 2 });
    expect(loader.getReady(0).map((r) => r.id)).toEqual([ids[0], ids[2], ids[5], ids[7]]);
  });

  it('cancels a load no longer wanted, aborting its task, and starts the next', async () => {
    const { set, calls } = delayed(base);
    const loader = createParticleTileLoader(set, { concurrency: 1 });
    let changes = 0;
    loader.subscribe(() => changes++);
    loader.setWanted([{ id: ids[0], frame: 0 }]);
    // its task's first update gives the loader its abort
    await flush();
    loader.setWanted([{ id: ids[1], frame: 0 }]);
    expect(loader.getState(ids[0], 0)).toBe('loading');
    expect(calls.tiles).toEqual([`${ids[0]}@0`]);

    await step();
    expect(calls.aborted).toEqual([`${ids[0]}@0`]);
    expect(loader.getState(ids[0], 0)).toBe('none');
    expect(loader.getError(ids[0], 0)).toBeUndefined();
    expect(calls.tiles).toEqual([`${ids[0]}@0`, `${ids[1]}@0`]);

    await step();
    expect(loader.getState(ids[1], 0)).toBe('ready');
    expect(changes).toBeGreaterThan(2);
  });

  it('tells its listeners of wants that start a load, not of wants alone (ready tiles in another order)', async () => {
    const { set } = delayed(base);
    const loader = createParticleTileLoader(set, { concurrency: 4 });
    const wants = (order: number[]) => loader.setWanted(order.map((id) => ({ id, frame: 0 })));
    wants([ids[0], ids[1]]);
    const idle = loader.whenIdle();
    await step();
    await step();
    await idle;
    let changes = 0;
    loader.subscribe(() => changes++);
    // the same tiles ready, in another order, or fewer of them: no state changed
    wants([ids[1], ids[0]]);
    wants([ids[1]]);
    expect(changes).toBe(0);
    // what is idle resolves at once still
    await loader.whenIdle();
    // a tile not loaded yet: its load starts, a change
    wants([ids[1], ids[2]]);
    expect(changes).toBe(1);
    expect(loader.getState(ids[2], 0)).toBe('loading');
  });

  it('keeps what is ready until released, then lets go of the unused targets', async () => {
    const { set, calls } = delayed(base);
    const loader = createParticleTileLoader(set);
    loader.setWanted([
      { id: ids[0], frame: 0 },
      { id: ids[1], frame: 0 },
    ]);
    await step();
    loader.setWanted([]);
    expect(loader.getState(ids[0], 0)).toBe('ready');
    expect(loader.getUnusedTargets()).toEqual([]);

    loader.release([{ id: ids[0], frame: 0 }]);
    expect(loader.getState(ids[0], 0)).toBe('none');
    expect(loader.getTile(ids[0], 0)).toBeUndefined();
    // the other tile still uses both
    loader.releaseTargets([1, 2]);
    expect(loader.targets.size).toBe(2);

    loader.release([{ id: ids[1], frame: 0 }]);
    expect(loader.getUnusedTargets()).toEqual([1, 2]);
    loader.releaseTargets([2]);
    expect([...loader.targets.keys()]).toEqual([1]);
    expect(loader.getLoaded()).toMatchObject({ tiles: 0, particles: 0, targets: 1 });

    // wanted again: the tile loaded again, the released target too, the kept one not
    loader.setWanted([{ id: ids[0], frame: 0 }]);
    await step();
    expect(loader.getState(ids[0], 0)).toBe('ready');
    expect(calls.targets.sort()).toEqual([1, 2, 2]);
  });

  it('marks a failed load, and loads it again where wanted after a release', async () => {
    const { set, fail } = delayed(base);
    const loader = createParticleTileLoader(set);
    fail.add(`${ids[3]}@0`);
    loader.setWanted([{ id: ids[3], frame: 0 }]);
    await step();
    await loader.whenIdle();
    expect(loader.getState(ids[3], 0)).toBe('failed');
    expect(String(loader.getError(ids[3], 0))).toContain('failed');
    expect(loader.getTile(ids[3], 0)).toBeUndefined();

    fail.clear();
    loader.release([{ id: ids[3], frame: 0 }]);
    loader.setWanted([{ id: ids[3], frame: 0 }]);
    await step();
    expect(loader.getState(ids[3], 0)).toBe('ready');
  });

  it('loads a trajectory frame wanted ahead, the frames apart', async () => {
    const frames = [createList(), createList()];
    const trajectory = { frameCount: 2, representative: frames[0], getFrameAtIndex: (i: number) => frames[i] };
    const { set, calls } = delayed(ParticleTiles.fromTrajectory(trajectory, [2, 2, 2]));
    const loader = createParticleTileLoader(set, { concurrency: 1 });
    loader.setWanted([
      { id: ids[0], frame: 0 },
      { id: ids[0], frame: 1 },
      { id: ids[0], frame: 2 }, // past the last frame: ignored
    ]);
    await step();
    expect(loader.getState(ids[0], 0)).toBe('ready');
    expect(loader.getState(ids[0], 1)).toBe('loading');
    await step();
    expect(loader.getState(ids[0], 1)).toBe('ready');
    expect(calls.tiles).toEqual([`${ids[0]}@0`, `${ids[0]}@1`]);
    expect(loader.getReady(1)).toEqual([{ id: ids[0], frame: 1 }]);
    expect(loader.getTile(ids[0], 1)).not.toBe(loader.getTile(ids[0], 0));
    // targets shared by the frames
    expect(calls.targets.sort()).toEqual([1, 2]);
    loader.dispose();
    expect(loader.getState(ids[0], 0)).toBe('none');
  });

  it("keeps a coarse level's node apart from the tiles of other frames (its id past the grid's cells)", async () => {
    const frames = [createList(), createList()];
    const trajectory = { frameCount: 2, representative: frames[0], getFrameAtIndex: (i: number) => frames[i] };
    const tiled = ParticleTiles.fromTrajectory(trajectory, [2, 2, 2]);
    // a node as a set's levels add one: the id after the 2³ grid's cells, a tile of the set of its own
    const node = 8;
    const tiles = new Map(tiled.tiles);
    tiles.set(node, tiled.tiles.get(ids[0])!);
    const getTile = (id: number, frame: number) => tiled.provider.getTile(id === node ? ids[0] : id, frame);
    const { set } = delayed({ ...tiled, tiles, provider: { ...tiled.provider, getTile } });
    const loader = createParticleTileLoader(set, { concurrency: 2 });
    // node 8 at frame 0 and tile 0 at frame 1: one key were a frame's keys the cells' alone
    loader.setWanted([
      { id: node, frame: 0 },
      { id: 0, frame: 1 },
    ]);
    await step();
    expect(loader.getState(node, 0)).toBe('ready');
    expect(loader.getState(0, 1)).toBe('ready');
    expect(loader.getState(0, 0)).toBe('none');
    expect(loader.getReady()).toEqual([
      { id: node, frame: 0 },
      { id: 0, frame: 1 },
    ]);
    loader.dispose();
  });
});
