/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Goldens of loading the animation fixtures (`_test/fixtures/*.mvsj`): the MVS frames `generateStateTransition` makes
 * and the state snapshots `loadMVS` makes from them, as a hash of each frame's canonical JSON and its duration. They
 * pin the output while the load is optimized (`.v6/mvs-optimizations`): any change of a value, a key, a ref, a
 * version, or of which time steps merge into one frame shows as a different hash at that frame. `Math.random` (noise)
 * is seeded. Regenerate after an intended change of the output:
 *
 *   MVS_UPDATE_GOLDENS=1 pnpm test:unit packages/mvs/runtime/src/_test/load-animation.golden.test.ts
 */

import { createHash } from 'node:crypto';
import * as fs from 'node:fs';
import * as path from 'node:path';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { PluginCommands } from '@molstar/plugin/commands';
import { PluginContext } from '@molstar/plugin/context';
import { PluginSpec } from '@molstar/plugin/spec';
import { Task } from '@molstar/core/task';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import { MolViewSpec } from '@molstar/mvs/behavior';
import { Canvas3DParams } from '@molstar/graphics/canvas3d/canvas3d';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { generateStateTransition } from '@molstar/mvs/helpers/animation';
import { loadMVS } from '@molstar/mvs/load';
import { SubtreeMemo } from '@molstar/mvs/load-generic';
import { MVSRuntimeRegistry } from '@molstar/mvs/registry';

const Fixtures = ['features', 'motm1', 'exosome'];
const Update = process.env.MVS_UPDATE_GOLDENS === '1';
const FixtureDir = path.join(import.meta.dirname, 'fixtures');

/** [durationMs, hash] per frame */
type FrameGoldens = [number, string][];
/** [durationMs, hash of the state tree, hash of the camera and canvas (with `include_camera` / `include_canvas`)] */
type LoadedFrameGoldens = [number, string, string][];
interface Goldens {
  /** per snapshot: the MVS frames, `[]` without an animation */
  expanded: FrameGoldens[];
  /** per snapshot: the state tree of the snapshot itself and of its transition's frames */
  loaded: { base: string; frames: LoadedFrameGoldens }[];
}

function canonicalJson(value: unknown): string {
  if (value === null || typeof value !== 'object') return JSON.stringify(value) ?? 'undefined';
  if (Array.isArray(value)) return `[${value.map(canonicalJson).join(',')}]`;
  const keys = Object.keys(value)
    .filter((k) => (value as any)[k] !== undefined)
    .sort();
  return `{${keys.map((k) => `${JSON.stringify(k)}:${canonicalJson((value as any)[k])}`).join(',')}}`;
}

function hash(value: unknown) {
  return createHash('sha1').update(canonicalJson(value)).digest('hex').slice(0, 16);
}

/** The state tree without the root's `version`, which is the plugin's (a fresh id per plugin). */
function stateTreeHash(data: { tree: { transforms: any[] } } | undefined) {
  if (!data) return 'none';
  return hash(data.tree.transforms.map((t) => (t.ref === t.parent ? { ...t, version: undefined } : t)));
}

function mulberry32(seed: number) {
  return () => {
    seed |= 0;
    seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/** Enough of a canvas for the camera and canvas props of snapshots and frames (`createPluginStateSnapshotCamera`). */
function stubCanvas3D() {
  const snapshot = {
    mode: 'perspective',
    fov: Math.PI / 4,
    position: [0, 0, 100],
    up: [0, 1, 0],
    target: [0, 0, 0],
    radius: 20,
    radiusMax: 20,
    fog: 50,
    clipFar: true,
    minNear: 5,
    minFar: 0,
  };
  return {
    props: PD.getDefaultValues(Canvas3DParams),
    setProps() {},
    dispose() {},
    camera: { state: snapshot, getSnapshot: () => ({ ...snapshot }) },
  };
}

async function createPlugin() {
  const plugin = new PluginContext({ registry: [MVSRuntimeRegistry], behaviors: [PluginSpec.Behavior(MolViewSpec)] });
  await plugin.init();
  (plugin as any).canvas3d = stubCanvas3D();
  // (applying a snapshot would download the structures)
  const dispatch = plugin.commands.dispatch.bind(plugin.commands);
  plugin.commands.dispatch = ((cmd: any, params: any) =>
    cmd === PluginCommands.State.Snapshots.Apply ? Promise.resolve() : dispatch(cmd, params)) as any;
  return plugin;
}

async function computeGoldens(name: string): Promise<Goldens> {
  const data = MVSData.fromMVSJ(fs.readFileSync(path.join(FixtureDir, `${name}.mvsj`), 'utf8'));
  const snapshots = data.kind === 'multiple' ? data.snapshots : MVSData.stateToStates(data).snapshots;
  const plugin = await createPlugin();
  try {
    vi.spyOn(Math, 'random').mockImplementation(mulberry32(1));
    const expanded = await plugin.runTask(
      Task.create('expand', async (ctx) => {
        const result: FrameGoldens[] = [];
        for (let i = 0; i < snapshots.length; i++) {
          const transition = await generateStateTransition(ctx, snapshots[i], i, snapshots.length);
          result.push((transition?.frames ?? []).map(([tree, duration]) => [duration, hash(tree)]));
        }
        return result;
      }),
    );

    vi.spyOn(Math, 'random').mockImplementation(mulberry32(1));
    const fallbacks = SubtreeMemo.stats.fallbacks;
    await loadMVS(plugin, data, { doNotReportErrors: true });
    // (a replayed subtree never had to be loaded again: the rules of `load.md §4` hold for the fixtures)
    expect(SubtreeMemo.stats.fallbacks).toBe(fallbacks);
    const loaded = plugin.managers.snapshot.state.entries.toArray().map((entry) => ({
      base: stateTreeHash(entry.snapshot.data),
      frames: (entry.snapshot.transition?.frames ?? []).map(
        (f) =>
          [
            f.durationInMs,
            stateTreeHash(f.data),
            hash({ camera: f.camera, canvas3d: f.canvas3d }),
          ] as LoadedFrameGoldens[number],
      ),
    }));
    return { expanded, loaded };
  } finally {
    plugin.dispose();
  }
}

describe('MVS animation load goldens', () => {
  afterEach(() => {
    vi.restoreAllMocks();
  });

  for (const name of Fixtures) {
    it(name, { timeout: 120_000 }, async () => {
      const goldenPath = path.join(FixtureDir, 'goldens', `${name}.json`);
      const actual = await computeGoldens(name);
      if (Update) {
        fs.mkdirSync(path.dirname(goldenPath), { recursive: true });
        fs.writeFileSync(goldenPath, JSON.stringify(actual));
        return;
      }
      const expected: Goldens = JSON.parse(fs.readFileSync(goldenPath, 'utf8'));
      expect(actual.expanded.map((s) => s.length)).toEqual(expected.expanded.map((s) => s.length));
      for (let i = 0; i < expected.expanded.length; i++) {
        expect(actual.expanded[i], `${name}: expanded frames of snapshot ${i}`).toEqual(expected.expanded[i]);
      }
      expect(actual.loaded.length).toBe(expected.loaded.length);
      for (let i = 0; i < expected.loaded.length; i++) {
        expect(actual.loaded[i].base, `${name}: state of snapshot ${i}`).toBe(expected.loaded[i].base);
        expect(actual.loaded[i].frames, `${name}: transition frames of snapshot ${i}`).toEqual(
          expected.loaded[i].frames,
        );
      }
    });
  }
});
