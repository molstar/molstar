/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Goldens of loading the animation fixtures (`_test/fixtures/*.mvsj`): the MVS frames `generateStateTransition` makes
 * and the state snapshots `loadMVS` makes from them, with their durations (`testing.md §2`). A frame is a fingerprint:
 * a hash of its canonical JSON with the numbers left out (keys, strings, refs, which transforms there are), and the
 * numbers as a count and sums compared with a relative tolerance, since math functions (`sin`, `pow`) can differ in
 * the last bit between platforms. Versions, hashes of exact numbers, are replaced by their equality pattern (the n-th
 * distinct version of the snapshot). `Math.random` (noise) is seeded. Regenerate after an intended change of the output:
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

/** [hash of the JSON without numbers, how many numbers, their sum, their position-weighted sum, the sum of |x|] */
type Fingerprint = [string, number, number, number, number];
/** [durationMs, the MVS tree] per frame */
type FrameGoldens = [number, Fingerprint][];
/** [durationMs, the state tree, the camera and canvas (with `include_camera` / `include_canvas`)] per frame */
type LoadedFrameGoldens = [number, Fingerprint, Fingerprint][];
interface Goldens {
  /** per snapshot: the MVS frames, `[]` without an animation */
  expanded: FrameGoldens[];
  /** per snapshot: the state tree of the snapshot itself and of its transition's frames */
  loaded: { base: Fingerprint; frames: LoadedFrameGoldens }[];
}

/** Canonical JSON (sorted keys) with each finite number written as `#` and pushed to `numbers`. */
function canonicalJson(value: unknown, numbers: number[]): string {
  if (typeof value === 'number' && Number.isFinite(value)) {
    numbers.push(value);
    return '#';
  }
  if (value === null || typeof value !== 'object') return JSON.stringify(value) ?? 'undefined';
  if (Array.isArray(value)) return `[${value.map((v) => canonicalJson(v, numbers)).join(',')}]`;
  const keys = Object.keys(value)
    .filter((k) => (value as any)[k] !== undefined)
    .sort();
  return `{${keys.map((k) => `${JSON.stringify(k)}:${canonicalJson((value as any)[k], numbers)}`).join(',')}}`;
}

function fingerprint(value: unknown): Fingerprint {
  const numbers: number[] = [];
  const text = canonicalJson(value, numbers);
  let sum = 0;
  let weighted = 0;
  let abs = 0;
  for (let i = 0; i < numbers.length; i++) {
    sum += numbers[i];
    weighted += numbers[i] * ((i % 97) + 1);
    abs += Math.abs(numbers[i]);
  }
  return [createHash('sha1').update(text).digest('hex').slice(0, 16), numbers.length, sum, weighted, abs];
}

/** The same JSON but for the last bits of numbers. */
function sameFingerprint(a: Fingerprint, b: Fingerprint) {
  const tolerance = 1e-9 * (Math.max(a[4], b[4]) + 1);
  return (
    a[0] === b[0] &&
    a[1] === b[1] &&
    Math.abs(a[2] - b[2]) <= tolerance &&
    Math.abs(a[3] - b[3]) <= 97 * tolerance &&
    Math.abs(a[4] - b[4]) <= tolerance
  );
}

/**
 * The state tree with each version as its index among the distinct versions met so far in `versions` (shared by a
 * snapshot and its frames), and without the root's version (the plugin's, a fresh id per plugin).
 */
function stateTreeFingerprint(data: { tree: { transforms: any[] } } | undefined, versions: Map<string, number>) {
  if (!data) return fingerprint('none');
  return fingerprint(
    data.tree.transforms.map((t) => {
      if (t.ref === t.parent) return { ...t, version: undefined };
      let v = versions.get(t.version);
      if (v === undefined) versions.set(t.version, (v = versions.size));
      return { ...t, version: `v${v}` };
    }),
  );
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
          result.push((transition?.frames ?? []).map(([tree, duration]) => [duration, fingerprint(tree)]));
        }
        return result;
      }),
    );

    vi.spyOn(Math, 'random').mockImplementation(mulberry32(1));
    const fallbacks = SubtreeMemo.stats.fallbacks;
    await loadMVS(plugin, data, { doNotReportErrors: true });
    // (a replayed subtree never had to be loaded again: the rules of `load.md §4` hold for the fixtures)
    expect(SubtreeMemo.stats.fallbacks).toBe(fallbacks);
    const loaded = plugin.managers.snapshot.state.entries.toArray().map((entry) => {
      const versions = new Map<string, number>();
      return {
        base: stateTreeFingerprint(entry.snapshot.data, versions),
        frames: (entry.snapshot.transition?.frames ?? []).map(
          (f) =>
            [
              f.durationInMs,
              stateTreeFingerprint(f.data, versions),
              fingerprint({ camera: f.camera, canvas3d: f.canvas3d }),
            ] as LoadedFrameGoldens[number],
        ),
      };
    });
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
      const near = (a: number, b: number) => Math.abs(a - b) <= 1e-9 * (Math.abs(b) + 1);
      expect(actual.expanded.map((s) => s.length)).toEqual(expected.expanded.map((s) => s.length));
      for (let i = 0; i < expected.expanded.length; i++) {
        const frame = expected.expanded[i].findIndex(
          ([d, f], j) => !near(actual.expanded[i][j][0], d) || !sameFingerprint(actual.expanded[i][j][1], f),
        );
        expect(frame, `${name}: snapshot ${i}, first differing MVS frame`).toBe(-1);
      }
      expect(actual.loaded.length).toBe(expected.loaded.length);
      for (let i = 0; i < expected.loaded.length; i++) {
        expect(sameFingerprint(actual.loaded[i].base, expected.loaded[i].base), `${name}: state of snapshot ${i}`).toBe(
          true,
        );
        expect(actual.loaded[i].frames.length, `${name}: frames of snapshot ${i}`).toBe(
          expected.loaded[i].frames.length,
        );
        const frame = expected.loaded[i].frames.findIndex(
          ([d, state, view], j) =>
            !near(actual.loaded[i].frames[j][0], d) ||
            !sameFingerprint(actual.loaded[i].frames[j][1], state) ||
            !sameFingerprint(actual.loaded[i].frames[j][2], view),
        );
        expect(frame, `${name}: snapshot ${i}, first differing state frame`).toBe(-1);
      }
    });
  }
});
