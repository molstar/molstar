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
import { generateStateTransition } from '@molstar/mvs/helpers/animation';
import { loadMVS } from '@molstar/mvs/load';
import { MVSRuntimeRegistry } from '@molstar/mvs/registry';

const Fixtures = ['features', 'motm1', 'exosome'];
const Update = process.env.MVS_UPDATE_GOLDENS === '1';
const FixtureDir = path.join(import.meta.dirname, 'fixtures');

/** [durationMs, hash] per frame */
type FrameGoldens = [number, string][];
interface Goldens {
  /** per snapshot: the MVS frames, `[]` without an animation */
  expanded: FrameGoldens[];
  /** per snapshot: the state tree of the snapshot itself and of its transition's frames */
  loaded: { base: string; frames: FrameGoldens }[];
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

async function createPlugin() {
  const plugin = new PluginContext({ registry: [MVSRuntimeRegistry], behaviors: [PluginSpec.Behavior(MolViewSpec)] });
  await plugin.init();
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
    // (`keepCamera`: the plugin has no canvas)
    await loadMVS(plugin, data, { keepCamera: true, doNotReportErrors: true });
    const loaded = plugin.managers.snapshot.state.entries.toArray().map((entry) => ({
      base: stateTreeHash(entry.snapshot.data),
      frames: (entry.snapshot.transition?.frames ?? []).map(
        (f) => [f.durationInMs, stateTreeHash(f.data)] as [number, string],
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
