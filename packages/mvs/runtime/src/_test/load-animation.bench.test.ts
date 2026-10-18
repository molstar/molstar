/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Benchmark for loading MVS documents with animations into plugin state snapshots (`loadMVS` without applying a
 * snapshot, so nothing is downloaded). Fixtures are the `index.mvsj` of real MVSX files (`_test/fixtures`). Skipped
 * unless `MVS_BENCH=1`:
 *
 *   MVS_BENCH=1 pnpm test:unit packages/mvs/runtime/src/_test/load-animation.bench.test.ts
 *
 * `MVS_BENCH_FIXTURES=motm1,exosome` picks fixtures; `MVS_BENCH_PROFILE=<dir>` also writes a `.cpuprofile` of one
 * `loadMVS` per fixture. Each phase prints best-of-N time and bytes allocated by one run (`v8`
 * `total_allocated_bytes`); `retained` is the heap held by the loaded snapshot entries after GC.
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { Session } from 'node:inspector/promises';
import { getHeapStatistics, setFlagsFromString } from 'node:v8';
import { runInNewContext } from 'node:vm';
import { describe, it } from 'vitest';
import { PluginContext } from '@molstar/plugin/context';
import { PluginCommands } from '@molstar/plugin/commands';
import { PluginSpec } from '@molstar/plugin/spec';
import { Task } from '@molstar/core/task';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import { MolViewSpec } from '@molstar/mvs/behavior';
import { generateStateTransition } from '@molstar/mvs/helpers/animation';
import { loadMVS } from '@molstar/mvs/load';
import { MVSRuntimeRegistry } from '@molstar/mvs/registry';
import { convertMvsToMolstar } from '@molstar/mvs/tree/molstar/conversion';

const Enabled = process.env.MVS_BENCH === '1';
const Fixtures = (process.env.MVS_BENCH_FIXTURES ?? 'features,motm1,exosome').split(',');
const ProfileDir = process.env.MVS_BENCH_PROFILE;
const NRuns = Number(process.env.MVS_BENCH_RUNS ?? 3);

setFlagsFromString('--expose-gc');
const gc: () => void = runInNewContext('gc');

function allocated() {
  return getHeapStatistics().total_allocated_bytes ?? 0;
}

async function measure<T>(f: () => T | Promise<T>): Promise<{ result: T; ms: number; mb: number }> {
  let best = Infinity;
  let mb = 0;
  let result!: T;
  for (let i = 0; i < NRuns; i++) {
    gc();
    const a0 = allocated();
    const t0 = performance.now();
    result = await f();
    const dt = performance.now() - t0;
    if (i === 0) mb = (allocated() - a0) / 1e6;
    best = Math.min(best, dt);
  }
  return { result, ms: best, mb };
}

function readFixture(name: string) {
  return MVSData.fromMVSJ(fs.readFileSync(path.join(import.meta.dirname, 'fixtures', `${name}.mvsj`), 'utf8'));
}

async function createPlugin() {
  const plugin = new PluginContext({ registry: [MVSRuntimeRegistry], behaviors: [PluginSpec.Behavior(MolViewSpec)] });
  await plugin.init();
  // (applying a snapshot would download the structures: the bench measures building the snapshots only)
  const dispatch = plugin.commands.dispatch.bind(plugin.commands);
  plugin.commands.dispatch = ((cmd: any, params: any) =>
    cmd === PluginCommands.State.Snapshots.Apply ? Promise.resolve() : dispatch(cmd, params)) as any;
  return plugin;
}

// (`keepCamera`: the plugin has no canvas; the camera of a frame is computed from its nodes only, without the scene)
const LoadOptions = { doNotReportErrors: true, keepCamera: true };

function fmt(label: string, r: { ms: number; mb: number }, extra = '') {
  return `${label.padEnd(34)} ${r.ms.toFixed(0).padStart(7)} ms ${r.mb.toFixed(0).padStart(7)} MB alloc ${extra}`;
}

describe.skipIf(!Enabled)('MVS animation load bench', () => {
  for (const name of Fixtures) {
    it(name, { timeout: 600_000 }, async () => {
      const data = readFixture(name);
      const snapshots = data.kind === 'multiple' ? data.snapshots : [{ root: data.root, metadata: data.metadata }];
      const plugin = await createPlugin();
      const lines: string[] = [`== ${name}: ${snapshots.length} snapshots`];

      // 1. expanding the animations into MVS trees, one per frame
      const expand = await measure(() =>
        plugin.runTask(
          Task.create('expand', async (ctx) => {
            const all = [];
            for (let i = 0; i < snapshots.length; i++) {
              all.push(await generateStateTransition(ctx, snapshots[i] as any, i, snapshots.length));
            }
            return all;
          }),
        ),
      );
      const frames = expand.result.flatMap((t) => t?.frames ?? []);
      lines.push(fmt('generateStateTransition', expand, `(${frames.length} distinct frames)`));

      // 2. MVS tree -> Molstar tree, per frame
      const convert = await measure(() => frames.map((f) => convertMvsToMolstar(f[0], undefined)));
      lines.push(fmt('convertMvsToMolstar (frames)', convert));

      // 3. the whole load (validation, expansion, conversion, Molstar tree -> state snapshot, camera)
      const load = await measure(() => loadMVS(plugin, data, LoadOptions));
      lines.push(fmt('loadMVS', load));

      // retained: the snapshot entries the load leaves in the manager
      plugin.managers.snapshot.clear();
      gc();
      const h0 = process.memoryUsage().heapUsed;
      await loadMVS(plugin, data, LoadOptions);
      gc();
      const h1 = process.memoryUsage().heapUsed;
      lines.push(`${'retained by the snapshot entries'.padEnd(34)} ${((h1 - h0) / 1e6).toFixed(1).padStart(7)} MB`);

      // how much of each frame differs from its snapshot: transforms whose params are not shared with the base
      let nFrames = 0;
      let nTransforms = 0;
      let nOwnParams = 0;
      for (const entry of plugin.managers.snapshot.state.entries) {
        const base = new Map(entry.snapshot.data?.tree.transforms.map((t) => [t.ref, t.params]));
        for (const frame of entry.snapshot.transition?.frames ?? []) {
          nFrames++;
          for (const t of frame.data.tree.transforms) {
            nTransforms++;
            if (base.get(t.ref) !== t.params) nOwnParams++;
          }
        }
      }
      lines.push(
        `${'frames'.padEnd(34)} ${nFrames} frames, ${(nTransforms / Math.max(1, nFrames)).toFixed(1)} transforms per frame, ` +
          `${(nOwnParams / Math.max(1, nFrames)).toFixed(1)} with params of their own`,
      );

      if (ProfileDir) {
        fs.mkdirSync(ProfileDir, { recursive: true });
        const session = new Session();
        session.connect();
        await session.post('Profiler.enable');
        await session.post('Profiler.setSamplingInterval', { interval: 100 });
        await session.post('Profiler.start');
        await loadMVS(plugin, data, LoadOptions);
        const { profile } = await session.post('Profiler.stop');
        fs.writeFileSync(path.join(ProfileDir, `${name}.cpuprofile`), JSON.stringify(profile));
        session.disconnect();
      }

      plugin.dispose();
      console.log(lines.join('\n'));
    });
  }
});
