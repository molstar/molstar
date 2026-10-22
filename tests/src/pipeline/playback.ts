/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Headless bench of playback: what one animation frame costs once a scene is loaded. A real PluginContext (the default
// spec and the MVS behavior) with the stub canvas of many.ts (canvas-stub.ts), driven tick by tick on a synthetic
// 60 Hz clock through `managers.animation.tick(t, true)`: every tick awaits its frame, so the per-tick numbers are
// what the browser's animation loop would pay when it keeps up. Plan and results: .v6/animation-playback/README.md.
//
// usage (through scripts/pipeline-bench.mjs, which bundles this file):
//   node scripts/pipeline-bench.mjs playback [--runs 2]                        # every scene below
//   node scripts/pipeline-bench.mjs playback --only traj:* --ticks 120
//   node scripts/pipeline-bench.mjs playback --mvsx <story.mvsx> [--snapshot i]  # also an MVS story (local file;
//                                                                               the first snapshot with a transition)
//   node scripts/pipeline-bench.mjs playback --profile <dir> [--alloc]         # one profile per scene (played range)
//
// Scenes: snapshots:plain (1cbs, a transition of 3 states × 10 with a camera move), traj:nmr (1grm, 10 models, the
// mmCIF path), traj:synth (7qpd topology and --frames jittered coordinate frames, the topology + coordinates path;
// --traj-source <cif> for another), each with a representation set: preset (the default preset: components,
// cartoon and ball-and-stick), cartoon, ball-and-stick, gaussian-surface; +focus and +select add a focused residue
// and a selected chain, +overpaint an overpainted chain (a bundle layer on every component). mvs:<file> plays one snapshot's transition of a local MVSX story.
//
// Per scene: ticks played, frames applied (ticks whose state changed), ms per tick (median and max, and the median of
// the applied ticks), KB allocated per tick over the range (V8 total_allocated_bytes, Node 25+; read once per range,
// since a read itself allocates ~165 KB), the state events per tick (cells created, updated in place, recreated,
// removed), camera resets and restarted transitions, setProps on the canvas, log entries added, marks, whole-tree
// `setSnapshot` calls, the snapshot manager's `changed` emits, and the KB of render-object arrays whose value cell
// changed per applied tick (what the renderer would upload again; read outside the timed part).

import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { Session } from 'node:inspector/promises';
import path from 'node:path';
import v8 from 'node:v8';
import { Task } from '@molstar/core/task';
import { Coordinates, Time } from '@molstar/model/model/structure/coordinates';
import { type Model, StructureElement, Unit, type Structure } from '@molstar/model/model/structure';
import { OrderedSet } from '@molstar/core/data/int';
import { Color } from '@molstar/core/util/color';
import { setStructureOverpaint } from '@molstar/plugin/state/helpers/structure-overpaint';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import { MolViewSpec } from '@molstar/mvs/behavior';
import { loadMVSX } from '@molstar/mvs/components/formats';
import { loadMVS } from '@molstar/mvs/load';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { PluginSpec } from '@molstar/plugin/spec';
import { PluginStateObject as SO, PluginStateTransform } from '@molstar/plugin/state/objects';
import { AnimateModelIndex } from '@molstar/plugin/state/animation/built-in/model-index';
import { AnimateStateSnapshotTransition } from '@molstar/plugin/state/animation/built-in/state-snapshots';
import { PluginStateSnapshotManager } from '@molstar/plugin/state/manager/snapshots';
import { PluginState } from '@molstar/plugin/state';
import { TrajectoryFromModelAndCoordinates } from '@molstar/plugin/state/transforms/structure/hierarchy';
import { graphicsStub, stubCanvas } from './canvas-stub.js';
import { ROOT } from './common.js';
import { gc } from './stage.js';

const args = process.argv.slice(2);
const arg = (n: string, d?: string) => {
  const i = args.indexOf(`--${n}`);
  return i >= 0 ? args[i + 1] : d;
};
const list = (n: string, d = '') => (arg(n, d) || '').split(',').filter(Boolean);
const runs = +arg('runs', '1')!;
const ticksArg = arg('ticks');
const profileDir = arg('profile');
const allocMode = args.includes('--alloc');
const meshMode = args.includes('--mesh');
const nFrames = +arg('frames', '20')!;
const trajSource = arg('traj-source', path.join(ROOT, 'data/examples/7qpd.fw2.cif'))!;
const mvsx = arg('mvsx');
const mvsSnapshotArg = arg('snapshot');
const only = list('only');
const skip = list('skip');
const TickMs = 1000 / 60;

const matches = (name: string, pattern: string) =>
  pattern.endsWith('*') || pattern.endsWith(':') ? name.startsWith(pattern.replace(/\*$/, '')) : name === pattern;
const enabled = (name: string) =>
  !skip.some((p) => matches(name, p)) && (!only.length || only.some((p) => matches(name, p)));

const hasAllocCounter = typeof (v8.getHeapStatistics() as any).total_allocated_bytes === 'number';
if (!hasAllocCounter) console.error('warning: Node 25+ is needed for the KB column (v8 total_allocated_bytes)');
const allocated = () => (hasAllocCounter ? ((v8.getHeapStatistics() as any).total_allocated_bytes as number) : NaN);

const graphics = graphicsStub(meshMode);

type Reprs = 'preset' | 'cartoon' | 'ball-and-stick' | 'gaussian-surface';
type Scene = { name: string; reprs?: Reprs; focus?: boolean; select?: boolean; overpaint?: boolean };

async function createPlugin() {
  const spec = DefaultPluginSpec();
  spec.behaviors = [...(spec.behaviors ?? []), PluginSpec.Behavior(MolViewSpec)];
  const plugin = new PluginContext(spec);
  await plugin.init();
  const c = stubCanvas(graphics);
  (plugin as any).canvas3d = c.canvas;
  return { plugin, c };
}

/** a residue of the first polymer unit */
function residueLoci(s: Structure, k = 0) {
  const unit = s.units.find((u) => Unit.isAtomic(u) && u.polymerElements.length > 0) ?? s.units[0];
  const { index, offsets } = unit.model.atomicHierarchy.residueAtomSegments;
  const rI = index[unit.elements[0]] + k;
  const start = OrderedSet.findPredecessorIndex(unit.elements as any, offsets[rI]);
  const end = OrderedSet.findPredecessorIndex(unit.elements as any, offsets[rI + 1]);
  return StructureElement.Loci(s, [{ unit, indices: OrderedSet.ofBounds(start, Math.max(start + 1, end)) as any }]);
}
const unitLoci = (s: Structure) => {
  const unit = s.units.find((u) => Unit.isAtomic(u) && u.polymerElements.length > 0) ?? s.units[0];
  return StructureElement.Loci(s, [{ unit, indices: OrderedSet.ofBounds(0, unit.elements.length) as any }]);
};

async function addRepresentations(plugin: PluginContext, structureCell: any, reprs: Reprs) {
  if (reprs === 'preset') {
    await plugin.builders.structure.representation.applyPreset(structureCell, 'auto');
    return;
  }
  const component = await plugin.builders.structure.tryCreateComponentStatic(structureCell, 'all');
  if (component) await plugin.builders.structure.representation.addRepresentation(component, { type: reprs });
}

async function addInteractive(plugin: PluginContext, scene: Scene) {
  const s = plugin.managers.structure.hierarchy.current.structures[0]?.cell.obj?.data as Structure | undefined;
  if (!s) return;
  if (scene.select) plugin.managers.interactivity.lociSelects.select({ loci: unitLoci(s) });
  if (scene.focus) plugin.managers.structure.focus.setFromLoci(residueLoci(s, 3));
  if (scene.overpaint) {
    const components = plugin.managers.structure.hierarchy.current.structures[0].components;
    await setStructureOverpaint(plugin, components, Color(0xff0000), async (s) => unitLoci(s));
  }
}

/** an MD-like trajectory: the topology of `source` and `n` frames of its coordinates, each jittered by up to 0.5 Å */
function syntheticCoordinates(model: Model, n: number) {
  const { x, y, z } = model.atomicConformation;
  const count = model.atomicHierarchy.atoms._rowCount;
  let seed = 7;
  const rand = () => ((seed = (seed * 16807) % 2147483647) / 2147483647 - 0.5) * 1.0;
  const frames = [];
  for (let f = 0; f < n; f++) {
    const fx = new Float32Array(count);
    const fy = new Float32Array(count);
    const fz = new Float32Array(count);
    for (let i = 0; i < count; i++) {
      fx[i] = x[i] + rand();
      fy[i] = y[i] + rand();
      fz[i] = z[i] + rand();
    }
    frames.push({
      elementCount: count,
      time: Time(f, 'step'),
      x: fx,
      y: fy,
      z: fz,
      xyzOrdering: { isIdentity: true },
    });
  }
  return Coordinates.create(frames, Time(1, 'step'), Time(0, 'step'));
}

const BenchCoordinates = PluginStateTransform.BuiltIn({
  name: 'bench-synthetic-coordinates',
  display: { name: 'Synthetic coordinates' },
  from: SO.Root,
  to: SO.Molecule.Coordinates,
  params: { key: (await import('@molstar/core/util/param-definition')).ParamDefinition.Text('') },
})({
  apply({ params }) {
    const coords = syntheticStore.get(params.key)!;
    return new SO.Molecule.Coordinates(coords, { label: 'Synthetic', description: `${coords.frames.length} frames` });
  },
});
const syntheticStore = new Map<string, Coordinates>();

async function loadTrajectory(plugin: PluginContext, scene: Scene) {
  const B = plugin.builders;
  if (scene.name.startsWith('traj:nmr')) {
    const text = readFileSync(path.join(ROOT, 'data/examples/1grm_updated.cif'), 'utf8');
    const data = await B.data.rawData({ data: text, label: '1grm' });
    const trajectory = await B.structure.parseTrajectory(data, 'mmcif');
    await B.structure.hierarchy.applyPreset(trajectory, 'default', { representationPreset: 'empty' as any });
  } else {
    const binary = trajSource.endsWith('.bcif');
    const raw = binary ? new Uint8Array(readFileSync(trajSource)) : readFileSync(trajSource, 'utf8');
    const data = await B.data.rawData({ data: raw as any, label: path.basename(trajSource) });
    const topology = await B.structure.parseTrajectory(data, 'mmcif');
    const topModel = await B.structure.createModel(topology);
    const key = `${trajSource}:${nFrames}`;
    if (!syntheticStore.has(key)) syntheticStore.set(key, syntheticCoordinates(topModel.data!, nFrames));
    const coords = await plugin.build().toRoot().apply(BenchCoordinates, { key }).commit();
    const trajectory = await plugin
      .build()
      .toRoot()
      .apply(
        TrajectoryFromModelAndCoordinates,
        { modelRef: topModel.ref, coordinatesRef: coords.ref },
        { dependsOn: [topModel.ref, coords.ref] },
      )
      .commit();
    await B.structure.hierarchy.applyPreset(trajectory, 'default', { representationPreset: 'empty' as any });
  }
  const ref = plugin.managers.structure.hierarchy.current.structures[0];
  const structure = ref.properties?.cell ?? ref.cell;
  await addRepresentations(plugin, structure, scene.reprs ?? 'preset');
  await addInteractive(plugin, scene);
}

/** 1cbs with three looks (the preset, chain-id colors, element colors), played as a transition of 30 frames of 100 ms */
async function loadPlainSnapshots(plugin: PluginContext, c: ReturnType<typeof stubCanvas>) {
  const B = plugin.builders;
  const data = await B.data.rawData({
    data: new Uint8Array(readFileSync(path.join(ROOT, 'data/examples/1cbs_full.bcif'))),
    label: '1cbs',
  });
  const trajectory = await B.structure.parseTrajectory(data, 'mmcif');
  await B.structure.hierarchy.applyPreset(trajectory, 'default');
  const components = () => plugin.managers.structure.hierarchy.current.structures[0].components;
  const looks: PluginState.Snapshot['data'][] = [];
  looks.push(plugin.state.data.getSnapshot());
  await plugin.managers.structure.component.updateRepresentationsTheme(components(), { color: 'chain-id' });
  looks.push(plugin.state.data.getSnapshot());
  await plugin.managers.structure.component.updateRepresentationsTheme(components(), { color: 'element-symbol' });
  looks.push(plugin.state.data.getSnapshot());
  const camera = c.canvas.camera.getSnapshot();
  const frames = Array.from({ length: 30 }, (_, i) => ({ durationInMs: 100, data: looks[Math.floor(i / 10)]! }));
  const snapshot: PluginState.Snapshot = {
    id: 'plain' as any,
    data: looks[0],
    camera: {
      current: { ...camera, position: [camera.position[0] + 20, camera.position[1], camera.position[2]] as any },
      transitionStyle: 'animate',
      transitionDurationInMs: 400,
    },
    transition: { frames },
  };
  plugin.managers.snapshot.add(PluginStateSnapshotManager.Entry(snapshot, { name: 'plain' }));
  plugin.managers.snapshot.setCurrent('plain');
  return { durationMs: 3000 };
}

async function loadMvs(plugin: PluginContext) {
  const data = new Uint8Array(readFileSync(mvsx!));
  const parsed = await plugin.runTask(Task.create('load', (ctx) => loadMVSX(plugin, ctx, data)));
  await loadMVS(plugin, parsed.mvsData as MVSData, { sourceUrl: parsed.sourceUrl, keepCamera: true });
  await plugin.managers.animation.stop();
  // the snapshot's start without its autoplay: the bench starts the transition itself
  const entries = plugin.managers.snapshot.state.entries;
  const index = mvsSnapshotArg ? +mvsSnapshotArg : entries.findIndex((e) => !!e.snapshot.transition);
  const entry = entries.get(index)!;
  plugin.managers.snapshot.setCurrent(entry.snapshot.id);
  const { transition } = entry.snapshot;
  await plugin.state.setSnapshot({ ...entry.snapshot, transition: transition && { ...transition, autoplay: false } });
  // the checks a snapshot apply runs over the snapshot and all its frames (S4): the median of 5 calls
  const checks: number[] = [];
  for (let k = 0; k < 5; k++) {
    const t0 = performance.now();
    PluginState.validateSnapshotTransformers(entry.snapshot);
    plugin.state.reportUnregisteredNames(entry.snapshot);
    checks.push(performance.now() - t0);
  }
  const frames = entry.snapshot.transition?.frames ?? [];
  return {
    durationMs: frames.reduce((a, f) => a + f.durationInMs, 0),
    frames: frames.length,
    index,
    checksMs: median(checks),
  };
}

type Tick = {
  ms: number;
  applied: boolean;
  created: number;
  updated: number;
  recreated: number;
  removed: number;
  changedKB: number;
};
type SceneResult = {
  scene: string;
  info: string;
  ticks: number;
  applied: number;
  msMed: number;
  msMax: number;
  msAppliedMed: number;
  kbPerTick: number;
  mbTotal: number;
  createdPerTick: number;
  updatedPerTick: number;
  recreatedPerTick: number;
  removedPerTick: number;
  cameraReset: number;
  cameraRestarts: number;
  setProps: number;
  logEntries: number;
  marks: number;
  /** whole-tree `State.setSnapshot` calls */
  setSnapshot: number;
  /** KB of render-object arrays whose value cell changed, per applied tick (median): what would be uploaded again */
  changedKBApplied: number;
  /** `managers.snapshot.events.changed` emits */
  snapshotChanged: number;
};

const median = (xs: number[]) => {
  if (!xs.length) return NaN;
  const a = [...xs].sort((x, y) => x - y);
  return a.length % 2 ? a[(a.length - 1) / 2] : (a[a.length / 2 - 1] + a[a.length / 2]) / 2;
};

async function playScene(scene: Scene, session?: Session): Promise<SceneResult | undefined> {
  const { plugin, c } = await createPlugin();
  let info = '';
  let ticks = ticksArg ? +ticksArg : 240;
  if (scene.name === 'snapshots:plain') {
    const { durationMs } = await loadPlainSnapshots(plugin, c);
    ticks = ticksArg ? +ticksArg : Math.ceil(durationMs / TickMs) + 2;
    info = '1cbs, 30 frames × 100 ms (3 looks), camera move 400 ms';
  } else if (scene.name.startsWith('mvs:')) {
    const { durationMs, frames, index, checksMs } = await loadMvs(plugin);
    ticks = ticksArg ? +ticksArg : Math.min(600, Math.ceil(durationMs / TickMs) + 2);
    info = `${path.basename(mvsx!)} snapshot ${index}: ${frames} frames, ${(durationMs / 1000).toFixed(1)} s; snapshot checks ${checksMs.toFixed(2)} ms`;
  } else {
    await loadTrajectory(plugin, scene);
    const s = plugin.managers.structure.hierarchy.current.structures[0].cell.obj!.data as Structure;
    const trajs = plugin.managers.structure.hierarchy.current.trajectories;
    const traj = trajs[trajs.length - 1].cell.obj!.data;
    // two loops at 30 frames/s: a frame change every other tick, as at the default target fps on a 60 Hz display
    ticks = ticksArg ? +ticksArg : Math.ceil((2 * traj.frameCount * 1000) / 30 / TickMs);
    info = `${s.elementCount} atoms, ${traj.frameCount} frames, ${plugin.state.data.cells.size} cells`;
  }
  await gc();

  const cellEvents = { created: 0, updated: 0, recreated: 0, removed: 0, setSnapshot: 0, snapshotChanged: 0 };
  const subs = [
    plugin.state.data.events.object.created.subscribe(() => cellEvents.created++),
    plugin.state.data.events.object.updated.subscribe((e) =>
      e.action === 'recreate' ? cellEvents.recreated++ : cellEvents.updated++,
    ),
    plugin.managers.snapshot.events.changed.subscribe(() => cellEvents.snapshotChanged++),
    plugin.state.data.events.object.removed.subscribe(() => cellEvents.removed++),
  ];
  const counts0 = { ...c.counts };
  const setSnapshot = plugin.state.data.setSnapshot.bind(plugin.state.data);
  plugin.state.data.setSnapshot = ((...a: Parameters<typeof setSnapshot>) => {
    cellEvents.setSnapshot++;
    return setSnapshot(...a);
  }) as typeof setSnapshot;
  const log0 = plugin.log.entries.size;

  c.clock.now = 0;
  if (scene.name.startsWith('traj:')) {
    await plugin.managers.animation.play(AnimateModelIndex, {
      mode: { name: 'loop', params: { direction: 'forward' } },
      duration: { name: 'computed', params: { targetFps: 30 } },
    } as any);
  } else {
    await plugin.managers.animation.play(AnimateStateSnapshotTransition, {});
  }

  if (session && allocMode) {
    await session.post('HeapProfiler.startSampling', {
      samplingInterval: 16 * 1024,
      includeObjectsCollectedByMajorGC: true,
      includeObjectsCollectedByMinorGC: true,
    });
  } else if (session) await session.post('Profiler.start');
  const perTick: Tick[] = [];
  c.changedBytes();
  // one read per range: a read of the counter itself allocates ~165 KB (v8.getHeapStatistics)
  const alloc0 = allocated();
  for (let k = 0; k < ticks; k++) {
    const t = k * TickMs;
    c.clock.now = t;
    const e0 = { ...cellEvents };
    const t0 = performance.now();
    await plugin.managers.animation.tick(t, true);
    const t1 = performance.now();
    const created = cellEvents.created - e0.created;
    const updated = cellEvents.updated - e0.updated;
    const recreated = cellEvents.recreated - e0.recreated;
    const removed = cellEvents.removed - e0.removed;
    perTick.push({
      ms: t1 - t0,
      applied: created + updated + recreated + removed > 0,
      created,
      updated,
      recreated,
      removed,
      // outside the timed part (a walk over every render object's values)
      changedKB: c.changedBytes() / 1024,
    });
    if (!plugin.managers.animation.isAnimating) break;
  }
  const allocMB = (allocated() - alloc0) / 2 ** 20;
  let profile: any;
  if (session) {
    profile = allocMode
      ? (await session.post('HeapProfiler.stopSampling')).profile
      : (await session.post('Profiler.stop')).profile;
  }
  await plugin.managers.animation.stop();
  for (const s of subs) s.unsubscribe();

  const n = perTick.length;
  const r: SceneResult = {
    scene: scene.name,
    info,
    ticks: n,
    applied: perTick.filter((x) => x.applied).length,
    msMed: +median(perTick.map((x) => x.ms)).toFixed(2),
    msMax: +Math.max(...perTick.map((x) => x.ms)).toFixed(1),
    msAppliedMed: +median(perTick.filter((x) => x.applied).map((x) => x.ms)).toFixed(2),
    kbPerTick: +((allocMB * 1024) / n).toFixed(0),
    mbTotal: +allocMB.toFixed(1),
    createdPerTick: +(perTick.reduce((a, x) => a + x.created, 0) / n).toFixed(2),
    updatedPerTick: +(perTick.reduce((a, x) => a + x.updated, 0) / n).toFixed(2),
    recreatedPerTick: +(perTick.reduce((a, x) => a + x.recreated, 0) / n).toFixed(2),
    removedPerTick: +(perTick.reduce((a, x) => a + x.removed, 0) / n).toFixed(2),
    cameraReset: c.counts.cameraReset - counts0.cameraReset,
    cameraRestarts: c.counts.cameraRestarts - counts0.cameraRestarts,
    setProps: c.counts.setProps - counts0.setProps,
    logEntries: plugin.log.entries.size - log0,
    marks: c.counts.mark - counts0.mark,
    setSnapshot: cellEvents.setSnapshot,
    changedKBApplied: +median(perTick.filter((x) => x.applied).map((x) => x.changedKB)).toFixed(0),
    snapshotChanged: cellEvents.snapshotChanged,
  };
  if (profileDir && profile) {
    writeFileSync(
      path.join(profileDir, `${scene.name.replace(/[:/+]/g, '_')}.${allocMode ? 'heapprofile' : 'cpuprofile'}`),
      JSON.stringify(profile),
    );
  }
  plugin.dispose();
  console.error(
    `${r.scene.padEnd(34)} ${String(r.ticks).padStart(4)} ticks ${String(r.applied).padStart(4)} applied ` +
      `${r.msMed.toFixed(2).padStart(7)} ms med ${r.msMax.toFixed(1).padStart(7)} max ${r.msAppliedMed.toFixed(2).padStart(7)} ms/applied ` +
      `${String(r.kbPerTick).padStart(7)} KB/tick ${r.mbTotal.toFixed(1).padStart(7)} MB ` +
      `cells c/u/R/r ${r.createdPerTick}/${r.updatedPerTick}/${r.recreatedPerTick}/${r.removedPerTick} ` +
      `cam ${r.cameraReset}/${r.cameraRestarts} setProps ${r.setProps} log ${r.logEntries} marks ${r.marks} ` +
      `setSnapshot ${r.setSnapshot} changed ${r.snapshotChanged} arrays ${r.changedKBApplied} KB/applied  ${r.info}`,
  );
  return r;
}

const scenes: Scene[] = [
  { name: 'snapshots:plain' },
  ...(mvsx ? [{ name: `mvs:${path.basename(path.dirname(mvsx))}/${path.basename(mvsx, '.mvsx')}` }] : []),
  { name: 'traj:nmr:preset', reprs: 'preset' },
  { name: 'traj:nmr:preset+focus+select', reprs: 'preset', focus: true, select: true },
  { name: 'traj:nmr:preset+focus', reprs: 'preset', focus: true },
  { name: 'traj:synth:preset', reprs: 'preset' },
  { name: 'traj:synth:preset+focus+select', reprs: 'preset', focus: true, select: true },
  { name: 'traj:synth:preset+overpaint', reprs: 'preset', overpaint: true },
  { name: 'traj:synth:cartoon', reprs: 'cartoon' },
  { name: 'traj:synth:ball-and-stick', reprs: 'ball-and-stick' },
  { name: 'traj:synth:gaussian-surface', reprs: 'gaussian-surface' },
];

let session: Session | undefined;
if (profileDir) {
  mkdirSync(profileDir, { recursive: true });
  session = new Session();
  session.connect();
  if (allocMode) await session.post('HeapProfiler.enable');
  else {
    await session.post('Profiler.enable');
    await session.post('Profiler.setSamplingInterval', { interval: 100 });
  }
}

const all: Record<string, SceneResult[]> = {};
for (let i = 0; i < runs; i++) {
  if (runs > 1) console.error(`# run ${i + 1}`);
  for (const scene of scenes) {
    if (!enabled(scene.name)) continue;
    const r = await playScene(scene, session);
    if (r) (all[scene.name] ??= []).push(r);
  }
}
// the run with the median per-tick time of each scene
const rows = Object.values(all).map((rs) => [...rs].sort((a, b) => a.msMed - b.msMed)[Math.floor((rs.length - 1) / 2)]);
const out = { bench: 'playback', mesh: meshMode, runs, frames: nFrames, rows };
const outFile = arg('out');
if (outFile) writeFileSync(outFile, JSON.stringify(out, null, 1));
else console.log(JSON.stringify(out));
