/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ModelPeriodicity } from '@molstar/model/model/structure/model/periodic';
import { PluginCommands } from '@molstar/plugin/commands';
import type { PeriodicProps } from '@molstar/plugin/state/helpers/periodic';
import type { PluginContext } from '@molstar/plugin/context';
import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { derive } from '@molstar/prism-ui/reactive/state';
import { merge } from 'rxjs';
import { structureNames } from './naming.js';
import { sameItems } from './state.js';
import { watchedState } from './watched.js';

// The models of a structure as the trajectory area and the trajectory bar read them (state.md §4, row 2): a structure
// built from a trajectory with several frames (an NMR ensemble, a simulation) shows one of them, its model index.

/** The name of the animation that steps through the frames (`built-in.animate-model-index`). */
export const ModelIndexAnimation = 'built-in.animate-model-index';

/** A structure with several models. */
export interface TrajectoryInfo {
  /** The structure's state ref. */
  readonly ref: string;
  /** The structure's name for people. */
  readonly label: string;
  /** The model's state ref (the cell whose `modelIndex` is the frame). */
  readonly model: string;
  readonly count: number;
  /** The frame shown, from 0. */
  readonly index: number;
  /**
   * What its frames are called: `frames` for coordinates added to a model (a simulation's DCD, XTC, …), `models` for a
   * file's own models (an NMR ensemble).
   */
  readonly unit: 'models' | 'frames';
  /**
   * A particle list read from a trajectory (a Simularium file), whose frames the canvas' bar shows as a structure's:
   * `ref` and `model` are the list's, its frames moved through `particles.*`. Absent: a structure's.
   */
  readonly kind?: 'particles';
}

/** The transformer that adds coordinates (a simulation's frames) to a model. */
// (MolViewSpec's `trajectory-with-coordinates` too: the same trajectory, built from an MVS coordinates node)
const FromCoordinates = new Set(['trajectory-from-model-and-coordinates', 'trajectory-with-coordinates']);

/** The structures that have several frames, in the order of the scene. */
export function trajectoryInfos(plugin: PluginContext): TrajectoryInfo[] {
  const list = plugin.managers.structure.hierarchy.current.structures;
  const names = structureNames(list);
  const out: TrajectoryInfo[] = [];
  for (const s of list) {
    const count = s.model?.trajectory?.cell.obj?.data.frameCount ?? 0;
    if (!s.model || count < 2) continue;
    // Read from the state, not the hierarchy's copy: the frame is the model cell's current param.
    const cell = plugin.state.data.cells.get(s.model.cell.transform.ref) ?? s.model.cell;
    out.push({
      ref: s.cell.transform.ref,
      label: names.get(s.cell.transform.ref) ?? s.cell.transform.ref,
      model: cell.transform.ref,
      count,
      index: (cell.transform.params as { modelIndex: number }).modelIndex,
      unit: FromCoordinates.has(s.model.trajectory?.cell.transform.transformer.definition.name ?? '')
        ? 'frames'
        : 'models',
    });
  }
  return out;
}

/** A structure whose frames are in a simulation box (`trajectory.periodic`). */
export interface PeriodicInfo {
  readonly ref: string;
  readonly label: string;
  /** The model's state ref (the cell whose `periodic` params say how the box is shown) */
  readonly model: string;
  /** Its `ModelFromTrajectory` params now */
  readonly params: PeriodicProps & { modelIndex: number };
  /** It shows its periodic images (the cells around: its structure's type is `symmetry`) */
  readonly images: boolean;
}

const ModelFromTrajectoryName = 'model-from-trajectory';

/** The structures whose model is a frame in a simulation box. */
export function periodicInfos(plugin: PluginContext): PeriodicInfo[] {
  const list = plugin.managers.structure.hierarchy.current.structures;
  const names = structureNames(list);
  const out: PeriodicInfo[] = [];
  for (const s of list) {
    if (!s.model) continue;
    const cell = plugin.state.data.cells.get(s.model.cell.transform.ref) ?? s.model.cell;
    if (cell.transform.transformer.definition.name !== ModelFromTrajectoryName || !cell.obj) continue;
    if (!ModelPeriodicity.isApplicable(cell.obj.data)) continue;
    out.push({
      ref: s.cell.transform.ref,
      label: names.get(s.cell.transform.ref) ?? s.cell.transform.ref,
      model: cell.transform.ref,
      params: cell.transform.params as PeriodicInfo['params'],
      images: (s.cell.transform.params as { type?: { name?: string } } | undefined)?.type?.name === 'symmetry',
    });
  }
  return out;
}

function samePeriodic(a: ReadonlyMap<string, PeriodicInfo>, b: ReadonlyMap<string, PeriodicInfo>) {
  if (a.size !== b.size) return false;
  for (const [ref, p] of a) {
    const q = b.get(ref);
    if (!q || q.label !== p.label || q.params.periodic !== p.params.periodic || q.images !== p.images) return false;
    if (q.params.periodicCenter?.name !== p.params.periodicCenter?.name) return false;
  }
  return true;
}

const ParticleListFromTrajectory = 'ms-plugin.particle-list-from-trajectory';
/** The animation that steps a particle trajectory (the particles feature's, `ParticlesRegistry`). */
export const ParticleTrajectoryAnimation = 'built-in.animate-particle-trajectory';

/** The particle lists read from a trajectory with several frames (Simularium), after the structures. */
export function particleTrajectoryInfos(plugin: PluginContext): TrajectoryInfo[] {
  const out: TrajectoryInfo[] = [];
  for (const l of plugin.managers.particles.hierarchy.current.lists) {
    const cell = plugin.state.data.cells.get(l.cell.transform.ref) ?? l.cell;
    if (cell.transform.transformer.id !== ParticleListFromTrajectory) continue;
    const parent = plugin.state.data.cells.get(cell.transform.parent)?.obj;
    const count = (parent?.data as { frameCount?: number } | undefined)?.frameCount ?? 0;
    if (count < 2) continue;
    out.push({
      ref: cell.transform.ref,
      label: cell.obj?.label ?? cell.transform.ref,
      model: cell.transform.ref,
      count,
      index: (cell.transform.params as { frameIndex?: number }).frameIndex ?? 0,
      unit: 'frames',
      kind: 'particles',
    });
  }
  return out;
}

/** Everything with frames the canvas' bar can show: the structures', then the particle trajectories'. */
export function frameInfos(plugin: PluginContext): TrajectoryInfo[] {
  return [...trajectoryInfos(plugin), ...particleTrajectoryInfos(plugin)];
}

/** A particle trajectory's frames are being played: any, or (`ref`, a list's) that list's. */
export function isPlayingParticleFrames(plugin: PluginContext, ref?: string): boolean {
  const animation = plugin.managers.animation;
  if (!animation.isAnimating || animation.current?.anim.name !== ParticleTrajectoryAnimation) return false;
  if (ref === undefined) return true;
  const refs = (animation.current.paramValues as { refs?: readonly string[] }).refs;
  return !refs?.length || refs.includes(ref);
}

/** The frames are being played: any, or (`model`, a model cell's ref) that model's. */
export function isPlayingFrames(plugin: PluginContext, model?: string): boolean {
  const animation = plugin.managers.animation;
  if (!animation.isAnimating || animation.current?.anim.name !== ModelIndexAnimation) return false;
  if (model === undefined) return true;
  const refs = (animation.current.paramValues as { refs?: readonly string[] }).refs;
  return !refs?.length || refs.includes(model);
}

/** Every structure's frames play (a play with no `structure`: the animation's `refs` empty). */
export function isPlayingAll(plugin: PluginContext): boolean {
  if (!isPlayingFrames(plugin)) return false;
  const refs = (plugin.managers.animation.current!.paramValues as { refs?: readonly string[] }).refs;
  return !refs?.length;
}

const NoneShown: ReadonlySet<string> = new Set();

/** The structures (of `list`) whose frames play: one, every one (a play of all), or none. */
export function playingStructures(
  plugin: PluginContext,
  list: readonly TrajectoryInfo[] = trajectoryInfos(plugin),
): ReadonlySet<string> {
  if (!isPlayingFrames(plugin)) return NoneShown;
  const out = new Set<string>();
  for (const t of list) if (isPlayingFrames(plugin, t.model)) out.add(t.ref);
  return out;
}

/** What plays: the structures whose frames play, or the particle list whose do. */
function playingFrames(plugin: PluginContext): ReadonlySet<string> {
  if (!isPlayingParticleFrames(plugin)) return playingStructures(plugin);
  const out = new Set<string>();
  for (const t of particleTrajectoryInfos(plugin)) if (isPlayingParticleFrames(plugin, t.ref)) out.add(t.ref);
  return out;
}

/** The structure whose frames play (the first, when every one plays), or none. */
export function playingStructure(plugin: PluginContext): string | undefined {
  return playingStructures(plugin).values().next().value;
}

const sameSet = (a: ReadonlySet<string>, b: ReadonlySet<string>) => {
  if (a.size !== b.size) return false;
  for (const x of a) if (!b.has(x)) return false;
  return true;
};

// Moving through frames is viewing, not editing (sequence.md §7.3): the operations that do it (`trajectory.frame`, the
// frame animation) are not history steps, so the history's watch for changes made outside the API leaves them alone.
const viewing = new WeakMap<PluginContext, number>();

/** Runs `work`, whose changes to the data tree are viewing (see above). */
export async function whileViewing<T>(plugin: PluginContext, work: () => Promise<T>): Promise<T> {
  viewing.set(plugin, (viewing.get(plugin) ?? 0) + 1);
  try {
    return await work();
  } finally {
    const n = (viewing.get(plugin) ?? 1) - 1;
    if (n > 0) viewing.set(plugin, n);
    else viewing.delete(plugin);
  }
}

/**
 * Sets the params of the cell `ref` (a frame to show) as viewing (`whileViewing`), with no "Updated X in N ms" log
 * line per cell it updates: about 9 a seek, which a scrub would repeat until no other line was left in the Log's
 * Transforms (perf-audit-2 10/P14).
 */
export function seekWhileViewing(plugin: PluginContext, ref: string, params: object): Promise<void> {
  const state = plugin.state.data;
  const tree = state.build().to(ref).update(params);
  return whileViewing(plugin, () =>
    PluginCommands.State.Update(plugin, { state, tree, options: { doNotLogTiming: true } }),
  );
}

/** The data tree is changing as the person views it: frames stepped through, or an animation playing. */
export function isViewing(plugin: PluginContext): boolean {
  return viewing.has(plugin) || plugin.managers.animation.isAnimating;
}

/** The same structure with the same frames (any frame shown). */
const sameShape = (a: TrajectoryInfo, b: TrajectoryInfo) =>
  a.ref === b.ref && a.label === b.label && a.model === b.model && a.count === b.count && a.unit === b.unit;
const sameShapes = (a: readonly TrajectoryInfo[], b: readonly TrajectoryInfo[]) => sameItems(a, b, sameShape);

const sameInfo = (a: TrajectoryInfo, b: TrajectoryInfo) =>
  a.ref === b.ref &&
  a.label === b.label &&
  a.model === b.model &&
  a.count === b.count &&
  a.index === b.index &&
  a.unit === b.unit;

export interface TrajectorySource {
  /** Every structure with several frames. */
  readonly frames: ReadonlyState<readonly TrajectoryInfo[]>;
  /**
   * The same by structure ref, as of when the structures last changed (its `index` may be old: the frame shown is
   * `indexByRef`'s): what changes rarely (names, counts, units), apart from what changes every frame.
   */
  readonly byRef: ReadonlyState<ReadonlyMap<string, TrajectoryInfo>>;
  /** The frame each structure shows (from 0), by structure ref. */
  readonly indexByRef: ReadonlyState<ReadonlyMap<string, number>>;
  /** The one in view: the focused structure's, else the first (its `index` as `byRef`'s). */
  readonly shown: ReadonlyState<TrajectoryInfo | undefined>;
  /** The frames are playing. */
  readonly playing: ReadonlyState<boolean>;
  /** The structures whose frames play: one, every one (a play of all), or none. */
  readonly playingRefs: ReadonlyState<ReadonlySet<string>>;
  /**
   * The structures whose Structure-panel row is on screen (the row tells, while its view is active): there the row has
   * the controls, and the canvas' bar steps back (ui/shell.tsx).
   */
  readonly rowsShown: ReadonlyState<ReadonlySet<string>>;
  /** A row of `structure` is on screen (`true`) or no longer (`false`); counted, rows may be made twice. */
  markRowShown(structure: string, on: boolean): void;
  /** The structures in a simulation box, by structure ref: how their frames show it (`trajectory.periodic`). */
  readonly periodicByRef: ReadonlyState<ReadonlyMap<string, PeriodicInfo>>;
}

const sources = new WeakMap<PluginContext, TrajectorySource>();

/** The plugin's trajectory source, one for every reader (the scene's rows, the canvas' bar): one read a frame. */
export function trajectorySourceOf(plugin: PluginContext): TrajectorySource {
  let s = sources.get(plugin);
  if (!s) {
    s = createTrajectorySource(plugin);
    sources.set(plugin, s);
  }
  return s;
}

/**
 * The states listen to the plugin only while something subscribes to them (a bar that shows), and read once per
 * animation frame however many events a transaction (or a playing trajectory) fires.
 */
export function createTrajectorySource(plugin: PluginContext): TrajectorySource {
  const hierarchy = plugin.managers.structure.hierarchy;
  const frames = watchedState(
    () => frameInfos(plugin),
    () =>
      merge(
        hierarchy.behaviors.selection,
        plugin.managers.particles.hierarchy.behaviors.selection,
        plugin.state.data.events.changed,
      ),
    { compare: (a, b) => sameItems(a, b, sameInfo) },
  );
  const pick = (list: readonly TrajectoryInfo[]) => {
    const root = plugin.managers.structure.focus.current?.loci.structure.root;
    const focused = root && hierarchy.current.structures.find((s) => s.cell.obj?.data.root === root);
    return list.find((f) => f.ref === focused?.cell.transform.ref) ?? list[0];
  };
  // What changes with the structures (which have frames, their names, counts and units), apart from the frame each
  // shows (which changes every frame of a play): readers of the one do not follow the other. One check a frame here.
  const shapes = derive(frames, (list) => list, sameShapes);
  const byRef = derive(shapes, (list) => new Map(list.map((t) => [t.ref, t])));
  const indexByRef = derive(frames, (list) => new Map(list.map((t) => [t.ref, t.index])));
  // Follows the structures (rare) and the focus (rare): no second wait.
  const shown = watchedState(
    () => pick(shapes.value),
    () => merge(shapes, plugin.managers.structure.focus.behaviors.current),
    { compare: (a, b) => (a && b ? sameShape(a, b) : a === b), coalesce: 'none', derived: true },
  );
  const playing = watchedState(
    () => isPlayingFrames(plugin) || isPlayingParticleFrames(plugin),
    () => plugin.managers.animation.events.updated,
    {
      coalesce: 'none',
    },
  );
  // Which structures play changes with the animation (played, stopped) and with the structures (one added or removed),
  // not with the frame shown: it does not follow `frames`, which changes every frame of a play.
  const playingRefs = watchedState(
    () => playingFrames(plugin),
    () =>
      merge(
        plugin.managers.animation.events.updated,
        hierarchy.behaviors.selection,
        plugin.managers.particles.hierarchy.behaviors.selection,
      ),
    { coalesce: 'none', compare: sameSet },
  );
  // Changes with the structures and with a periodic param, read once a frame; equal while those stay (a playing
  // trajectory changes the frame, not them).
  const periodicByRef = watchedState(
    () => new Map(periodicInfos(plugin).map((p) => [p.ref, p])),
    () => merge(hierarchy.behaviors.selection, plugin.state.data.events.changed),
    { compare: samePeriodic },
  );
  const counts = new Map<string, number>();
  const rowsShown = observableState<ReadonlySet<string>>(NoneShown);
  const markRowShown = (structure: string, on: boolean) => {
    const n = (counts.get(structure) ?? 0) + (on ? 1 : -1);
    if (n > 0) counts.set(structure, n);
    else counts.delete(structure);
    const next = new Set(counts.keys());
    if (!sameSet(next, rowsShown.value)) rowsShown.set(next);
  };
  return {
    periodicByRef,
    frames,
    byRef,
    indexByRef,
    shown,
    playing,
    playingRefs,
    rowsShown: rowsShown.asReadonly,
    markRowShown,
  };
}
