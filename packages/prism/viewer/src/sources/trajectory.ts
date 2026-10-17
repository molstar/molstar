/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

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
}

/** The transformer that adds coordinates (a simulation's frames) to a model. */
const FromCoordinates = 'trajectory-from-model-and-coordinates';

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
      unit: s.model.trajectory?.cell.transform.transformer.definition.name === FromCoordinates ? 'frames' : 'models',
    });
  }
  return out;
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
    () => trajectoryInfos(plugin),
    () => merge(hierarchy.behaviors.selection, plugin.state.data.events.changed),
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
    () => isPlayingFrames(plugin),
    () => plugin.managers.animation.events.updated,
    {
      coalesce: 'none',
    },
  );
  // Which structures play changes with the animation (played, stopped) and with the structures (one added or removed),
  // not with the frame shown: it does not follow `frames`, which changes every frame of a play.
  const playingRefs = watchedState(
    () => playingStructures(plugin),
    () => merge(plugin.managers.animation.events.updated, hierarchy.behaviors.selection),
    { coalesce: 'none', compare: sameSet },
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
