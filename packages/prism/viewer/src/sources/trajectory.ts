/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import type { ReadonlyState } from '@molstar/prism-ui/core';
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

/** The structure whose frames play (the first, when every one plays), or none. */
export function playingStructure(plugin: PluginContext): string | undefined {
  if (!isPlayingFrames(plugin)) return undefined;
  return trajectoryInfos(plugin).find((t) => isPlayingFrames(plugin, t.model))?.ref;
}

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
  /** The one the bar follows: the focused structure's, else the first. */
  readonly shown: ReadonlyState<TrajectoryInfo | undefined>;
  /** The frames are playing. */
  readonly playing: ReadonlyState<boolean>;
  /** The structure whose frames play (the canvas' chip follows it), or none. */
  readonly playingRef: ReadonlyState<string | undefined>;
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
  // Follows `frames` (already once per frame) and the focus (rare): no second wait.
  const shown = watchedState(
    () => pick(frames.value),
    () => merge(frames, plugin.managers.structure.focus.behaviors.current),
    { compare: (a, b) => (a && b ? sameInfo(a, b) : a === b), coalesce: 'none', derived: true },
  );
  const playing = watchedState(
    () => isPlayingFrames(plugin),
    () => plugin.managers.animation.events.updated,
    {
      coalesce: 'none',
    },
  );
  const playingRef = watchedState(
    () => playingStructure(plugin),
    () => merge(plugin.managers.animation.events.updated, frames),
    { coalesce: 'none' },
  );
  return { frames, shown, playing, playingRef };
}
