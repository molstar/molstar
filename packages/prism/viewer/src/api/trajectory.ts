/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { PluginContext } from '@molstar/plugin/context';
import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { AnimateModelIndex } from '@molstar/plugin/state/animation/built-in/model-index';
import { ModelFromTrajectory } from '@molstar/plugin/state/transforms/structure/hierarchy';
import { type Periodic, PolymerCenter } from '@molstar/plugin/state/helpers/periodic';
import {
  isPlayingAll,
  isPlayingFrames,
  periodicInfos,
  type TrajectoryInfo,
  trajectoryInfos,
  seekWhileViewing,
} from '../sources/trajectory.js';
import { ApiArea, type ApiContext, ApiResult } from './area.js';
import { StructureParam, OptionalParam } from './params.js';
import { isFailure, labelOf, loadedStructures, oneStructure } from './scope.js';

// `trajectory` (api.md §2): the models of a structure built from a trajectory (an NMR ensemble, a simulation).
// Its frame is the model index of the structure's model. Moving through frames is viewing, as the classic's arrows are:
// no history step (sequence.md §7.3), and playing is the plugin's animation of the model index. A structure with one
// model has no frames: it is not listed, and the operations say `not-found` when no structure has any.
// `trajectory.periodic` says how a simulation's box is shown: molecules whole (the
// default), in the box (centered on the polymers or not), without jumps between frames, or as in the file. A setting,
// so a history step. Operations: get, frame, play, stop, periodic.

/** What `trajectory.get` gives back. */
export interface TrajectoryValue {
  readonly structures: readonly TrajectoryInfo[];
  readonly playing: boolean;
}

/** Whether any structure has several models (the operations' `when`: the command menu greys them out otherwise). */
const hasFrames = (ctx: { plugin: PluginContext }): true | string =>
  trajectoryInfos(ctx.plugin).length ? true : 'No structure has several models.';

/** Whether any structure's frames are in a simulation box (`trajectory.periodic`'s `when`). */
const hasBox = (ctx: { plugin: PluginContext }): true | string =>
  periodicInfos(ctx.plugin).length ? true : 'No structure is in a simulation box.';

const PeriodicModes: [Periodic, string][] = [
  ['whole', 'Whole'],
  ['box', 'In the box'],
  ['nojump', 'No jumps'],
  ['off', 'As in the file'],
];
const PeriodicWords: Record<Periodic, string> = {
  whole: 'its molecules whole',
  box: 'its molecules in the box',
  nojump: 'its molecules without jumps between frames',
  off: 'its coordinates as in the file',
};

const noFrames = () =>
  ApiResult.fail('not-found', 'No structure has several models.', 'Load an ensemble or a trajectory first.');

export class TrajectoryApi extends ApiArea {
  readonly id = 'trajectory';

  readonly get = this.operation({
    id: 'trajectory.get',
    title: 'Describe the trajectories',
    params: {},
    read: true,
    async run(_, ctx) {
      const structures = trajectoryInfos(ctx.plugin);
      const playing = isPlayingFrames(ctx.plugin);
      if (!structures.length)
        return ApiResult.ok('No structure has several models.', {}, { value: { structures, playing } });
      const lines = structures.map((s) => `${s.label}: model ${s.index + 1} of ${s.count}`);
      return ApiResult.ok(`${lines.join('\n')}${playing ? '\n(playing)' : ''}`, {}, { value: { structures, playing } });
    },
  });

  readonly frame = this.operation({
    id: 'trajectory.frame',
    // One operation per state (api.md §1 "Words"): the model a structure shows; `structure.model` was another.
    aliases: ['structure.model'],
    title: 'Show a frame',
    params: {
      structure: StructureParam('A structure by label or state ref; empty: the one with frames in view'),
      index: OptionalParam(PD.Numeric(0, { min: 0, step: 1 }, { label: 'Model', description: 'The frame, from 0' })),
      by: OptionalParam(
        PD.Numeric(
          1,
          { step: 1 },
          { label: 'Move by', description: 'Frames to move by (negative: back), wrapping around' },
        ),
      ),
    },
    when: hasFrames,
    // The classic's arrows are viewing (sequence.md §7.3): the frame is not a step, and Undo leaves it where it is.
    history: false,
    // A scrub or a held arrow seeks once a frame: not a record entry each (journal.md §1, as `snapshot.seek`).
    record: false,
    async run(v, ctx) {
      const { structure } = v;
      const { index, by } = v as { index?: number; by?: number };
      if (index !== undefined && by !== undefined)
        return ApiResult.fail('invalid', 'Give an `index` or `by`, not both.');
      if (index === undefined && by === undefined) {
        return ApiResult.fail('invalid', 'Give the frame to show (`index`, from 0) or how far to move (`by`).');
      }
      const found = pickFramed(ctx, structure);
      if (isFailure(found)) return found;
      const s = trajectoryInfos(ctx.plugin).find((t) => t.ref === found.cell.transform.ref);
      if (!s) return ApiResult.fail('blocked', `${labelOf(ctx, found)} has one model, no frames to move through.`);
      if (index !== undefined && index >= s.count) {
        return ApiResult.fail('invalid', `${s.label} has ${s.count} models.`, `Give 0 to ${s.count - 1}.`);
      }
      // As `UpdateTrajectory` moves (the classic's arrows): past the end is the start.
      const to = index !== undefined ? index : (((s.index + by!) % s.count) + s.count) % s.count;
      if (to === s.index) return ApiResult.already(`${s.label} already shows model ${to + 1} of ${s.count}.`);
      // The model cell, not `UpdateTrajectory`: that moves every trajectory in the scene.
      const cell = ctx.plugin.state.data.cells.get(s.model)!;
      await seekWhileViewing(ctx.plugin, s.model, { ...(cell.transform.params as object), modelIndex: to });
      return ApiResult.ok(
        `${s.label} shows model ${to + 1} of ${s.count}.`,
        { nodes: [s.ref] },
        { value: { ...s, index: to } },
      );
    },
  });

  readonly play = this.operation({
    id: 'trajectory.play',
    title: 'Play the frames',
    params: {
      structure: StructureParam('A structure by label or state ref; empty: every structure with frames'),
      mode: PD.Select<'loop' | 'palindrome' | 'once'>(
        'loop',
        PD.arrayToOptions(['loop', 'palindrome', 'once'] as const),
        { label: 'Mode', description: 'Loop, back and forth (palindrome), or once through' },
      ),
      seconds: PD.Numeric(
        5,
        { min: 1, max: 120, step: 0.5 },
        { label: 'Seconds per pass', description: 'Seconds for a pass through the frames' },
      ),
    },
    when: hasFrames,
    // The animation steps the model index of one structure's model (`structure`), or of every trajectory in the scene
    // (empty, as the classic's): viewing, no step. One plays at a time: playing another stops the one before.
    history: false,
    // The plugin's animation loop runs on the canvas's frames.
    modules: ['canvas'],
    async run({ structure, mode, seconds }, ctx) {
      const all = trajectoryInfos(ctx.plugin);
      if (!all.length) return noFrames();
      let one: TrajectoryInfo | undefined;
      if (structure) {
        const found = oneStructure(ctx, structure);
        if (isFailure(found)) return found;
        one = all.find((t) => t.ref === found.cell.transform.ref);
        if (!one) return ApiResult.fail('blocked', `${labelOf(ctx, found)} has one model, no frames to play.`);
        if (isPlayingFrames(ctx.plugin, one.model)) return ApiResult.already(`${one.label}'s ${one.unit} are playing.`);
      } else if (isPlayingAll(ctx.plugin)) return ApiResult.already('The models are playing.');
      await ctx.plugin.managers.animation.play(AnimateModelIndex, {
        mode:
          mode === 'palindrome'
            ? { name: 'palindrome' as const, params: {} }
            : { name: mode, params: { direction: 'forward' as const } },
        duration: { name: 'fixed', params: { durationInS: seconds } },
        refs: one ? [one.model] : [],
      });
      return ApiResult.ok(one ? `Playing ${one.label}'s ${one.unit}.` : 'Playing the models.', {
        view: true,
        ...(one ? { nodes: [one.ref] } : {}),
      });
    },
  });

  readonly periodic = this.operation({
    id: 'trajectory.periodic',
    title: 'Show the periodic box',
    params: {
      structure: StructureParam('A structure by label or state ref; empty: the one in a simulation box'),
      mode: PD.Select<Periodic>('whole', PeriodicModes, {
        label: 'Molecules',
        description:
          'whole: each molecule split across the box joined again; box: also moved into the box; nojump: kept from jumping a box length between frames (they leave the box as they diffuse); off: as in the file',
      }),
      center: OptionalParam(
        PD.Select<'box' | 'polymers'>(
          'box',
          [
            ['box', 'The box'],
            ['polymers', 'The polymers'],
          ],
          { label: 'Centered on', description: 'With box: the box as it is, or the protein (polymers) at its center' },
        ),
      ),
    },
    when: hasBox,
    async run(v, ctx) {
      const { structure, mode } = v;
      const { center } = v as { center?: 'box' | 'polymers' };
      const all = periodicInfos(ctx.plugin);
      if (!all.length) return ApiResult.fail('not-found', 'No structure is in a simulation box.');
      let info = all.length === 1 && !structure ? all[0] : undefined;
      if (!info) {
        if (!structure) {
          return ApiResult.fail(
            'ambiguous',
            `Which structure? ${all.length} are in a simulation box.`,
            `Name one: ${all.map((s) => `${s.label} (${s.ref})`).join(', ')}.`,
          );
        }
        const found = oneStructure(ctx, structure);
        if (isFailure(found)) return found;
        info = all.find((s) => s.ref === found.cell.transform.ref);
        if (!info) return ApiResult.fail('blocked', `${labelOf(ctx, found)} is not in a simulation box.`);
      }
      const current = info.params;
      const currentCenter = current.periodicCenter?.name === 'selection' ? 'polymers' : 'box';
      const nextCenter = center ?? currentCenter;
      if ((current.periodic ?? 'whole') === mode && nextCenter === currentCenter) {
        return ApiResult.already(`${info.label} already shows ${PeriodicWords[mode]}.`);
      }
      const periodicCenter =
        nextCenter === 'polymers'
          ? { name: 'selection' as const, params: { script: PolymerCenter } }
          : { name: 'box' as const, params: {} };
      await ctx.plugin
        .build()
        .to(info.model)
        .update(ModelFromTrajectory, (p) => ({ ...p, periodic: mode, periodicCenter }))
        .commit();
      const centered = mode === 'box' && nextCenter === 'polymers' ? ', centered on the polymers' : '';
      return ApiResult.ok(`${info.label} shows ${PeriodicWords[mode]}${centered}.`, { nodes: [info.ref] });
    },
  });

  readonly stop = this.operation({
    id: 'trajectory.stop',
    title: 'Stop playing the frames',
    params: {},
    history: false,
    async run(_, ctx) {
      if (!isPlayingFrames(ctx.plugin)) return ApiResult.already('The models are not playing.');
      await ctx.plugin.managers.animation.stop();
      return ApiResult.ok('Stopped playing the models.', { view: true });
    },
  });
}

/** The named structure, else the one with frames (the focused one among several), else why not. */
function pickFramed(
  ctx: ApiContext,
  name: string,
): StructureRef | ReturnType<typeof noFrames> | ReturnType<typeof oneStructure> {
  if (name) return oneStructure(ctx, name);
  const framed = new Set(trajectoryInfos(ctx.plugin).map((t) => t.ref));
  const list = loadedStructures(ctx).filter((s) => framed.has(s.cell.transform.ref));
  if (!list.length) return noFrames();
  if (list.length === 1) return list[0];
  const focus = ctx.plugin.managers.structure.focus.current?.loci.structure.root;
  const focused = focus && list.find((s) => s.cell.obj?.data.root === focus);
  if (focused) return focused;
  return ApiResult.fail(
    'ambiguous',
    `Which structure? ${list.length} have several models.`,
    `Name one: ${list.map((s) => `${labelOf(ctx, s)} (${s.cell.transform.ref})`).join(', ')}.`,
  );
}
