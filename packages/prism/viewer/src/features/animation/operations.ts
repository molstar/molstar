/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { areAnimationPropsEqual } from '@molstar/graphics/geo/geometry/animation';
import type { PluginContext } from '@molstar/plugin/context';
import {
  clearStructureWiggle,
  setStructureWiggleFromUncertainty,
} from '@molstar/plugin/state/helpers/structure-wiggle';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import { type ApiAreaContext, ApiArea, type ApiContext, type ApiFailure, ApiResult } from '../../api/area.js';
import { OptionalParam } from '../../api/params.js';
import { checkValues } from '../../api/run.js';
import {
  type AnimationInfo,
  type AnimationState,
  animationInfos,
  createAnimationSource,
  paramsOf,
  shortName,
  uncertaintyScale,
  wiggleCells,
  wiggleOf,
  type AnimationSource,
  type WiggleMode,
  type WiggleType,
} from './source.js';

// The `animation` feature's operations (api.md §2): the plugin's animation manager as `play` (an animation by
// name and its params) and `stop`, and the structures' procedural animation as `wiggle`. Playing is a presentation, not
// an edit: no history step (sequence.md §7.3), and the history leaves the state it changes alone (sources/trajectory.ts,
// `isViewing`). Wiggle changes representation state, so it is a step, and Undo restores it.

/** What `animation.get` gives back. */
export type AnimationList = AnimationState;

type RegisteredAnimation = PluginContext['managers']['animation']['animations'][number];

/** The registered animations, as `name (Title)` for a hint. */
const listed = (list: readonly RegisteredAnimation[]) =>
  list.map((a) => `${shortName(a.name)} (${a.display.name})`).join(', ');

/**
 * The animation a name means (api.md §2): its full id exactly (`built-in.animate-camera-spin`), else the id
 * without its namespace (`animate-camera-spin`), else its display name ("Camera Spin"), the last two in any
 * case. A name several animations answer to is `ambiguous`, one none does `not-found`. No name is the one playing, else
 * the only one registered, else `invalid`: the manager's "current" is the first registered until something plays, so it
 * does not say what a person chose.
 */
export function pickAnimation(plugin: PluginContext, given: string | undefined): RegisteredAnimation | ApiFailure {
  const manager = plugin.managers.animation;
  const all = manager.animations;
  if (!all.length) return ApiResult.fail('not-found', 'No animations are registered.');
  const name = given?.trim();
  if (!name) {
    const playing = manager.isAnimating ? manager.current?.anim : undefined;
    if (playing) return playing;
    if (all.length === 1) return all[0];
    return ApiResult.fail('invalid', `Which animation? ${all.length} are registered.`, `Name one: ${listed(all)}.`);
  }
  const exact = all.find((a) => a.name === name);
  if (exact) return exact;
  const lower = name.toLowerCase();
  for (const key of [(a: RegisteredAnimation) => shortName(a.name), (a: RegisteredAnimation) => a.display.name]) {
    const found = all.filter((a) => key(a).toLowerCase() === lower);
    if (found.length === 1) return found[0];
    if (found.length > 1) {
      return ApiResult.fail(
        'ambiguous',
        `Which '${name}'? ${found.length} animations have that name.`,
        `Give its full id: ${found.map((a) => a.name).join(', ')}.`,
      );
    }
  }
  return ApiResult.fail('not-found', `There is no animation '${name}'.`, `Animations: ${listed(all)}.`);
}

/** `play`'s params: the animation by name (read when the operation is described, so it lists what is registered). */
const PlayParams = (ctx: ApiAreaContext) => ({
  name: OptionalParam(
    PD.Text('', {
      description: `The animation: its name without the namespace (animate-camera-spin), its title (Camera Spin) or its full id; absent: the one playing. Registered: ${listed(ctx.plugin.managers.animation.animations) || 'none'}`,
    }),
  ),
  params: PD.Value<Record<string, unknown>>(
    {},
    { description: "Some of that animation's params (`animation.get` lists them); the rest are its defaults" },
  ),
});

const WiggleParams = {
  mode: OptionalParam(
    PD.Select<WiggleMode>(
      'dynamics',
      [
        ['dynamics', 'Dynamics (B-factor)'],
        ['uncertainty', 'Model uncertainty'],
        ['off', 'Off'],
      ],
      {
        description: 'Move atoms by noise (dynamics), by their B-factor or RMSF (uncertainty), or stop; absent: as now',
      },
    ),
  ),
  amplitude: OptionalParam(
    PD.Numeric(1, { min: 0.01, max: 5, step: 0.01 }, { description: 'How far the atoms move; absent: as now, else 1' }),
  ),
  type: OptionalParam(
    PD.Select<WiggleType>(
      'position',
      [
        ['position', 'Position'],
        ['group', 'Group'],
      ],
      { description: 'position: nearby atoms move together; group: each atom (group) on its own' },
    ),
  ),
  speed: OptionalParam(PD.Numeric(7, { min: 0, max: 10, step: 0.1 }, { description: 'How fast the atoms move' })),
  frequency: OptionalParam(
    PD.Numeric(0.2, { min: 0.01, max: 2, step: 0.01 }, { description: 'Lower: nearby atoms move more together' }),
  ),
  tumble: OptionalParam(
    PD.Numeric(
      0,
      { min: 0, max: 10, step: 0.1 },
      { description: 'How far whole copies (instances) tumble, in Å; 0: off' },
    ),
  ),
  tumbleSpeed: OptionalParam(PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, { description: 'How fast they tumble' })),
  tumbleFrequency: OptionalParam(
    PD.Numeric(0.2, { min: 0, max: 2, step: 0.01 }, { description: 'Lower: nearby copies tumble more together' }),
  ),
};

/** The animation's params with the values given: a mapped param takes its option's defaults for what is missing. */
function animationValues(
  params: PD.Params,
  given: Record<string, unknown>,
): { values: Record<string, unknown> } | ApiFailure {
  const checked = checkValues(params, given);
  if (ApiResult.isFailure(checked)) return checked;
  for (const [key, value] of Object.entries(given)) {
    const p = params[key];
    if (!p || p.type !== 'mapped' || !value || typeof value !== 'object') continue;
    const { name, params: own } = value as { name?: unknown; params?: unknown };
    const names = p.select.options.map((o) => o[0]);
    if (typeof name !== 'string' || !names.includes(name)) {
      return ApiResult.fail('invalid', `'${key}' takes a name and params.`, `Names: ${names.join(', ')}.`);
    }
    checked.values[key] = { name, params: mergeOwn(p.map(name).defaultValue, own) };
  }
  return checked;
}

const mergeOwn = (defaults: unknown, given: unknown) =>
  defaults && typeof defaults === 'object' && !Array.isArray(defaults) && given && typeof given === 'object'
    ? { ...defaults, ...given }
    : (given ?? defaults);

export class AnimationApi extends ApiArea {
  readonly id = 'animation';

  #source: AnimationSource | undefined;

  /** The animations and what plays, as a state the views read; made on first use, released with the area. */
  get state(): ReadonlyState<AnimationState> {
    return (this.#source ??= createAnimationSource(this.ctx.plugin)).state;
  }

  /** An animation's params, by its full id (with their defaults), for the forms. */
  paramsOf(id: string): PD.Params | undefined {
    return paramsOf(this.ctx.plugin, id);
  }

  dispose(): void {
    this.#source?.dispose();
    this.#source = undefined;
    // An animation playing holds the plugin's frame loop.
    this.ctx.plugin.managers.animation.stop().catch(() => {});
  }

  readonly get = this.operation({
    id: 'animation.get',
    title: 'List the animations',
    params: {},
    read: true,
    async run(_, ctx) {
      const animations = animationInfos(ctx.plugin);
      const manager = ctx.plugin.managers.animation;
      const state: AnimationList = {
        animations,
        ...(manager.isAnimating && manager.current ? { playing: shortName(manager.current.anim.name) } : {}),
        wiggle: wiggleOf(ctx.plugin),
      };
      const line = (a: AnimationInfo) =>
        `${a.name}: ${a.title}${a.available === true ? '' : ` (not now: ${a.available})`}${state.playing === a.name ? ' (playing)' : ''}`;
      return ApiResult.ok(
        animations.length ? animations.map(line).join('\n') : 'No animations are registered.',
        {},
        { value: state },
      );
    },
  });

  readonly play = this.operation({
    id: 'animation.play',
    title: 'Play an animation',
    params: PlayParams,
    // A presentation, not an edit (sequence.md §7.3).
    history: false,
    // The plugin ticks its animations on the canvas's frames.
    modules: ['canvas'],
    async run({ name, params }, ctx) {
      const manager = ctx.plugin.managers.animation;
      const anim = pickAnimation(ctx.plugin, name);
      if (ApiResult.isFailure(anim)) return anim;
      const can = anim.canApply?.(ctx.plugin);
      if (can && !can.canApply) {
        return ApiResult.fail(
          'blocked',
          `${anim.display.name} cannot play now: ${can.reason ?? 'nothing to animate'}.`,
        );
      }
      const checked = animationValues(anim.params(ctx.plugin) as PD.Params, params);
      if (ApiResult.isFailure(checked)) return checked;
      const now = manager.current;
      if (manager.isAnimating && now?.anim === anim && same(now.paramValues, checked.values)) {
        return ApiResult.already(`${anim.display.name} is playing.`);
      }
      await manager.play(anim, checked.values);
      return ApiResult.ok(`Playing ${anim.display.name}.`, { view: true });
    },
  });

  readonly stop = this.operation({
    id: 'animation.stop',
    title: 'Stop the animation',
    params: {},
    history: false,
    async run(_, ctx) {
      const manager = ctx.plugin.managers.animation;
      if (!manager.isAnimating) return ApiResult.already('No animation is playing.');
      const name = manager.current?.anim.display.name ?? 'The animation';
      await manager.stop();
      return ApiResult.ok(`Stopped ${name}.`, { view: true });
    },
  });

  readonly wiggle = this.operation({
    id: 'animation.wiggle',
    title: 'Wiggle',
    params: WiggleParams,
    // It is representation state (the animation params of every representation, the uncertainty layers): a step.
    // Each value absent is kept as it is.
    async run({ mode, amplitude, type, speed, frequency, tumble, tumbleSpeed, tumbleFrequency }, ctx) {
      const plugin = ctx.plugin;
      const structures = plugin.managers.structure.hierarchy.current.structures;
      if (!structures.length) return ApiResult.fail('not-found', 'Nothing is loaded.', 'Load a structure first.');
      const now = wiggleOf(plugin);
      const nextMode = mode ?? now.mode;
      const nextAmplitude = amplitude ?? (nextMode === now.mode && now.amplitude > 0 ? now.amplitude : 1);
      const manager = plugin.managers.structure.component;
      const options = manager.state.options;
      const animation = {
        ...options.animation,
        ...(type !== undefined ? { wiggleMode: type } : {}),
        ...(speed !== undefined ? { wiggleSpeed: speed } : {}),
        ...(frequency !== undefined ? { wiggleFrequency: frequency } : {}),
        ...(tumble !== undefined ? { tumbleAmplitude: tumble } : {}),
        ...(tumbleSpeed !== undefined ? { tumbleSpeed } : {}),
        ...(tumbleFrequency !== undefined ? { tumbleFrequency } : {}),
        // As the classic's Procedural Animation: dynamics is the scene's amplitude (every representation); uncertainty
        // each component's layers. One mode at a time.
        wiggleAmplitude: nextMode === 'dynamics' ? nextAmplitude : 0,
      };
      // Uncertainty: the amplitude is the layers' largest value, so the same amplitude is the same wiggle.
      const layers =
        nextMode === 'uncertainty'
          ? !(now.mode === 'uncertainty' && Math.abs(uncertaintyScale(plugin) - nextAmplitude) < 1e-6)
          : wiggleCells(plugin).length > 0;
      const kept = areAnimationPropsEqual(options.animation, animation);
      if (kept && !layers) {
        return ApiResult.already(
          nextMode === 'off'
            ? animation.tumbleAmplitude > 0
              ? 'Nothing wiggles; the copies tumble so already.'
              : 'Nothing wiggles.'
            : `The structures already wiggle so (${nextMode}, amplitude ${nextAmplitude}).`,
        );
      }
      const components = structures.flatMap((s) => s.components);
      if (!kept) await manager.setOptions({ ...options, animation });
      if (nextMode !== 'uncertainty' && layers) await clearStructureWiggle(plugin, components);
      if (nextMode === 'uncertainty' && layers) {
        await setStructureWiggleFromUncertainty(plugin, components, nextAmplitude);
        if (!wiggleCells(plugin).length) {
          return ApiResult.fail(
            'blocked',
            'The structures have no uncertainty to wiggle by: their B-factors do not vary.',
            'Use mode dynamics.',
          );
        }
      }
      const tumbling = animation.tumbleAmplitude > 0 ? `; the copies tumble (${animation.tumbleAmplitude} Å)` : '';
      if (nextMode === 'off') {
        return ApiResult.ok(
          now.mode !== 'off'
            ? `Stopped the wiggle${tumbling}.`
            : tumbling
              ? `The copies tumble (${animation.tumbleAmplitude} Å).`
              : 'Stopped the tumble.',
          { nodes: nodesOf(ctx) },
        );
      }
      const started = now.mode === nextMode ? 'Changed' : 'Started';
      return ApiResult.ok(`${started} the wiggle (${nextMode}, amplitude ${nextAmplitude})${tumbling}.`, {
        nodes: nodesOf(ctx),
      });
    },
  });
}

const nodesOf = (ctx: ApiContext) =>
  ctx.plugin.managers.structure.hierarchy.current.structures.map((s) => s.cell.transform.ref);

const same = (a: unknown, b: unknown) => JSON.stringify(a) === JSON.stringify(b);
