/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type ReadonlyState, mapState, observableState } from '@molstar/prism-ui/core';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { combineLatest, distinctUntilChanged, map, of, switchMap } from 'rxjs';
import type { ViewerApi } from '../../api/area.js';
import type { AnimationApi } from './operations.js';
import { sameWiggle, type WiggleState, type WiggleType } from './source.js';

// The Animation card's view model (state.md §4, a feature's own contract; it could move to `src/model/` with the others
// if another UI wants it): what the card's sections read, over the animation area's state and `params`, and what they
// do, through the person's API, never the plugin. DOM-free.

export interface AnimationModel {
  /** The choice of animation and its params, as one mapped param (`animation: { name, params }`); none registered: none. */
  readonly form: ParamsFormModel<any> | undefined;
  /**
   * Builds `form` again when the animations or their params changed since (one registered later, a structure loaded
   * that Unwind can target), keeping the choice and the values whose param is the same; the Play section calls it when
   * it shows.
   */
  refresh(): void;
  /** Why the chosen animation cannot play now ("No trajectory to animate"); undefined: it can. */
  readonly blocked: ReadonlyState<string | undefined>;
  /** Any animation plays (the Play button is Stop). */
  readonly isPlaying: ReadonlyState<boolean>;
  play(): void;
  stop(): void;
  /** The procedural animation: mode (dynamics, uncertainty), strength, type, speed, frequency, and the tumble. */
  readonly wiggleForm: ParamsFormModel<any>;
  /** The form from the scene's wiggle (the Wiggle section, whenever the wiggle changes). */
  refreshWiggle(): void;
  readonly wiggle: ReadonlyState<WiggleState>;
  applyWiggle(): void;
  clearWiggle(): void;
}

const WiggleFormParams = {
  mode: PD.Select<'dynamics' | 'uncertainty'>(
    'dynamics',
    [
      ['dynamics', 'Dynamics (B-factor)'],
      ['uncertainty', 'Model uncertainty'],
    ],
    { description: 'Move atoms by noise, or by how uncertain their position is (B-factor or RMSF)' },
  ),
  amplitude: PD.Numeric(1, { min: 0.01, max: 5, step: 0.01 }, { label: 'Amplitude' }),
  type: PD.Select<WiggleType>(
    'position',
    [
      ['position', 'Position'],
      ['group', 'Group'],
    ],
    { label: 'Type', description: 'Position: nearby atoms move together. Group: each atom (group) on its own.' },
  ),
  speed: PD.Numeric(7, { min: 0, max: 10, step: 0.1 }, { label: 'Speed' }),
  frequency: PD.Numeric(
    0.2,
    { min: 0.01, max: 2, step: 0.01 },
    { label: 'Frequency', description: 'Lower: nearby atoms move more together (position type)' },
  ),
  tumble: PD.Group(
    {
      amplitude: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }, { label: 'Amplitude', description: 'In Å; 0: off' }),
      speed: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, { label: 'Speed' }),
      frequency: PD.Numeric(0.2, { min: 0, max: 2, step: 0.01 }, { label: 'Frequency' }),
    },
    { label: 'Tumble', description: 'Whole copies (instances) turning in place', isExpanded: false },
  ),
};

type WiggleFormValues = PD.Values<typeof WiggleFormParams>;

/** The form's values from the scene's wiggle (off: the mode and amplitude the form has, the rest the scene's). */
function wiggleValues(w: WiggleState, form: WiggleFormValues): WiggleFormValues {
  return {
    mode: w.mode === 'off' ? form.mode : w.mode,
    amplitude: w.mode === 'off' ? form.amplitude : w.amplitude,
    type: w.type,
    speed: w.speed,
    frequency: w.frequency,
    tumble: { amplitude: w.tumble, speed: w.tumbleSpeed, frequency: w.tumbleFrequency },
  };
}

const Axes = ['x', 'y', 'z'] as const;

/**
 * A vector param as three numbers labeled "Axis x", "Axis y", "Axis z" (the form's own vector field has no labels on
 * its numbers). The values are `{ x, y, z }` in the form and `[x, y, z]` for the animation: `joinVectors` converts.
 */
function splitVectors(params: PD.Params): PD.Params {
  const out: PD.Params = {};
  for (const [key, p] of Object.entries(params)) {
    if (p.type !== 'vec3') {
      out[key] = p;
      continue;
    }
    const name = p.label ?? key.replace(/^./, (c) => c.toUpperCase());
    const [x, y, z] = p.defaultValue;
    const range = { min: p.min, max: p.max, step: p.step ?? 0.1 };
    out[key] = PD.Group(
      Object.fromEntries(
        Axes.map((a, i) => [a, PD.Numeric([x, y, z][i], range, { label: `${name} ${a}`, description: p.description })]),
      ),
      { isFlat: true },
    );
  }
  return out;
}

function joinVectors(params: PD.Params, values: Record<string, unknown>): Record<string, unknown> {
  const out = { ...values };
  for (const [key, p] of Object.entries(params)) {
    const v = values[key] as { x: number; y: number; z: number } | undefined;
    if (p.type === 'vec3' && v && !Array.isArray(v)) out[key] = [v.x, v.y, v.z];
  }
  return out;
}

/** A param's definition as text, to compare (undefined when it cannot be written: then it never compares equal). */
function json(value: unknown): string | undefined {
  try {
    return JSON.stringify(value);
  } catch {
    return undefined;
  }
}

export function createAnimationModel(api: ViewerApi, area: AnimationApi): AnimationModel {
  const current = observableState<ParamsFormModel<any> | undefined>(undefined, { compare: Object.is });
  // What the form was built from: each animation's params (as the form has them) and their key, so a refresh with
  // nothing changed keeps the form, its folds and values.
  let built: { key: string | undefined; groups: Record<string, PD.Params> } | undefined;
  const refresh = () => {
    const animations = area.state.value.animations;
    const groups: Record<string, PD.Params> = {};
    for (const a of animations) groups[a.id] = splitVectors(area.paramsOf(a.id) ?? {});
    const key = json(animations.map((a) => [a.id, a.title, groups[a.id]]));
    if (key !== undefined && built?.key === key) return;
    const previous = current.value?.values.value.animation as
      | { name: string; params: Record<string, unknown> }
      | undefined;
    const before = built?.groups;
    built = { key, groups };
    // The first that can play now is the choice to start with.
    const first = animations.find((a) => a.available === true) ?? animations[0];
    if (!first) return current.set(undefined);
    const animation = PD.MappedStatic(
      first.id,
      Object.fromEntries(animations.map((a) => [a.id, PD.Group(groups[a.id], { isFlat: true })])),
      {
        label: 'Animation',
        isEssential: true,
        options: animations.map((a) => [a.id, a.title] as [string, string]),
      },
    );
    let values: { animation: { name: string; params: Record<string, unknown> } } | undefined;
    if (previous && groups[previous.name] && before?.[previous.name]) {
      // A value goes on where its param is as it was (a Select whose options changed starts again from its default).
      const was = before[previous.name];
      const params = { ...PD.getDefaultValues(groups[previous.name]) } as Record<string, unknown>;
      for (const [k, p] of Object.entries(groups[previous.name])) {
        const same = was[k] && json(was[k]) !== undefined && json(was[k]) === json(p);
        if (same && k in previous.params) params[k] = previous.params[k];
      }
      values = { animation: { name: previous.name, params } };
    }
    current.set(new ParamsFormModel({ animation }, values as never));
  };
  refresh();
  const chosen = () => (current.value?.values.value as { animation?: { name: string } } | undefined)?.animation?.name;

  const reason = () => {
    const name = chosen();
    const info = name ? area.state.value.animations.find((a) => a.id === name) : undefined;
    return info && info.available !== true ? info.available : undefined;
  };
  // Both what is available and the choice change it; as `mapState`, nothing is held until someone subscribes.
  const blocked = combineLatest([
    area.state,
    current.asReadonly.pipe(switchMap((form) => (form ? form.values : of(undefined)))),
  ]).pipe(map(reason), distinctUntilChanged()) as ReadonlyState<string | undefined>;
  Object.defineProperty(blocked, 'value', { get: reason, enumerable: true });

  const wiggleForm = new ParamsFormModel(WiggleFormParams);
  const refreshWiggle = () => wiggleForm.load(wiggleValues(area.state.value.wiggle, wiggleForm.values.value));
  refreshWiggle();
  return {
    get form() {
      return current.value;
    },
    refresh,
    blocked,
    isPlaying: mapState(area.state, (s) => s.playing !== undefined, Object.is),
    play() {
      const v = current.value?.values.value as
        | { animation?: { name: string; params: Record<string, unknown> } }
        | undefined;
      if (!v?.animation) return;
      void api.animation?.play({
        name: v.animation.name,
        params: joinVectors(area.paramsOf(v.animation.name) ?? {}, v.animation.params),
      });
    },
    stop() {
      void api.animation?.stop();
    },
    wiggleForm,
    refreshWiggle,
    wiggle: mapState(area.state, (s) => s.wiggle, sameWiggle),
    applyWiggle() {
      const v = wiggleForm.values.value;
      void api.animation?.wiggle({
        mode: v.mode,
        amplitude: v.amplitude,
        type: v.type,
        speed: v.speed,
        frequency: v.frequency,
        tumble: v.tumble.amplitude,
        tumbleSpeed: v.tumble.speed,
        tumbleFrequency: v.tumble.frequency,
      });
    },
    clearWiggle() {
      void api.animation?.wiggle({ mode: 'off', tumble: 0 });
    },
  };
}
