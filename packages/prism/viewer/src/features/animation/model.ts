/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type ReadonlyState, mapState } from '@molstar/prism-ui/core';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { combineLatest, distinctUntilChanged, map } from 'rxjs';
import type { ViewerApi } from '../../api/area.js';
import type { AnimationApi } from './operations.js';
import type { WiggleState } from './source.js';

// The Animation card's view model (state.md §4, a feature's own contract; it could move to `src/model/` with the others
// if another UI wants it): what the card's sections read, over the animation area's state and `params`, and what they
// do, through the person's API, never the plugin. DOM-free.

export interface AnimationModel {
  /** The choice of animation and its params, as one mapped param (`animation: { name, params }`); none registered: none. */
  readonly form: ParamsFormModel<any> | undefined;
  /** Why the chosen animation cannot play now ("No trajectory to animate"); undefined: it can. */
  readonly blocked: ReadonlyState<string | undefined>;
  /** Any animation plays (the Play button is Stop). */
  readonly isPlaying: ReadonlyState<boolean>;
  play(): void;
  stop(): void;
  /** The procedural animation: mode (dynamics, uncertainty) and its strength. */
  readonly wiggleForm: ParamsFormModel<any>;
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
};

export function createAnimationModel(api: ViewerApi, area: AnimationApi): AnimationModel {
  const animations = area.state.value.animations;
  // The first that can play now is the choice to start with.
  const first = animations.find((a) => a.available === true) ?? animations[0];
  const form = first
    ? new ParamsFormModel({
        animation: PD.MappedStatic(
          first.name,
          Object.fromEntries(
            animations.map((a) => [a.name, PD.Group((area.paramsOf(a.name) ?? {}) as PD.Params, { isFlat: true })]),
          ),
          {
            label: 'Animation',
            isEssential: true,
            options: animations.map((a) => [a.name, a.title] as [string, string]),
          },
        ),
      })
    : undefined;
  const chosen = () => (form?.values.value as { animation?: { name: string } } | undefined)?.animation?.name;

  const reason = () => {
    const name = chosen();
    const info = name ? area.state.value.animations.find((a) => a.name === name) : undefined;
    return info && info.available !== true ? info.available : undefined;
  };
  // Both what is available and the choice change it; as `mapState`, nothing is held until someone subscribes.
  const blocked = (
    form
      ? combineLatest([area.state, form.values]).pipe(map(reason), distinctUntilChanged())
      : area.state.pipe(map(reason), distinctUntilChanged())
  ) as ReadonlyState<string | undefined>;
  Object.defineProperty(blocked, 'value', { get: reason, enumerable: true });

  const wiggleForm = new ParamsFormModel(WiggleFormParams);
  return {
    form,
    blocked,
    isPlaying: mapState(area.state, (s) => s.playing !== undefined, Object.is),
    play() {
      const v = form?.values.value as { animation?: { name: string; params: Record<string, unknown> } } | undefined;
      if (!v?.animation) return;
      void api.animation?.play({ name: v.animation.name, params: v.animation.params });
    },
    stop() {
      void api.animation?.stop();
    },
    wiggleForm,
    wiggle: mapState(
      area.state,
      (s) => s.wiggle,
      (a, b) => a.mode === b.mode && a.amplitude === b.amplitude,
    ),
    applyWiggle() {
      const v = wiggleForm.values.value;
      void api.animation?.wiggle({ mode: v.mode, amplitude: v.amplitude });
    },
    clearWiggle() {
      void api.animation?.wiggle({ mode: 'off' });
    },
  };
}
