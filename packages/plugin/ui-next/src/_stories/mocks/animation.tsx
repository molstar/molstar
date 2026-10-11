/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { map } from 'rxjs';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { ParamsForm } from '../../components/params/params-form.js';
import { Show } from '../../dom/control-flow.js';
import { domView } from '../../dom/view.js';
import { ParamsFormModel } from '../../forms/params-form.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ShellModel } from '../../shell/shell.js';

const Animations = [
  ['trajectory', 'Animate trajectory'],
  ['spin', 'Camera spin'],
  ['rock', 'Camera rock'],
  ['snapshots', 'State snapshots'],
  ['unwind', 'Assembly unwind'],
] as [any, string][];

// The plugin's animations (`plugin.managers.animation`), each with its own params, as Mol*'s animation panel has them.
const AnimationParams = {
  animation: PD.MappedStatic(
    'trajectory',
    {
      trajectory: PD.Group({
        mode: PD.Select('loop', [
          ['loop', 'Loop'],
          ['palindrome', 'Palindrome'],
          ['once', 'Once'],
        ]),
        duration: PD.Numeric(5, { min: 1, max: 60, step: 0.5 }, { label: 'Duration (s)' }),
      }),
      spin: PD.Group({
        speed: PD.Numeric(1, { min: 0.1, max: 5, step: 0.1 }, { label: 'Speed' }),
        direction: PD.Select('cw', [
          ['cw', 'Clockwise'],
          ['ccw', 'Counter-clockwise'],
        ]),
      }),
      rock: PD.Group({
        angle: PD.Numeric(15, { min: 1, max: 90, step: 1 }, { label: 'Angle (°)' }),
        speed: PD.Numeric(0.3, { min: 0.05, max: 2, step: 0.05 }, { label: 'Speed' }),
      }),
      snapshots: PD.Group({ duration: PD.Numeric(2, { min: 0.5, max: 10, step: 0.5 }, { label: 'Per snapshot (s)' }) }),
      unwind: PD.Group({ duration: PD.Numeric(3, { min: 0.5, max: 10, step: 0.5 }, { label: 'Duration (s)' }) }),
    },
    { label: 'Animation', isEssential: true, options: [...Animations] },
  ),
};

/** The Animation card (rail ▶): choose an animation, its options, play or stop. */
export function registerMockAnimation(shell: ShellModel, post: (text: string) => void) {
  const form = new ParamsFormModel(AnimationParams);
  const playing = observableState(false);
  const label = () => {
    const name = (form.values.value as any).animation.name as string;
    return Animations.find(([v]) => v === name)?.[1] ?? name;
  };
  shell.registerOverlays([
    {
      id: 'animation',
      title: 'Animation',
      render: domView(() => (
        <div style="display: grid; gap: 8px; padding: 4px 12px 12px">
          <ParamsForm form={form} />
          <div style="display: flex; justify-content: flex-end">
            <Button
              variant="primary"
              size="sm"
              icon={
                <Show when={playing} fallback={() => <Icon name="play" />}>
                  {() => <Icon name="square" />}
                </Show>
              }
              onClick={() => {
                playing.set(!playing.value);
                post(playing.value ? `Playing: ${label()}` : 'Animation stopped');
              }}
            >
              {playing.asReadonly.pipe(map((p) => (p ? 'Stop' : 'Play')))}
            </Button>
          </div>
        </div>
      )),
    },
  ]);
  return { playing: playing.asReadonly };
}
