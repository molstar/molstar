/** @jsxImportSource @molstar/plugin-prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { map } from 'rxjs';
import { ActionFooter } from '../../components/action-footer/action-footer.js';
import { ParamsForm } from '../../components/params/params-form.js';
import { Show } from '../../dom/control-flow.js';
import { domView } from '../../dom/view.js';
import { ParamsFormModel } from '../../forms/params-form.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ShellModel } from '../../shell/shell.js';
import { SectionList } from '../../views/shell/section-list.js';

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
  // The card is a list of sections, like Export: playing an animation, procedural animation (the wiggle), and
  // whatever extensions add.
  shell.registerOverlays([
    {
      id: 'animation',
      title: 'Animation',
      render: domView(({ ctx }) => <SectionList shell={ctx.shell} slot="animation" focus={ctx.focus} />),
    },
  ]);
  const wiggle = observableState<string | undefined>(undefined);
  const wiggleForm = new ParamsFormModel({
    by: PD.Select('dynamics', [
      ['dynamics', 'Dynamics (B-factor)'],
      ['uncertainty', 'Model uncertainty'],
    ]),
  });
  shell.registerSections([
    {
      id: 'play',
      slot: 'animation',
      title: 'Play',
      icon: 'play',
      order: 0,
      render: domView(() => (
        <div style="display: grid; gap: 8px">
          <ParamsForm form={form} />
          {/* Play, or Stop while playing: the footer's action follows. */}
          <Show
            when={playing}
            fallback={() => (
              <ActionFooter
                icon="play"
                label="Play"
                onReset={() => form.reset()}
                onApply={() => {
                  playing.set(true);
                  post(`Playing: ${label()}`);
                }}
              />
            )}
          >
            {() => (
              <ActionFooter
                icon="square"
                label="Stop"
                onReset={() => form.reset()}
                onApply={() => {
                  playing.set(false);
                  post('Animation stopped');
                }}
              />
            )}
          </Show>
        </div>
      )),
    },
    {
      id: 'procedural',
      slot: 'animation',
      title: 'Procedural animation',
      icon: 'orbit',
      keywords: ['wiggle', 'dynamics', 'uncertainty'],
      order: 1,
      // Mol*'s wiggle: motion from B-factors (dynamics) or from model uncertainty; Clear removes it.
      render: domView(() => (
        <div style="display: grid; gap: 8px">
          <ParamsForm form={wiggleForm} />
          <ActionFooter
            icon="orbit"
            label="Apply"
            info={wiggle.asReadonly.pipe(map((w) => (w ? `Wiggling by ${w}` : '')))}
            secondary={{
              icon: 'x',
              label: 'Clear',
              onClick: () => {
                wiggle.set(undefined);
                post('Wiggle cleared');
              },
            }}
            onApply={() => {
              const by = (wiggleForm.values.value as any).by as string;
              wiggle.set(by);
              post(`Wiggle: ${by}`);
            }}
          />
        </div>
      )),
    },
  ]);
  return { playing: playing.asReadonly };
}
