/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { domView, Icon, ParamsForm, Show } from '@molstar/prism-ui';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { onActive, onCleanup, requireScope } from '@molstar/prism-ui/dom/scope';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { ShellSection } from '@molstar/prism-ui/shell/shell';
import type { SceneNode, SceneTrajectory } from '../../model/scene.js';
import { classes } from './_generated/trajectory-row.module.css.js';
import { SceneRow } from './scene-row.js';
import { SceneRowSlot } from './scene-row-slot.js';

// A structure's frames in the Structure panel (an NMR ensemble's models, a simulation's frames): a row of its own, the
// first under the structure, with the frame, a scrubber, ‹ ›, and play; picked, the inspector has the same at full
// width with how it plays (mode, seconds per pass). While the row is on screen the canvas' bar steps back to a compact
// chip (shown while frames play); otherwise the canvas has the full bar (ui/shell.tsx).

const Capital = { models: 'Models', frames: 'Frames' } as const;

/**
 * The trajectory row as a section of the tree's slot (`SceneRowSlot`), so the tree hosts it as it hosts an extension's
 * row (picked, its controls in the inspector); not registered in the shell: the tree puts it first under the structure.
 */
export const TrajectoryRowSection: ShellSection = {
  id: 'scene.trajectory',
  slot: SceneRowSlot,
  title: 'Models frames trajectory',
  when: (c) => !!(c.subject as SceneNode | undefined)?.trajectory,
  render: domView(({ ctx }) => {
    const trajectory = (ctx.subject as SceneNode).trajectory;
    return trajectory ? <TrajectoryRow trajectory={trajectory} /> : null;
  }),
};

function TrajectoryRow(props: { trajectory: SceneTrajectory }) {
  const t = props.trajectory;
  // On screen while its view is active (the Structure tab shown, the structure open, the row in the tree's window): the
  // canvas' bar steps back then (ui/shell.tsx). Hidden, collapsed, searched away, scrolled out or expanded over: not.
  const scope = requireScope('TrajectoryRow');
  let shown = false;
  const show = (on: boolean) => {
    if (on === shown) return;
    shown = on;
    t.markShown(on);
  };
  show(scope.active);
  onActive(show);
  onCleanup(() => show(false));
  return (
    <SceneRow
      label={Capital[t.unit]}
      description={derive([t.index, t.count], (i, c) => `${i + 1} / ${c}`)}
      inline={() => <Scrubber trajectory={t} compact />}
      body={() => <TrajectoryControls trajectory={t} />}
    />
  );
}

/** The scrubber, ‹ › and play: on the row (`compact`, out of the tree's keyboard) and in the inspector. */
function Scrubber(props: { trajectory: SceneTrajectory; compact?: boolean }) {
  const t = props.trajectory;
  const one = t.unit === 'models' ? 'model' : 'frame';
  const tab = props.compact ? -1 : 0;
  return (
    <span class={props.compact ? `${classes.scrubber} ${classes.compact}` : classes.scrubber}>
      <input
        class={classes.range}
        type="range"
        min="0"
        max={derive(t.count, (c) => String(Math.max(0, c - 1)))}
        value={derive(t.index, String)}
        tabIndex={tab}
        style={{
          '--msp-progress': derive([t.index, t.count], (i, c) => `${c > 1 ? (100 * i) / (c - 1) : 0}%`),
        }}
        aria-label={`${Capital[t.unit].slice(0, -1)} shown`}
        // Said as people count them ("Model 3 of 10"), not the 0-based value.
        aria-valuetext={derive([t.index, t.count], (i, c) => `${Capital[t.unit].slice(0, -1)} ${i + 1} of ${c}`)}
        onInput={(e) => t.seek(Number((e.currentTarget as HTMLInputElement).value))}
      />
      <button type="button" class={classes.button} tabIndex={tab} title={`Previous ${one}`} onClick={() => t.step(-1)}>
        <Icon name="chevron-left" size={14} />
      </button>
      <button
        type="button"
        class={classes.button}
        tabIndex={tab}
        aria-pressed={t.playing}
        title={derive(t.playing, (p) => (p ? 'Pause' : `Play the ${t.unit}`))}
        onClick={() => t.togglePlay()}
      >
        <Show when={t.playing} fallback={() => <Icon name="play" size={14} />}>
          {() => <Icon name="pause" size={14} />}
        </Show>
      </button>
      <button type="button" class={classes.button} tabIndex={tab} title={`Next ${one}`} onClick={() => t.step(1)}>
        <Icon name="chevron-right" size={14} />
      </button>
    </span>
  );
}

const PlayParams = {
  mode: PD.Select<'loop' | 'palindrome' | 'once'>(
    'loop',
    [
      ['loop', 'Loop'],
      ['palindrome', 'Back and forth'],
      ['once', 'Once'],
    ],
    { label: 'Play' },
  ),
  seconds: PD.Numeric(5, { min: 1, max: 120, step: 0.5 }, { label: 'Seconds per pass' }),
};

/** The inspector's controls (its header says the frame shown): the scrubber at full width, and how it plays. */
function TrajectoryControls(props: { trajectory: SceneTrajectory }) {
  const t = props.trajectory;
  const form = new ParamsFormModel(
    PlayParams,
    { mode: t.mode.value, seconds: t.seconds.value },
    {
      onChange: (v) => {
        t.mode.set(v.mode);
        t.seconds.set(v.seconds);
      },
    },
  );
  return (
    <div class={classes.controls}>
      <Scrubber trajectory={t} />
      <ParamsForm form={form} />
    </div>
  );
}
