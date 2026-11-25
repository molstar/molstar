/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Heading, Note, story } from '../../../_stories/support.js';
import { type Child, onCleanup } from '../../../dom/index.js';
import { installTooltips } from '../../tooltip/tooltip.js';
import type { PlaybackLabels } from '../labels.js';
import { MuteButton, PlaybackDeck, StepCard as StepCardView, StepperRow, TimelineRow } from '../playback.js';
import { type Step, Stepper as StepperRef } from '../stepper-model.js';
import { DefaultSpeeds, Durationless as DurationlessRef, Timeline as TimelineRef } from '../timeline-model.js';

// The playback controls over the reference models (playback §3): stepped, timed (by time, by count) and durationless,
// a step's card, a deck of all three. The models tick while playing, and stop when the story goes.

const meta: Meta = { title: 'Components/Playback' };
export default meta;

const Tour: readonly Step[] = [
  {
    title: 'Overview',
    description:
      '# Overview\n\nThe whole assembly, colored by chain. Use **next** to step through, or [jump to the site](#3).',
  },
  {
    title: 'The binding pocket',
    description: 'A deep cleft between the two domains, lined with *hydrophobic* residues.',
  },
  { title: 'The ligand', description: '## The ligand\n\nHeld by three hydrogen bonds and a salt bridge.' },
  { title: 'Conservation', description: 'Residues colored from variable (light) to conserved (dark).' },
  { title: 'Summary', description: 'Back to the whole assembly.' },
];

const many: readonly Step[] = Array.from({ length: 32 }, (_, i) => ({
  title: `Scene ${i + 1}: ${['Overview', 'Pocket', 'Ligand', 'Interface', 'Surface', 'Cavity'][i % 6]}`,
}));

/** A stage-like backdrop for the floating controls, tooltips on. */
function Backdrop(props: { children: Child; width?: number }) {
  return (
    <div
      ref={(el: HTMLElement) => onCleanup(installTooltips(el))}
      style="display: grid; gap: var(--msui-space-4); align-content: start; padding: var(--msui-space-6); min-height: 100vh; background: var(--msui-surface-sunken);"
    >
      <div style={`display: grid; gap: var(--msui-space-4); width: ${props.width ?? 360}px;`}>{props.children}</div>
    </div>
  );
}

/** A model the story owns: stopped when the story goes. */
function owned<T extends { dispose(): void }>(model: T): T {
  onCleanup(() => model.dispose());
  return model;
}

/** An app's words for what one frame is (a structure's models). */
const ModelLabels: Partial<PlaybackLabels> = {
  frames: 'Models',
  frame: 'Model',
  previousFrame: 'Previous model',
  nextFrame: 'Next model',
  countedValue: (frame, count) => `model ${frame} of ${count}`,
};

/** Stepped: dots up to seven; past that "i / N" opens a filterable list. A mute at the end. */
export const Stepper: StoryObj = {
  render: () =>
    story(() => {
      const muted = new State(false);
      return (
        <Backdrop>
          <Heading>Five steps</Heading>
          <PlaybackDeck>
            <StepperRow
              model={owned(new StepperRef({ steps: Tour, intervalMs: 2000 }))}
              actions={() => <MuteButton muted={muted} onClick={() => muted.set(!muted.value)} />}
            />
          </PlaybackDeck>
          <Heading>Thirty-two scenes</Heading>
          <PlaybackDeck>
            <StepperRow model={owned(new StepperRef({ steps: many, index: 11 }))} label="Scenes" />
          </PlaybackDeck>
        </Backdrop>
      );
    }),
};

/** Timed (a length in seconds: the scrubber is a time) and counted (frames: the scrubber is an index). */
export const Timeline: StoryObj = {
  render: () =>
    story(() => (
      <Backdrop>
        <Heading>Timed</Heading>
        <PlaybackDeck>
          <TimelineRow model={owned(new TimelineRef({ duration: 2.4 }))} label="Transition" />
        </PlaybackDeck>
        <Heading>By count, with a speed</Heading>
        <PlaybackDeck>
          <TimelineRow
            model={owned(new TimelineRef({ count: 50, speeds: DefaultSpeeds }))}
            label="1ABC"
            labels={ModelLabels}
          />
        </PlaybackDeck>
      </Backdrop>
    )),
};

/**
 * Durationless: no length, so no scrubber; play, what it ran (time, steps), one step while paused, restart, a speed.
 * Running, paused, steps only, a spin with only play and its time.
 */
export const Durationless: StoryObj = {
  render: () =>
    story(() => {
      const running = owned(new DurationlessRef({ speeds: DefaultSpeeds }));
      running.play();
      const paused = owned(new DurationlessRef({ speeds: DefaultSpeeds }));
      for (let i = 0; i < 340; i++) paused.stepOnce?.();
      const steps = owned(new DurationlessRef({ elapsed: false, restart: false }));
      const spin = owned(new DurationlessRef({ steps: false, stepOnce: false, restart: false }));
      spin.play();
      return (
        <Backdrop>
          <Heading>Running</Heading>
          <PlaybackDeck>
            <TimelineRow model={running} label="Simulation" icon="orbit" />
          </PlaybackDeck>
          <Heading>Paused at a step</Heading>
          <PlaybackDeck>
            <TimelineRow model={paused} label="Simulation" icon="orbit" />
          </PlaybackDeck>
          <Heading>Steps only</Heading>
          <PlaybackDeck>
            <TimelineRow model={steps} label="Stream" icon="refresh-cw" />
          </PlaybackDeck>
          <Heading>A spin: play and its time</Heading>
          <PlaybackDeck>
            <TimelineRow model={spin} label="Spin" icon="rotate-3d" />
          </PlaybackDeck>
        </Backdrop>
      );
    }),
};

/** A step's card: its title, its position, its description (markdown; a heading that repeats the title is dropped). */
export const StepCard: StoryObj = {
  render: () =>
    story(() => (
      <Backdrop>
        <StepCardView model={owned(new StepperRef({ steps: Tour }))} label="Tour" />
        <Note>One step: the title and the description, nothing to step through.</Note>
        <StepCardView model={owned(new StepperRef({ steps: Tour.slice(0, 1) }))} label="Tour" />
      </Backdrop>
    )),
};

/** The deck with all three kinds, a row each. */
export const Deck: StoryObj = {
  render: () =>
    story(() => (
      <Backdrop width={400}>
        <PlaybackDeck>
          <StepperRow model={owned(new StepperRef({ steps: Tour }))} />
          <TimelineRow model={owned(new TimelineRef({ duration: 1.6 }))} label="Transition" />
          <TimelineRow model={owned(new TimelineRef({ count: 20 }))} label="1ABC" labels={ModelLabels} />
          <TimelineRow model={owned(new DurationlessRef({ steps: false, speeds: DefaultSpeeds }))} label="Spin" />
        </PlaybackDeck>
      </Backdrop>
    )),
};
