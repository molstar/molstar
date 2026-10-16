/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { MolstarWordmark } from '../../../components/icon/wordmark.js';
import { GlassWordmark } from '../glass-wordmark.js';
import { story } from '../../../testing/story-root.js';
import { PrismFieldDefaults, type PrismFieldOptions, PrismSplash } from '../prism-splash.js';

type Args = PrismFieldOptions;

const meta: Meta<Args> = {
  title: 'Components/Splash',
  args: { ...PrismFieldDefaults },
  argTypes: {
    rate: { control: { type: 'range', min: 1, max: 60, step: 1 }, description: 'Specks emitted per second' },
    speed: { control: { type: 'range', min: 4, max: 120, step: 1 }, description: 'Their speed (px/s)' },
    beam: { control: { type: 'range', min: 40, max: 600, step: 10 }, description: 'The stream before the prism (px)' },
    width: {
      control: { type: 'range', min: 0, max: 60, step: 1 },
      description: "The stream's half height at its start",
    },
    spread: { control: { type: 'range', min: 2, max: 60, step: 1 }, description: "A ray's half angle (degrees)" },
    reach: { control: { type: 'range', min: 60, max: 800, step: 10 }, description: 'How far the rays go (px)' },
    alpha: { control: { type: 'range', min: 0.05, max: 1, step: 0.01 }, description: 'The most a speck shows' },
    size: { control: { type: 'range', min: 0.3, max: 4, step: 0.1 }, description: "A speck's radius (px)" },
    variety: { control: { type: 'range', min: 0, max: 1, step: 0.05 }, description: 'How much their sizes vary' },
    bob: { control: { type: 'range', min: 0, max: 1, step: 0.05 }, description: 'Bubbles swell and shrink' },
    twinkle: { control: { type: 'range', min: 0, max: 1, step: 0.05 }, description: 'How much they twinkle' },
    swirl: { control: { type: 'range', min: 0, max: 2, step: 0.05 }, description: 'The pointer stirs an eddy' },
    colors: { control: 'boolean', description: "The rays in the logo's colors" },
    intro: {
      control: { type: 'range', min: 0, max: 6, step: 0.1 },
      description: 'Seconds the stream rushes in at the start',
    },
    pulse: { control: { type: 'range', min: 0, max: 1, step: 0.05 }, description: 'The light in waves' },
    curl: { control: { type: 'range', min: 0, max: 2, step: 0.05 }, description: 'The specks curl as they go' },
    breathe: { control: { type: 'range', min: 0, max: 1, step: 0.05 }, description: 'The rays breathe' },
    interactive: { control: 'boolean', description: 'It answers the pointer (mouse, pen, touch)' },
    radius: { control: { type: 'range', min: 20, max: 300, step: 5 }, description: "The pointer's reach (px)" },
    push: { control: { type: 'range', min: 0, max: 3000, step: 20 }, description: 'How hard it pushes specks aside' },
    wake: { control: { type: 'range', min: 0, max: 2, step: 0.05 }, description: 'How much of its stroke it lends' },
    follow: { control: { type: 'range', min: 0, max: 1, step: 0.05 }, description: 'More light toward the pointer' },
    burst: { control: { type: 'range', min: 0, max: 200, step: 2 }, description: 'Specks a press releases' },
    seed: { control: { type: 'number' } },
  },
};
export default meta;

/** The empty canvas as the shell shows it: the wordmark in the middle, the field under it. */
function Canvas(props: { theme: 'light' | 'dark'; options: Partial<PrismFieldOptions> }) {
  const markBox = (
    <div style="position: absolute; left: round(down, calc(50% - 85.5px), 1px); top: round(down, calc(50% - 33px), 1px); color: var(--msp-color-text-subtle)">
      <MolstarWordmark scale={3} flash="reveal" prism={false} />
    </div>
  ) as HTMLElement;
  const root = (
    <div
      data-msp-theme={props.theme}
      style={`position: relative; height: 520px; overflow: hidden; border-radius: var(--msp-radius-lg); background: ${props.theme === 'dark' ? '#000' : '#fff'}`}
    >
      <PrismSplash options={props.options} mark={() => markBox} pointer={() => root} />
      {markBox}
    </div>
  ) as HTMLElement;
  return root;
}

/**
 * The empty canvas's splash: light through the Mol* star, leaving it along its arms in the logo's colors. Move over it,
 * swipe through the rays, click or tap for a burst; tune it with the controls. It is meant to be barely there.
 */
export const Splash: StoryObj<Args> = {
  render: (args) =>
    story(() => (
      <div style="display: grid; gap: var(--msp-space-5); padding: var(--msp-space-5)">
        <Canvas theme="light" options={args} />
        <Canvas theme="dark" options={args} />
      </div>
    )),
};

/**
 * The wordmark as glass without the light: Home's heading. It assembles as it shows; move over it, press it. It draws
 * nothing while it is whole and nothing moves.
 */
export const Glass: StoryObj<Args> = {
  render: () =>
    story(() => (
      <div style="display: grid; gap: var(--msp-space-5); padding: 48px">
        {(['light', 'dark'] as const).map((theme) => (
          <div
            data-msp-theme={theme}
            style="display: flex; gap: 48px; align-items: end; padding: 32px; border-radius: var(--msp-radius-lg); background: var(--msp-color-panel)"
          >
            <GlassWordmark />
            <GlassWordmark scale={2} />
          </div>
        ))}
      </div>
    )),
};
