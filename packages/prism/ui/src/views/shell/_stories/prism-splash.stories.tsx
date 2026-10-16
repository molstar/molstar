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
    spread: { control: { type: 'range', min: 2, max: 60, step: 1 }, description: "The fan's half angle (degrees)" },
    reach: { control: { type: 'range', min: 60, max: 800, step: 10 }, description: 'How far the fan goes (px)' },
    alpha: { control: { type: 'range', min: 0.05, max: 1, step: 0.01 }, description: 'The most a speck shows' },
    size: { control: { type: 'range', min: 0.3, max: 3, step: 0.1 }, description: "A speck's radius (px)" },
    tail: { control: { type: 'range', min: 0, max: 30, step: 1 }, description: 'A streak behind a colored speck (px)' },
    colors: { control: 'boolean', description: "The fan in the logo's colors" },
    seed: { control: { type: 'number' } },
  },
};
export default meta;

/** The empty canvas as the shell shows it: the wordmark in the middle, the field under it. */
function Canvas(props: { theme: 'light' | 'dark'; options: Partial<PrismFieldOptions> }) {
  const markBox = (
    <div style="position: absolute; left: round(down, calc(50% - 85.5px), 1px); top: round(down, calc(50% - 33px), 1px); color: var(--msp-color-text-subtle)">
      <MolstarWordmark scale={3} flash="reveal" />
    </div>
  ) as HTMLElement;
  return (
    <div
      data-msp-theme={props.theme}
      style={`position: relative; height: 520px; overflow: hidden; border-radius: var(--msp-radius-lg); background: ${props.theme === 'dark' ? '#000' : '#fff'}`}
    >
      <PrismSplash options={props.options} mark={() => markBox} />
      {markBox}
    </div>
  );
}

/**
 * The empty canvas's splash: light through Prism's mark, splitting into the logo's colors. Tune it with the controls;
 * it is meant to be barely there.
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
