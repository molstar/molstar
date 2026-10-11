/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../testing/story-root.js';
import { observableState } from '../../reactive/observable-state.js';
import { Select, type SelectOption } from './select.js';

const meta: Meta = { title: 'Components/Select' };
export default meta;

const panel =
  'display: grid; gap: var(--msp-space-3); width: 300px; padding: var(--msp-space-4); background: var(--msp-color-panel); color: var(--msp-color-text); font-size: var(--msp-font-size-sm);';
const readout = 'color: var(--msp-color-text-muted);';

const representations: SelectOption[] = [
  ['cartoon', 'Cartoon'],
  ['ball-and-stick', 'Ball & Stick'],
  ['spacefill', 'Spacefill'],
  ['molecular-surface', 'Molecular Surface'],
];

const colorThemes: SelectOption[] = [
  'atom-id',
  'carbohydrate-symbol',
  'chain-id',
  'element-index',
  'element-symbol',
  'entity-id',
  'entity-source',
  'hydrophobicity',
  'illumination',
  'model-index',
  'molecule-type',
  'occupancy',
  'operator-name',
  'partial-charge',
  'polymer-id',
  'residue-name',
  'secondary-structure',
  'sequence-id',
  'shape-group',
  'uncertainty',
  'uniform',
].map((v) => [v, v.replace(/-/g, ' ').replace(/^./, (c) => c.toUpperCase())] as const);

/** The small and the medium select side by side, each with its value below. */
export const Sizes: StoryObj = {
  render: () =>
    story(() => {
      const sm = observableState('cartoon');
      const md = observableState('spacefill');
      return (
        <div style={panel}>
          <Select
            size="sm"
            aria-label="Small"
            options={representations}
            value={sm.asReadonly}
            onChange={(v) => sm.set(v)}
          />
          <span style={readout}>value: {sm.asReadonly}</span>
          <Select
            size="md"
            aria-label="Medium"
            options={representations}
            value={md.asReadonly}
            onChange={(v) => md.set(v)}
          />
          <span style={readout}>value: {md.asReadonly}</span>
        </div>
      );
    }),
};

/** A long option list: opened, it scrolls, with the current option marked; arrows, Home and End move. */
export const LongList: StoryObj = {
  render: () =>
    story(() => {
      const value = observableState('chain-id');
      return (
        <div style={panel}>
          <Select
            aria-label="Color theme"
            options={colorThemes}
            value={value.asReadonly}
            onChange={(v) => value.set(v)}
          />
          <span style={readout}>value: {value.asReadonly}</span>
        </div>
      );
    }),
};

/**
 * Floating (`over`): the options open over what is below instead of pushing it down, for a short settings row; flips
 * above when there is no room below.
 */
export const Over: StoryObj = {
  render: () =>
    story(() => {
      const value = observableState('chain-id');
      const row =
        'display: grid; grid-template-columns: 96px minmax(0, 1fr); align-items: center; gap: var(--msp-space-3);';
      return (
        <div style={panel}>
          <div style={row}>
            <span style={readout}>Color</span>
            <Select
              over
              size="sm"
              aria-label="Color theme"
              options={colorThemes}
              value={value.asReadonly}
              onChange={(v) => value.set(v)}
            />
          </div>
          <p style={`${readout} margin: 0`}>This line stays where it is while the options are open.</p>
          <span style={readout}>value: {value.asReadonly}</span>
        </div>
      );
    }),
};
