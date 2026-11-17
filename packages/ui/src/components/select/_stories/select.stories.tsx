/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../../_stories/support.js';
import { state } from '../../../state.js';
import { Select, type SelectOption } from '../select.js';

const meta: Meta = { title: 'Components/Select' };
export default meta;

// Room above: a list opened over its field centers the current option on it, and moves only at the window's edge.
const panel = 'width: 320px; margin: 128px var(--msui-space-4) var(--msui-space-4);';

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

const grouped: SelectOption[] = [
  ['uniform', 'Uniform'],
  ['element-symbol', 'Element Symbol', 'Atom Property'],
  ['occupancy', 'Occupancy', 'Atom Property'],
  ['chain-id', 'Chain Id', 'Chain Property'],
  ['entity-id', 'Entity Id', 'Chain Property'],
  ['secondary-structure', 'Secondary Structure', 'Residue Property'],
];

/**
 * In a form (`size="sm"`): each select a control in its row; opened (the default placement), the options float over it,
 * the current option on the field's text, and the rows after it stay where they are. A long list scrolls the current
 * option to its middle.
 */
export const InAForm: StoryObj = {
  render: () =>
    story(() => {
      const repr = state('cartoon');
      const color = state('chain-id');
      return (
        <div class="msui-surface" style={panel}>
          <div class="msui-form">
            <span class="msui-label">Type</span>
            <Select size="sm" aria-label="Type" options={representations} value={repr} onChange={(v) => repr.set(v)} />
            <span class="msui-label">Color</span>
            <Select
              size="sm"
              aria-label="Color theme"
              options={colorThemes}
              value={color}
              onChange={(v) => color.set(v)}
            />
            <span class="msui-label">Value</span>
            <span class="msui-value">
              {repr} · {color}
            </span>
          </div>
        </div>
      );
    }),
};

/** Grouped options: each group a fold (a chevron in its slot), the current option's open, its options a level in. */
export const Groups: StoryObj = {
  render: () =>
    story(() => {
      const value = state('entity-id');
      return (
        <div class="msui-surface" style={panel}>
          <div class="msui-form">
            <span class="msui-label">Color</span>
            <Select size="sm" aria-label="Color theme" options={grouped} value={value} onChange={(v) => value.set(v)} />
          </div>
        </div>
      );
    }),
};

/**
 * Inline (`over={false}`), the other placement: the options open under the field, in the flow, pushing the rows after
 * it down (a list meant to stay in a panel's flow).
 */
export const Inline: StoryObj = {
  render: () =>
    story(() => {
      const value = state('chain-id');
      return (
        <div class="msui-surface" style={panel}>
          <div class="msui-form">
            <span class="msui-label">Color</span>
            <Select
              size="sm"
              over={false}
              aria-label="Color theme"
              options={colorThemes}
              value={value}
              onChange={(v) => value.set(v)}
            />
            <span class="msui-label">Below</span>
            <span class="msui-value">This row moves down while the options are open.</span>
            <span class="msui-label">Value</span>
            <span class="msui-value">{value}</span>
          </div>
        </div>
      );
    }),
};

/** The two sizes, as Button's: `md` (default) a row's height, on its own; `sm` a control, in a form row. */
export const Sizes: StoryObj = {
  render: () =>
    story(() => {
      const md = state('cartoon');
      const sm = state('spacefill');
      return (
        <div class="msui-surface" style={panel}>
          <div class="msui-form">
            <span class="msui-label">md (24px)</span>
            <Select aria-label="Medium" options={representations} value={md} onChange={(v) => md.set(v)} />
            <span class="msui-label">sm (20px)</span>
            <Select size="sm" aria-label="Small" options={representations} value={sm} onChange={(v) => sm.set(v)} />
          </div>
        </div>
      );
    }),
};
