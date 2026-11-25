/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { type Child, onCleanup } from '../../../dom/index.js';
import { Heading, Note, story } from '../../../_stories/support.js';
import { state } from '../../../state.js';
import { installTooltips } from '../../tooltip/tooltip.js';
import { Segmented, type SegmentedOption } from '../segmented.js';

const meta: Meta = { title: 'Components/Segmented' };
export default meta;

/** A surface with tooltips: what each story is in. */
function Panel(props: { width?: number; children: Child }) {
  return (
    <div
      ref={(el: HTMLElement) => onCleanup(installTooltips(el))}
      class="msui-surface"
      style={`display: grid; gap: var(--msui-space-4); width: ${props.width ?? 320}px; margin: var(--msui-space-4); padding: var(--msui-space-3); border-radius: var(--msui-radius-lg);`}
    >
      {props.children}
    </div>
  );
}

const types: SegmentedOption<string>[] = [
  { value: 'cartoon', label: 'Cartoon', icon: 'cartoon' },
  { value: 'ball-and-stick', label: 'Ball & stick', icon: 'ball-and-stick' },
  { value: 'spacefill', label: 'Spacefill', icon: 'spacefill' },
  { value: 'molecular-surface', label: 'Surface', icon: 'surface' },
  { value: 'gaussian-surface', label: 'Gaussian surface', icon: 'volume' },
  { value: 'putty', label: 'Putty', icon: 'worm' },
  { value: 'line', label: 'Line', icon: 'line-squiggle' },
];

/**
 * A representation's type: the types as their icons, the chosen one raised in the well; the heading names it. Hover a
 * segment for its name; Tab reaches the chosen one, the arrows, Home and End move and pick.
 */
export const RepresentationType: StoryObj = {
  render: () =>
    story(() => {
      const type = state('ball-and-stick');
      return (
        <Panel>
          <Segmented label="Type" options={types} value={type} iconOnly />
        </Panel>
      );
    }),
};

/** Text segments: the labels themselves, equal while they fit. With icons and text, the icon leads. */
export const Text: StoryObj = {
  render: () =>
    story(() => {
      const quality = state('auto');
      const view = state('ortho');
      return (
        <Panel>
          <Segmented
            label="Quality"
            value={quality}
            options={[
              { value: 'auto', label: 'Auto' },
              { value: 'low', label: 'Low' },
              { value: 'medium', label: 'Medium' },
              { value: 'high', label: 'High' },
            ]}
          />
          <Segmented
            label="Camera"
            value={view}
            options={[
              { value: 'ortho', label: 'Orthographic', icon: 'box' },
              { value: 'persp', label: 'Perspective', icon: 'axis-3d' },
            ]}
          />
        </Panel>
      );
    }),
};

/**
 * The sizes: `md` (segments a row high in a well a header high, 16px icons) on its own; `sm` (segments a control high in
 * a well a row high, 12px icons) in a form row, its label in the label column.
 */
export const Sizes: StoryObj = {
  render: () =>
    story(() => {
      const type = state('cartoon');
      const align = state('left');
      return (
        <Panel width={360}>
          <Heading>md</Heading>
          <Segmented aria-label="Type" options={types.slice(0, 4)} value={type} iconOnly />
          <Heading>sm, in a form row</Heading>
          <div class="msui-form">
            <span class="msui-label">Label position</span>
            <Segmented
              size="sm"
              aria-label="Label position"
              value={align}
              options={[
                { value: 'left', label: 'Left' },
                { value: 'center', label: 'Center' },
                { value: 'right', label: 'Right' },
              ]}
            />
            <span class="msui-label">Type</span>
            <Segmented size="sm" aria-label="Type" options={types.slice(0, 4)} value={type} iconOnly />
          </div>
        </Panel>
      );
    }),
};

/**
 * Types that do not apply to this structure: greyed, never picked, skipped by the arrows; one says why in its tooltip.
 */
export const Disabled: StoryObj = {
  render: () =>
    story(() => {
      const type = state('cartoon');
      const options: SegmentedOption<string>[] = types.map((o) =>
        o.value === 'gaussian-surface'
          ? { ...o, disabled: 'Gaussian surface: too many atoms' }
          : o.value === 'putty'
            ? { ...o, disabled: true }
            : o,
      );
      return (
        <Panel>
          <Segmented label="Type" options={options} value={type} iconOnly />
        </Panel>
      );
    }),
};

/**
 * In a narrow panel (a docked side panel): the segments shrink to squares and stay whole; the heading's value wraps
 * under its label when they do not fit on one line.
 */
export const NarrowPanel: StoryObj = {
  render: () =>
    story(() => {
      const type = state('gaussian-surface');
      return (
        <Panel width={200}>
          <Segmented label="Representation type" options={types} value={type} iconOnly />
          <Note>200px wide: seven 24px segments and the well's padding fit in its 176px.</Note>
        </Panel>
      );
    }),
};
