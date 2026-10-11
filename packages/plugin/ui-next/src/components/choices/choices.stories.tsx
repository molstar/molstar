/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../_stories/story-root.js';
import { observableState } from '../../reactive/observable-state.js';
import { Choices } from './choices.js';

const meta: Meta = { title: 'Components/Choices' };
export default meta;

const panel =
  'display: grid; gap: var(--msp-space-3); width: 300px; padding: var(--msp-space-4); background: var(--msp-color-panel); color: var(--msp-color-text); font-size: var(--msp-font-size-sm);';
const readout = 'color: var(--msp-color-text-muted);';

/** Quick styles with their icons: picking one marks it as the one in effect. */
export const QuickStyles: StoryObj = {
  render: () =>
    story(() => {
      const active = observableState<string | undefined>('cartoon');
      return (
        <div style={panel}>
          <Choices
            label="Quick styles"
            choices={[
              { id: 'default', label: 'Default', icon: 'molstar' },
              { id: 'cartoon', label: 'Cartoon', icon: 'cartoon' },
              { id: 'spacefill', label: 'Spacefill', icon: 'spacefill' },
              { id: 'surface', label: 'Surface', icon: 'surface' },
            ]}
            active={active.asReadonly}
            onPick={(id) => active.set(id)}
          />
          <span style={readout}>active: {active.asReadonly}</span>
        </div>
      );
    }),
};

/** Without a label or an active state: plain one-click actions, the pick echoed below. */
export const Plain: StoryObj = {
  render: () =>
    story(() => {
      const picked = observableState('');
      return (
        <div style={panel}>
          <Choices
            choices={[
              { id: 'dynamics', label: 'Dynamics' },
              { id: 'uncertainty', label: 'Uncertainty' },
              { id: 'clear', label: 'Clear' },
            ]}
            onPick={(id) => picked.set(id)}
          />
          <span style={readout}>picked: {picked.asReadonly}</span>
        </div>
      );
    }),
};
