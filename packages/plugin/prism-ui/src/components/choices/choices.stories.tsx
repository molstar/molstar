/** @jsxImportSource @molstar/plugin-prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../_stories/story-root.js';
import { onCleanup } from '../../dom/scope.js';
import { observableState } from '../../reactive/observable-state.js';
import { installTooltips } from '../tooltip/tooltip.js';
import { Choices } from './choices.js';

const meta: Meta = { title: 'Components/Choices' };
export default meta;

const panel =
  'display: grid; gap: var(--msp-space-3); width: 300px; padding: var(--msp-space-4); background: var(--msp-color-panel); color: var(--msp-color-text); font-size: var(--msp-font-size-sm);';
const readout = 'color: var(--msp-color-text-muted);';

/** Quick styles (as the Structure tab has them): picking one marks it as the one in effect. */
export const QuickStyles: StoryObj = {
  render: () =>
    story(() => {
      const active = observableState<string | undefined>('cartoon');
      return (
        <div style={panel}>
          <Choices
            label="Quick styles"
            choices={[
              { id: 'default', label: 'Default' },
              { id: 'cartoon', label: 'Cartoon' },
              { id: 'spacefill', label: 'Spacefill' },
              { id: 'surface', label: 'Surface' },
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

/** Narrow: the choices wrap into rows of equal widths; a label that still does not fit gives way (full on hover). */
export const Narrow: StoryObj = {
  render: () =>
    story(() => (
      <div ref={(el: HTMLElement) => onCleanup(installTooltips(el))} style={`${panel} width: 180px;`}>
        <Choices
          label="Representation"
          choices={[
            { id: 'default', label: 'Default' },
            { id: 'cartoon', label: 'Cartoon' },
            { id: 'spacefill', label: 'Spacefill' },
            { id: 'surface', label: 'Molecular surface' },
          ]}
          onPick={() => {}}
        />
      </div>
    )),
};
