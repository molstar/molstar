/** @jsxImportSource @molstar/plugin-prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../_stories/story-root.js';
import { observableState } from '../../reactive/observable-state.js';
import { TextInput } from './text-input.js';

const meta: Meta = { title: 'Components/Text input' };
export default meta;

const panel =
  'display: grid; gap: var(--msp-space-3); width: 300px; padding: var(--msp-space-4); background: var(--msp-color-panel); color: var(--msp-color-text); font-size: var(--msp-font-size-sm);';
const readout = 'color: var(--msp-color-text-muted);';

/** The small and the medium field, with a placeholder; the value follows the typing. */
export const Sizes: StoryObj = {
  render: () =>
    story(() => {
      const sm = observableState('');
      const md = observableState('1tqn');
      return (
        <div style={panel}>
          <TextInput
            size="sm"
            placeholder="Small"
            aria-label="Small"
            value={sm.asReadonly}
            onInput={(v) => sm.set(v)}
          />
          <span style={readout}>value: {sm.asReadonly}</span>
          <TextInput
            size="md"
            placeholder="Medium"
            aria-label="Medium"
            value={md.asReadonly}
            onInput={(v) => md.set(v)}
          />
          <span style={readout}>value: {md.asReadonly}</span>
        </div>
      );
    }),
};

/** Enter submits: the submitted value is echoed below, apart from the live one. */
export const Submit: StoryObj = {
  render: () =>
    story(() => {
      const submitted = observableState('');
      return (
        <div style={panel}>
          <TextInput placeholder="PDB ID, then Enter" aria-label="PDB ID" onSubmit={(v) => submitted.set(v)} />
          <span style={readout}>submitted: {submitted.asReadonly}</span>
        </div>
      );
    }),
};
