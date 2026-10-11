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
import { Markdown } from './markdown.js';

const meta: Meta = { title: 'Components/Markdown' };
export default meta;

const panel =
  'display: grid; gap: var(--msp-space-3); width: 300px; padding: var(--msp-space-4); background: var(--msp-color-panel); color: var(--msp-color-text); font-size: var(--msp-font-size-sm);';
const readout = 'color: var(--msp-color-text-muted);';

const source = `## Binding site

This snapshot shows the **heme** pocket of *CYP3A4*, colored by \`element-symbol\`.

- Select the ligand
- Focus it with [the camera](#focus-ligand)
- Compare with [the apo state](#snapshot-2)

1. Load the structure
2. Apply the style

See the [PDB entry](https://www.rcsb.org/structure/1TQN) for details.`;

/** A help text with a heading, a list, code and links: `#` links go to `onLink`, web links open in a new tab. */
export const HelpText: StoryObj = {
  render: () =>
    story(() => {
      const link = observableState('none');
      return (
        <div style={panel}>
          <Markdown source={source} onLink={(href) => link.set(href)} />
          <span style={readout}>link: {link.asReadonly}</span>
        </div>
      );
    }),
};
