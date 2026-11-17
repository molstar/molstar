/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../../_stories/support.js';
import { state } from '../../../state.js';
import { Markdown } from '../markdown.js';

const meta: Meta = { title: 'Components/Markdown' };
export default meta;

// The text at the inset from the surface's edge (styles §5).
const panel = (width: number) =>
  `display: grid; gap: var(--msui-space-3); width: ${width}px; padding: var(--msui-inset); border-radius: var(--msui-radius-lg); background: var(--msui-surface); color: var(--msui-text);`;
const readout = 'color: var(--msui-text-muted);';

const source = `## Binding site

This snapshot shows the **heme** pocket of *CYP3A4*, colored by \`element-symbol\`.

- Select the ligand
- Focus it with [the camera](#focus-ligand)
- Compare with [the apo state](#snapshot-2)

1. Load the structure
2. Apply the style

| Residue | Atoms |
|:--|--:|
| **HEM** | 43 |
| CYS 442 | 6 |

\`\`\`
select resn HEM
\`\`\`

---

See the [PDB entry](https://www.rcsb.org/structure/1TQN) for details.`;

/** A help text with a heading, a list, code and links: `#` links go to `onLink`, web links open in a new tab. */
export const HelpText: StoryObj = {
  render: () =>
    story(() => {
      const link = state('none');
      return (
        <div style={panel(300)}>
          <Markdown source={source} onLink={(href) => link.set(href)} />
          <span style={readout}>link: {link.asReadonly}</span>
        </div>
      );
    }),
};

/** Reading text (`size="md"`, 14/20): a help page or a long description. */
export const Reading: StoryObj = {
  render: () =>
    story(() => (
      <div style={panel(420)}>
        <Markdown source={source} size="md" />
      </div>
    )),
};
