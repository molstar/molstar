/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/react-vite';
import { QueryStory } from './support/seqstar-query-runtime.js';

/**
 * SeqQL's runtime: a gallery of example queries beside a mounted
 * viewer, over cases 1, 3, 7, 13 and 14 and a barnase family bound to 1BRS, each with its S-expression, what `sel`
 * expands to, and a line of English. Clicking one runs it on the viewer's set (`@molstar/seq-core/model/query`): the result
 * becomes the viewer's selection, or a highlight with the toggle, and the summary gives its groups, positions and
 * systems, each run's labels with its slice; Previous and Next step through the result's groups one at a time (a
 * search hit's score with it). "Show me" focuses on the result's extent (or the group stepped to), switching to the
 * view that shows it when the current one shows none (`showExtents`). The box takes any query; Ctrl+Enter runs it, and
 * issues show under the line and column they point at, the runtime's included, an empty result's diagnosis among them
 * (case 3's `(values :of cons :above 1)`). Every base symbol runs: the scans (`letters`, `positions`, `search` with the
 * story's letter-class tables), the position, system and group tests among them.
 */
const meta = { title: 'Seq*/Query and search/SeqQL', parameters: { layout: 'fullscreen' } } satisfies Meta;
export default meta;
type Story = StoryObj<typeof meta>;

export const Gallery: Story = {
  name: 'Runtime gallery',
  render: () => <QueryStory />,
};
