/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = P53Views.find((v) => v.id === 'gene-region-ideogram')!;

const meta: Meta = { title: 'Examples/p53/Gene region on chr17' };
export default meta;

/**
 * The gene's region with chromosome 17's ideogram above it (GRCh38's bands, the centromere two triangles): the region's
 * window on the chromosome, the rest of it washed (out of reach). Hover a band for its name; drag the window, click
 * outside it, double-click. The SVG beside it includes the ideogram, its window an outline titled with its positions.
 */
export const GeneRegionOnChr17: StoryObj = { render: () => exampleStory(view, { navigators: 'include' }) };
