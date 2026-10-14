/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from './support/seqstar-story.js';
import { AtlasPage } from './support/seqstar-atlas.js';

/**
 * Performance of the Canvas2D backend on large documents. The glyph atlas on a large lazy alignment:
 * frame times with it on and off, and a magnified crop of the same letters either way.
 */
const meta = {
  title: 'Seq*/Performance/Large documents',
  parameters: { layout: 'fullscreen', docs: { source: { type: 'code' } } },
} satisfies Meta;
export default meta;
type Story = StoryObj<typeof meta>;

export const GlyphAtlas: Story = {
  name: 'Glyph atlas',
  render: () => <AtlasPage />,
};
