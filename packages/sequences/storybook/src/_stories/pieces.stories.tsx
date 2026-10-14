/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from './support/seqstar-story.js';
import { useMemo } from 'react';
import { flipped } from './support/seqstar-docs-views.js';
import f from './support/seqstar-format.module.css';
import { PrismPage } from './support/seqstar-prism.js';
import { ViewerCard } from './support/seqstar-viewer.js';

/**
 * Pieces: what a block shows of its coordinate system, a `range` of it,
 * laid out in a `direction`. A reversed piece is a layout map, never new data: every layer draws it mirrored, the ruler
 * descends, and hover and selection reach the same positions in every section that shows them.
 */
const meta = { title: 'Seq*/Views and layouts/Axes and strands', parameters: { layout: 'fullscreen' } } satisfies Meta;
export default meta;
type Story = StoryObj;

function Cards() {
  const doc = useMemo(flipped, []);
  return (
    <div className={f.stack}>
      <ViewerCard
        doc={doc}
        label="The construct, its reverse strand and the GFP fragment reversed"
        options={{ parts: ['SectionHeader'] }}
      />
      <ViewerCard doc={doc} label="The reverse strand as text" options={{ view: 'seq' }} />
    </div>
  );
}

export const Direction: Story = {
  name: 'Direction of a piece',
  render: () => (
    <PrismPage
      title="Case 5, flipped: pieces with range and direction"
      caption={
        <>
          Three sections over one insert. As written: the ruler counts up, EcoRI is cut at 4 and XhoI at 64. Reverse
          strand (direction −1): the same blocks laid out from the last nucleotide down, so the ruler counts down, the
          features mirror and each cut site sits on the other side of its nucleotide. The GFP fragment (range [9, 63],
          direction −1): one piece of the insert, numbered as the insert. Hover or select in any section: the same
          nucleotides light in the other two (one coordinate system, three layout maps). The letters are the forward
          strand's read right to left; complementing them is a letter map's.
        </>
      }
    >
      <Cards />
    </PrismPage>
  ),
};
