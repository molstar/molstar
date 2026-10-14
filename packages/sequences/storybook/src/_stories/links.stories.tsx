/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SectionViewport } from '@molstar/seq-viewer';
import { useLiveVersion, useSeqViewer } from '@molstar/seq-viewer/react';
import type { Meta, StoryObj } from '@storybook/react-vite';
import { useMemo } from 'react';
import { codingPair, linkStrands } from './support/seqstar-docs-views.js';
import f from './support/seqstar-format.module.css';
import { PrismPage } from './support/seqstar-prism.js';
import { ViewerCard } from './support/seqstar-viewer.js';

/**
 * Axis links: a section that `follow`s another. When either moves, the
 * other centers the image of the leader's center element (with the fraction of its cell) through the route between
 * them and takes the leader's cell times the route's scale. Hover and selection already crossed the route; following
 * couples the windows.
 */
const meta = { title: 'Seq*/Views and layouts/Axes and strands', parameters: { layout: 'fullscreen' } } satisfies Meta;
export default meta;
type Story = StoryObj;

/** One section's window as the follow sees it: the element at its center (1-based, with its fraction) and its cell. */
function Center({ section, unit }: { section: string; unit: string }) {
  const { context } = useSeqViewer();
  const vp = context.viewports.get(section);
  return vp ? <CenterOf vp={vp} unit={unit} /> : <div className={f.status}>{section}: no viewport</div>;
}
function CenterOf({ vp, unit }: { vp: SectionViewport; unit: string }) {
  const out = useMemo(() => new Float64Array(2), []);
  useLiveVersion(vp.live);
  const b = vp.centerElement(out);
  const at =
    b < 0 ? 'no center' : `${unit} ${(out[0] + 1).toFixed(2)} at the center, ${out[1].toFixed(2)} px a ${unit}`;
  return (
    <div className={f.status} data-follow={vp.ref}>
      {vp.ref}: {at} (last moved by {vp.live.current.origin})
    </div>
  );
}

function Cards() {
  const pair = useMemo(codingPair, []);
  const both = useMemo(linkStrands, []);
  return (
    <div className={f.stack}>
      <div style={{ maxWidth: 560 }}>
        <ViewerCard doc={pair} label="The protein follows the gene (a codon route, a third of the cell)">
          <Center section="pair.dna" unit="nucleotide" />
          <Center section="pair.protein" unit="residue" />
        </ViewerCard>
      </div>
      <ViewerCard doc={both} label="The reverse strand follows the forward one (a reversed piece)">
        <Center section="strands.forward" unit="nucleotide" />
        <Center section="strands.reverse" unit="nucleotide" />
      </ViewerCard>
    </div>
  );
}

export const Follow: Story = {
  name: 'Following sections',
  render: () => (
    <PrismPage
      title="Two views following each other: axis links"
      caption={
        <>
          A section with <code>follow</code> and the section it names set each other's window. Drag, zoom (Ctrl or ⌘ and
          the wheel, or the navigator's edges) or scroll either one: the other centers the same element. The protein and
          its gene (case 14) follow through the codon mapping, so the gene's cell is a third of the protein's and a
          residue sits over its codon; with the gene's center in the intron the protein waits at the nearest coding
          nucleotide's residue. The card is narrow so both sections scroll at that scale (the protein's 20 residues fit
          a wide one). The reverse strand of case 5's insert follows the forward one through the identity: dragging it
          right moves the forward one left. The lines under each card read each section's center and cell; the
          follower's last move is the <code>behavior</code>'s.
        </>
      }
    >
      <Cards />
    </PrismPage>
  ),
};
