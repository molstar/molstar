/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/react-vite';
import { useMemo } from 'react';
import { asText, exonsOnly, letterMapStrands } from './support/seqstar-docs-views.js';
import f from './support/seqstar-format.module.css';
import { PrismPage } from './support/seqstar-prism.js';
import { ViewerCard } from './support/seqstar-viewer.js';

/**
 * Letter maps (`letter_map`): a layer's table from the stored letter to
 * the one shown, applied where the layer draws a letter and where its theme reads one. The data keeps its letters:
 * search, mapping statuses and the model are the stored ones, and hovering a position names the letter shown and the
 * stored one. A map is a name (`complement`, `complement_rna`, `lower`) or single letters mapped; "complement" is one
 * table among others, and the core knows no strand: the strand flip is `direction: -1` and a map, composed by
 * the builder's `reverseComplement`.
 */
const meta = { title: 'Seq*/Views and layouts/Axes and strands', parameters: { layout: 'fullscreen' } } satisfies Meta;
export default meta;
type Story = StoryObj;

function Cards() {
  const flip = useMemo(letterMapStrands, []);
  const exons = useMemo(exonsOnly, []);
  const text = useMemo(asText, []);
  return (
    <div className={f.stack}>
      <ViewerCard
        doc={flip}
        label="Reverse complement, and mRNA (T shown as U) over the codons"
        options={{ parts: ['SectionHeader'] }}
      />
      <ViewerCard
        doc={exons}
        label="Exons only: their positions, the mRNA map and the codons across the intron's break"
        options={{ parts: ['SectionHeader'] }}
      />
      <ViewerCard doc={text} label="The reverse complement as text" />
    </div>
  );
}

export const Maps: Story = {
  name: 'Letter maps',
  render: () => (
    <PrismPage
      title="Letters through a map: the other strand, the mRNA, the exons"
      caption={
        <>
          Case 14's gene fragment three ways. <b>Reverse complement</b> is <code>direction: -1</code> (the layout runs
          from the highest position down) and <code>letter_map: "complement"</code> on its letters, the builder's{' '}
          <code>b.preset.reverseComplement(section)</code>: the ruler counts down, the letters are the other strand's,
          and their cells are colored by the letter shown (A is green, T blue); the exon labels mirror, and the
          emphasized start codon reads <code>CAT</code>, as its track shows it. <b>mRNA</b> is an inline map,{' '}
          <code>{'{ T: "U" }'}</code>, on the same sequence under its protein's codons. Hover a position: the label
          names the letter shown and, where it differs, the stored one. <b>Exons only</b> is <code>positions</code> over
          the exon annotation (<code>b.preset.exonsOnly</code>); the intron is a break and the codon spanning it stays
          one residue. The last card is the reverse complement as text.
        </>
      }
    >
      <Cards />
    </PrismPage>
  ),
};
