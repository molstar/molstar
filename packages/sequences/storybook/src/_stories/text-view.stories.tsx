/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from './support/seqstar-story.js';
import { useMemo, useState } from 'react';
import { PublicationAnnotations, PublicationControls } from './support/seqstar-annotations.js';
import { ExportCard } from './support/seqstar-export.js';
import { PrismPage } from './support/seqstar-prism.js';
import { type PublicationView, publicationTracksDoc } from './support/seqstar-publication.js';

/**
 * The text view with tracks drawn in full (seqspec:viewer#the-text-view): a track with a `height` takes that many of
 * the text's rows under every line, and draws its full encoding there, wrapping with the text.
 */
const meta = { title: 'Seq*/Text view', parameters: { layout: 'fullscreen' } } satisfies Meta;
export default meta;
type Story = StoryObj;

function Publication() {
  const doc = useMemo(publicationTracksDoc, []);
  const [view, setView] = useState<PublicationView>('text');
  return (
    <ExportCard
      doc={doc}
      label={`p53 (UniProt P04637) as ${view === 'text' ? 'text, its tracks in full' : 'tracks'}`}
      controls={(m) => <PublicationControls m={m} view={view} onView={setView} />}
      aside={(m) => <PublicationAnnotations m={m} />}
      note="Hydrophobicity (Kyte–Doolittle over 9 residues, from the sequence) as a line, 2 rows; a pLDDT-like confidence (illustrative) as an area in AlphaFold DB's bands, 1.5 rows; the Zn²⁺-binding residues as glyphs; UniProt's regions as labeled ranges; the DNA-binding domain's strands and helices (illustrative). One y scale for every line; the axes and the labels are the section style's track_axes and track_labels. Text and Tracks are two presets of one section of one document (the document is not rebuilt): an annotation added in one is in the other. The Annotations panel keeps two bookmark sets, each with its own ranges track: Above is ordered before the sequence's track, Below after it. Export the figure as SVG or PNG: it is of the view shown."
    />
  );
}

export const PublicationTracks: Story = {
  name: 'Publication tracks',
  render: () => (
    <PrismPage
      title="Publication tracks in the text view"
      caption="Tracks with a height draw in full under every line of text: a values line or area on one scale across lines (a line breaks at a line's end and starts again on the next), glyphs, labeled ranges, secondary structure. Without a height a track keeps its thin lane (an underline or a strip). The Tracks button draws the same document, tracks at their heights, on one line: two views of the same data. Select a range, give it a label and a color, choose Above or Below the sequence and a line style (solid, dashed, dotted, wavy, double, box or highlight) and Add from selection (the panel is beside the viewer): the section's annotations track over the sequence or under it shows it as an underline in the text and as a ranges row in the tracks."
    >
      <Publication />
    </PrismPage>
  ),
};
