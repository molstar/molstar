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
import {
  alignmentFigureDoc,
  FIGURE_PLACES,
  FIGURE_PRESETS,
  FIGURE_SECTION,
  type FigureView,
} from './support/seqstar-figure.js';
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
      note="Hydrophobicity (Kyte–Doolittle over 9 residues, from the sequence) as a line, 2 rows; a pLDDT-like confidence (illustrative) as an area in AlphaFold DB's bands, 1.5 rows; the Zn²⁺-binding residues as glyphs; UniProt's regions as labeled ranges; the DNA-binding domain's strands and helices (illustrative). One y scale for every line; the axes and the labels are the section style's track_axes and track_labels. Text and Tracks are two presets of one section of one document (the document is not rebuilt): an annotation added in one is in the other. The Annotations panel keeps four bookmark sets, each with its own track: Above is ordered before the sequence's track, Below after it; Frame outlines the selection from the sequence down to the secondary structure, titled with its label; Markers draws a glyph each. Export the figure as SVG or PNG: it is of the view shown."
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

function Figure() {
  const doc = useMemo(alignmentFigureDoc, []);
  const [view, setView] = useState<FigureView>('text');
  return (
    <ExportCard
      doc={doc}
      label={`Globins (PF00042, eight rows) as ${view === 'text' ? 'an ESPript-like figure' : 'tracks'}`}
      controls={(m) => (
        <PublicationControls m={m} view={view} onView={setView} section={FIGURE_SECTION} presets={FIGURE_PRESETS} />
      )}
      aside={(m) => <PublicationAnnotations m={m} section={FIGURE_SECTION} places={FIGURE_PLACES} />}
      note="The core's globin fixture (Pfam PF00042's seed, its first eight rows). Shading: the espript theme (black where a column's residues are identical, grey where the residue's group holds 70% of it; gaps not counted), the letters in the ink that reads on their cells (the style's text_cells). HBA_HUMAN's helices are approximate (2HHB chain A's, in UniProt numbering); the strands under the block are illustrative (globins have none). The numbers above each line are HBA_HUMAN's (the reference row's), each row's last residue number on a line at its right. The Annotations panel adds to the figure's sets: a Frame over the selected columns frames every row (over a selection on rows, those rows), a Box on a row outlines a row's selected residues, Markers above or below draw a glyph each. Text and Tracks are two presets of one section of one document: what you add shows in both."
    />
  );
}

export const AlignmentFigure: Story = {
  name: 'Alignment figure',
  render: () => (
    <PrismPage
      title="An alignment figure in the text view"
      caption="An alignment's rows as text, ESPript-style: a block of rows per line, each row named on the left and numbered on the right, conserved columns shaded, secondary structure above and below the block as cylinders and arrows with their labels, glyph markers, a box on one row's range, and titled frames over every row that stay open where a line breaks them. The Tracks button draws the same document as an alignment's tracks, the frames over the rows."
    >
      <Figure />
    </PrismPage>
  ),
};
