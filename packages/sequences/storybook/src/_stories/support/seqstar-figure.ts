/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Text view › Alignment figure story (`text-view.stories.tsx`): an ESPript-like figure of the core's globin fixture
// (`pf00042-seed-8.aln`, Pfam's PF00042 seed, its first eight rows): the rows as text (seqspec:viewer#rows-in-text) on
// the `espript` theme's cells, HBA_HUMAN's helices above the block (approximate: 2HHB chain A's, DSSP, in UniProt
// numbering) and an illustrative strand pair under it, glyph markers above and below, a box on one row's range and two
// titled frames over every row (seqspec:viewer#frames). The panel's places are the figure's bookmark sets: a frame, a
// box on a row, markers above and below. No React here (the smoke test imports it).

import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import SEED_ALN from '../../../../core/_test/fixtures/pf00042-seed-8.aln?raw';
import { figureFromRows, parseClustal } from '../../../../viewer/_test/figure-docs.js';

/** The section of the figure: the text view and the tracks view (`msa`) are two presets of it. */
export const FIGURE_SECTION = 'figure.main';
export const FIGURE_FRAMES = `${FIGURE_SECTION}.frames`;
export const FIGURE_BOXES = `${FIGURE_SECTION}.boxes`;
export const FIGURE_MARKERS_ABOVE = `${FIGURE_SECTION}.markers.above`;
export const FIGURE_MARKERS_BELOW = `${FIGURE_SECTION}.markers.below`;
/**
 * The Annotations panel's places: the frames (a selection over the columns frames every row; over rows, those rows),
 * a box on a row (a selection on one row's residues), markers above and below the block (a glyph each).
 */
export const FIGURE_PLACES = [
  { data: FIGURE_FRAMES, label: 'Frame', kind: 'frame' },
  { data: FIGURE_BOXES, label: 'Box on a row', kind: 'ranges' },
  { data: FIGURE_MARKERS_ABOVE, label: 'Markers above', kind: 'points' },
  { data: FIGURE_MARKERS_BELOW, label: 'Markers below', kind: 'points' },
] as const;
/** The section's two presets: `null` is the document's own (`figure`, the text), `msa` the tracks. */
export type FigureView = 'text' | 'tracks';
export const FIGURE_PRESETS: Readonly<Record<FigureView, string | null>> = { text: null, tracks: 'msa' };

/**
 * The figure: the shared builder's (viewer/_test/figure-docs.ts) over the fixture's text, the rows on the `espript`
 * scheme's cells (its dark variant on a dark skin).
 */
export function alignmentFigureDoc(): SeqViewSpecDoc {
  return figureFromRows(parseClustal(SEED_ALN), {
    scheme: true,
    title: 'Globins (PF00042): an ESPript-like figure',
  });
}
