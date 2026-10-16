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

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import SEED_ALN from '../../../../core/_test/fixtures/pf00042-seed-8.aln?raw';

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

/** CLUSTAL W's blocks, read back: a row per name (`ID/first-last`), its segments joined. */
export function figureRows(): { name: string; gapped: string }[] {
  const rows = new Map<string, string>();
  for (const line of SEED_ALN.split('\n')) {
    const m = /^(\S+)\s+([A-Za-z.-]+)\s+(\d+)$/.exec(line);
    if (m && !line.startsWith('CLUSTAL')) rows.set(m[1], (rows.get(m[1]) ?? '') + m[2]);
  }
  return [...rows].map(([name, gapped]) => ({ name, gapped }));
}

/** HBA_HUMAN's helices in UniProt numbering (2HHB chain A's, DSSP, +1; approximate), within the row's 27–137. */
const HBA_HELICES: readonly (readonly [number, number])[] = [
  [27, 36],
  [37, 43],
  [53, 72],
  [81, 89],
  [95, 112],
  [119, 137],
];

/** The figure; its sets hold the opening annotations, the panel adds to them. */
export function alignmentFigureDoc(): SeqViewSpecDoc {
  const rows = figureRows();
  const b = createSeqViewSpec({ title: 'Globins (PF00042): an ESPript-like figure' });
  const aln = b.data.alignment({
    ref: 'globins',
    label: 'PF00042 seed (eight rows)',
    value: rows.map((r) => r.gapped),
    labels: rows.map((r) => r.name),
  } as never);
  const ann = (ref: string, label: string, value: unknown) =>
    b.data.annotation({ ref, label, geometry: 'ranges', value } as never);
  const first = 27;
  const top = ann('ss.top', 'HBA_HUMAN helices', {
    target: HBA_HELICES.map(() => 'globins/0'),
    start: HBA_HELICES.map((h) => h[0] - first),
    stop: HBA_HELICES.map((h) => h[1] - first + 1),
    category: HBA_HELICES.map(() => 'helix'),
    label: HBA_HELICES.map((_, k) => `α${k + 1}`),
  });
  const bottom = ann('ss.bottom', 'BRUPA (illustrative)', {
    target: ['globins/3', 'globins/3', 'globins/3'],
    start: [18, 28, 44],
    stop: [24, 34, 64],
    category: ['strand', 'strand', 'helix'],
    label: ['β1', 'β2', 'α1'],
  });
  const above = ann('figure.main.markers.above', 'Markers above', [
    { target: 'globins/0', start: 32, stop: 33, id: 'm1', label: 'His58 (distal)', color: '#1f77b4', glyph: 'star' },
    {
      target: 'globins/0',
      start: 61,
      stop: 62,
      id: 'm2',
      label: 'His87 (proximal)',
      color: '#d62728',
      glyph: 'triangle',
    },
  ]);
  const below = ann('figure.main.markers.below', 'Markers below', [
    { target: 'globins', start: 14, stop: 15, id: 'm3', label: 'Conserved', color: '#2ca02c', glyph: 'circle' },
    { target: 'globins', start: 75, stop: 76, id: 'm4', label: 'Conserved', color: '#2ca02c', glyph: 'circle' },
  ]);
  const boxes = ann('figure.main.boxes', 'Boxes', [
    { target: 'globins/2', start: 10, stop: 18, id: 'b1', label: 'Motif', color: '#e3120b', line_style: 'box' },
  ]);
  const frames = ann('figure.main.frames', 'Frames', [
    { target: 'globins', start: 0, stop: 20, id: 'f1', label: 'N-segment', color: '#9467bd' },
    { target: 'globins', start: 88, stop: 112, id: 'f2', label: 'C–E region', color: '#ff7f0e', line_style: 'dashed' },
  ]);
  const view = b.view({ ref: 'figure', label: 'Figure', layout: 'text' });
  const sec = view.section({
    ref: FIGURE_SECTION,
    style: { preset: 'figure', text_cells: true, track_axes: false, track_labels: false },
  } as never);
  sec.block({ alignment: aln } as never);
  const glyph = { theme: 'field_glyph', params: { field: 'glyph' } };
  const color = { theme: 'field_color' };
  const style = { theme: 'field_line_style', params: { field: 'line_style' } };
  sec
    .track({ ref: 'ss.top.track', label: 'HBA_HUMAN (2HHB, approximate)', height: 1 })
    .secondary({ data: top, labels: true } as never);
  sec.track({ ref: 'markers.above', label: 'Markers above', height: 1 }).points({ data: above, glyph, color } as never);
  sec
    .track({ ref: 'frames', label: 'Regions' })
    .frame({ data: frames, span: { from: 'msa' }, color, line_style: style } as never);
  const msa = sec.track({ ref: 'msa', label: 'Globins' });
  // The `espript` scheme: its theme, and its dark variant (light cells) on a dark skin.
  msa.rows({ letter_map: { '-': '.' }, color: { theme: 'scheme', params: { name: 'espript' } } } as never);
  msa.ranges({ data: boxes, color, line_style: style } as never);
  sec.track({ ref: 'markers.below', label: 'Markers below', height: 1 }).points({ data: below, glyph, color } as never);
  sec
    .track({ ref: 'ss.bottom.track', label: 'BRUPA (illustrative)', height: 1 })
    .secondary({ data: bottom, labels: true } as never);
  return b.getDocument();
}
