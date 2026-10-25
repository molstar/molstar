/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// An ESPript-like figure of the core's globin fixture (`pf00042-seed-8.aln`, Pfam's PF00042 seed, its first eight rows
// in CLUSTAL W's layout), from its text: the rows as text (seqspec:viewer#rows-in-text) on `espript` cells, HBA_HUMAN's
// helices above the block (approximate, after 2HHB's), an illustrative strand pair under it, glyph markers above and
// below, a box on one row's range and two titled frames over the rows (seqspec:viewer#frames). No Node import: the
// Storybook's Alignment figure story builds its document here too (from the fixture's `?raw` text).

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';

/** The fixture's rows: name (`ID/first-last`), gapped residues, and each block's printed end numbers. */
export interface AlnRow {
  readonly name: string;
  readonly gapped: string;
  readonly ends: readonly number[];
}
/** CLUSTAL W's blocks read back (`text`): a row per name, its segments joined, the number printed after each. */
export function parseClustal(text: string): AlnRow[] {
  const rows = new Map<string, { gapped: string; ends: number[] }>();
  for (const line of text.split(/\r?\n/)) {
    const m = /^(\S+)\s+([A-Za-z.-]+)\s+(\d+)$/.exec(line);
    if (!m || line.startsWith('CLUSTAL')) continue;
    const r = rows.get(m[1]) ?? { gapped: '', ends: [] };
    r.gapped += m[2];
    r.ends.push(Number(m[3]));
    rows.set(m[1], r);
  }
  return [...rows].map(([name, r]) => ({ name, gapped: r.gapped, ends: r.ends }));
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

/** A frame's row in the figure's frames set (a bookmark's: target, range, id, label, color, line style). */
export interface FrameRow {
  readonly target: string;
  readonly start: number;
  readonly stop: number;
  readonly id: string;
  readonly label: string;
  readonly color: string;
  readonly line_style?: string;
}
/** The figure's two frames: one over every row's first 20 columns, one dashed over columns 88 to 112. */
export const FIGURE_FRAMES: readonly FrameRow[] = [
  { target: 'globins', start: 0, stop: 20, id: 'f1', label: 'N-segment', color: '#9467bd' },
  { target: 'globins', start: 88, stop: 112, id: 'f2', label: 'C–E region', color: '#ff7f0e', line_style: 'dashed' },
];

/** The figure's sets and tracks; `frames` in place of `FIGURE_FRAMES`; `style` over the section's. */
export interface FigureOptions {
  /** In place of `FIGURE_FRAMES`. */
  readonly frames?: readonly FrameRow[];
  /** Over the section's style. */
  readonly style?: Record<string, unknown>;
  /** A Newick tree over the rows. */
  readonly tree?: string;
  /** The rows' color as the `espript` scheme (its dark variant on a dark skin), not the theme. */
  readonly scheme?: boolean;
  /** The alignment's reference row (row 0 by default). */
  readonly reference?: number;
  /** In place of the box on row 2: the rows of the boxes set (a range each, on a row, with its line style). */
  readonly boxes?: readonly Record<string, unknown>[];
  /** Markers above as a `points` layer of the rows track too (a layer with no `behindText` there). */
  readonly rowPoints?: boolean;
  readonly title?: string;
}
/** The figure over `rows` (`parseClustal`'s), its sets holding the opening annotations. */
export function figureFromRows(rows: readonly AlnRow[], o: FigureOptions = {}): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: o.title ?? 'Globins, ESPript-like' });
  const aln = b.data.alignment({
    ref: 'globins',
    value: rows.map((r) => r.gapped),
    labels: rows.map((r) => r.name),
    ...(o.reference !== undefined ? { reference_row: o.reference } : {}),
  } as never);
  if (o.tree) (aln as unknown as { tree(p: unknown): void }).tree({ ref: 'globins_tree', value: o.tree });
  // The label is what a hover names the set by ("Region C–E region"): a reader's word, not the ref.
  const ann = (ref: string, label: string, value: unknown) =>
    b.data.annotation({ ref, label, geometry: 'ranges', value } as never);
  const first = 27;
  const top = ann('ss.top', 'Secondary structure', {
    target: HBA_HELICES.map(() => 'globins/0'),
    start: HBA_HELICES.map((h) => h[0] - first),
    stop: HBA_HELICES.map((h) => h[1] - first + 1),
    category: HBA_HELICES.map(() => 'helix'),
    label: HBA_HELICES.map((_, k) => `α${k + 1}`),
  });
  const bottom = ann('ss.bottom', 'Secondary structure', {
    target: ['globins/3', 'globins/3', 'globins/3'],
    start: [18, 28, 44],
    stop: [24, 34, 64],
    category: ['strand', 'strand', 'helix'],
    label: ['β1', 'β2', 'α1'],
  });
  const above = ann('figure.main.markers.above', 'Marker', [
    { target: 'globins/0', start: 32, stop: 33, id: 'm1', label: 'His59 (distal)', color: '#1f77b4', glyph: 'star' },
    {
      target: 'globins/0',
      start: 61,
      stop: 62,
      id: 'm2',
      label: 'His88 (proximal)',
      color: '#d62728',
      glyph: 'triangle',
    },
  ]);
  const below = ann('figure.main.markers.below', 'Marker', [
    { target: 'globins', start: 14, stop: 15, id: 'm3', label: 'Conserved', color: '#2ca02c', glyph: 'circle' },
    { target: 'globins', start: 75, stop: 76, id: 'm4', label: 'Conserved', color: '#2ca02c', glyph: 'circle' },
  ]);
  const boxes = ann(
    'figure.main.boxes',
    'Box',
    o.boxes ?? [
      { target: 'globins/2', start: 10, stop: 18, id: 'b1', label: 'Motif', color: '#e3120b', line_style: 'box' },
    ],
  );
  const frames = ann('figure.main.frames', 'Region', o.frames ?? FIGURE_FRAMES);
  const view = b.view({ ref: 'figure', label: 'Figure', layout: 'text' });
  const sec = view.section({
    ref: 'figure.main',
    style: { preset: 'figure', text_cells: true, track_axes: false, track_labels: false, ...o.style },
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
  const espript = o.scheme ? { theme: 'scheme', params: { name: 'espript' } } : { theme: 'espript' };
  msa.rows({ letter_map: { '-': '.' }, color: espript } as never);
  msa.ranges({ data: boxes, color, line_style: style } as never);
  if (o.rowPoints) msa.points({ data: above, glyph, color } as never);
  sec.track({ ref: 'markers.below', label: 'Markers below', height: 1 }).points({ data: below, glyph, color } as never);
  sec
    .track({ ref: 'ss.bottom.track', label: 'BRUPA (illustrative)', height: 1 })
    .secondary({ data: bottom, labels: true } as never);
  return b.getDocument();
}

/** Render props giving the figure a line per 60 columns, as the fixture's own blocks: one group of 60 a line. */
export function sixtyColumns(cell: number) {
  return { wrappedGroup: 60, wrappedMaxWidth: Math.ceil(60 * cell) + 2, wrappedGroupLabels: 'none' as const };
}

/** A tree over the eight rows (the storybook's, written for the story): leaves by the rows' names before the `/`. */
export const GLOBIN_NWK =
  '((HBA_HUMAN:0.12,(A0A096MFI1_POEFO:0.21,A0A1I7SGZ4_BURXY:0.33):0.08):0.1,((A0A010R001_9PEZI:0.25,A0A178ZLR6_9EURO:0.22):0.14,(A0A0N4SXE9_BRUPA:0.3,(A0A1Y1I094_KLENI:0.18,A0A2I0WAL8_9ASPA:0.27):0.09):0.11):0.1);';
