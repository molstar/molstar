/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SeqCoreConfig } from '@molstar/seqx-core/context';
import type { RailTint } from '@molstar/seqx-core/render';
import { FeatureViews, GenomeFormatsViews, P53Views } from '@molstar/seqx-examples';
import type { ExampleView } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// The rail looks story's composition, DOM-free: the story and its headless twin call these. A look is the render
// config's `navigator.rail.tint`; each cell is one rail of
// the one rail path (the ruler's, the rows', the lanes' band's) in a view of the examples.

/** A tint of the rail, with the sentence the story prints under it. */
export interface LookEntry {
  readonly look: RailTint;
  readonly title: string;
  readonly says: string;
}

export const Looks: readonly LookEntry[] = [
  {
    look: 'neutral',
    title: 'Neutral (the default)',
    says: 'A hairline the axis’s length and the window a heavier segment in a mid grey (the text ink at 55%), a round knob at each end where the box zooms; full ink and grown while used; at the fit the hairline’s weight and small knobs; pan-only rails (lanes) round caps.',
  },
  {
    look: 'accent',
    title: 'Accent (`navigator.rail.tint: accent`)',
    says: 'The same with the segment and knobs in a faint accent tint (the accent at 55%).',
  },
];

/** A cell: a rail in a view, and where the view opens (a function of the viewer; none: at its fit). */
export interface CellEntry {
  readonly id: 'ruler-zoomed' | 'ruler-fit' | 'rows' | 'lanes';
  readonly title: string;
  readonly view: ExampleView;
  readonly opening?: (viewer: SeqViewer) => void;
}

const p53 = P53Views.find((v) => v.id === 'tracks')!;

export const Cells: readonly CellEntry[] = [
  {
    id: 'ruler-zoomed',
    title: 'Ruler rail, zoomed in',
    view: p53,
    opening: (v) => v.view.setRange('/views/main/protein', 94, 200),
  },
  { id: 'ruler-fit', title: 'Ruler rail, at the fit', view: p53 },
  { id: 'rows', title: 'Rows rail (an alignment)', view: FeatureViews.find((v) => v.id === 'rows-rail')! },
  { id: 'lanes', title: 'Lanes band (GFF3, TP53)', view: GenomeFormatsViews.find((v) => v.id === 'gff-tp53')! },
];

/** The viewer config item that picks a look. */
export const lookConfig = (tint: RailTint) => [SeqCoreConfig.render, { navigator: { rail: { tint } } }] as const;
