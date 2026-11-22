/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SeqCoreConfig } from '@molstar/seqx-core/context';
import type { RailLook } from '@molstar/seqx-core/render';
import { FeatureViews, GenomeFormatsViews, P53Views } from '@molstar/seqx-examples';
import type { ExampleView } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// The rail looks story's composition, DOM-free: the story and its headless twin call these. A look is the render
// config's `navigator.rail.look` (internal: the story's alone until the author picks one); each cell is one rail of
// the one rail path (the ruler's, the rows', the lanes' band's) in a view of the examples.

/** A look: today's and the three variants, with the sentence the story prints under it. */
export interface LookEntry {
  readonly look: RailLook;
  readonly title: string;
  readonly says: string;
}

export const Looks: readonly LookEntry[] = [
  {
    look: 'current',
    title: 'Today',
    says: 'A dark line with the window a thicker segment and a square tick at each end; the ticks overhang the well’s rounded corners.',
  },
  {
    look: 'capsule',
    title: 'A: capsule',
    says: 'The window a rounded pill in a soft accent, no groove; its ends are the grips (no ticks), taking the ink and growing under the pointer.',
  },
  {
    look: 'groove',
    title: 'B: capsule in a groove',
    says: 'The same pill the thumb of a faint rounded groove the axis’s length, as a native scroll track is.',
  },
  {
    look: 'quiet',
    title: 'C: quiet until hovered',
    says: 'B thin and low at rest, in the full ink and grown while hovered or dragged; at the fit only the groove.',
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
export const lookConfig = (look: RailLook) => [SeqCoreConfig.render, { navigator: { rail: { look } } }] as const;
