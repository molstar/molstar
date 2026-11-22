/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SeqCoreConfig } from '@molstar/seqx-core/context';
import type { ValuesLook } from '@molstar/seqx-core/render';
import { P53Views } from '@molstar/seqx-examples';
import type { ExampleView } from '@molstar/seqx-examples';
import { GeneLength } from '@molstar/seqx-examples/p53/data';
import { GeneOffset } from '@molstar/seqx-examples/p53/region';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// The values looks story's composition, DOM-free: the story and its headless twin call these. A look is the render
// config's `values.look` (internal, for the story); each cell is the values layer of one of p53's views, zoomed out or in.

/** A look of the values layer, with the sentence the story prints under it. */
export interface LookEntry {
  readonly look: ValuesLook;
  readonly title: string;
  readonly says: string;
}

export const Looks: readonly LookEntry[] = [
  {
    look: 'today',
    title: 'Today',
    says: 'Hard columns of full ink from the baseline; the tallest touch the well’s top; a missing value a grey stub at the floor.',
  },
  {
    look: 'profile',
    title: 'A, profile',
    says: 'Below 4 px a value one stepped area: the fill at 50%, a crisp 1 px top line in the full color; from 4 px plain bars. 10% headroom, a missing value a gap, a hairline baseline.',
  },
  {
    look: 'bars',
    title: 'B, profile and bars',
    says: 'A, and from 4 px separate bars with a 1 px gap and slightly rounded tops.',
  },
  {
    look: 'quiet',
    title: 'C, B quieter',
    says: 'B with a lighter fill (35%) and a thinner top line (0.75 px); its bars at 60%.',
  },
];

/** A cell: a values layer in a view, and where the view opens (a function of the viewer; none: at its fit). */
export interface CellEntry {
  readonly id: 'conservation-out' | 'conservation-in' | 'bigwig-out' | 'bigwig-in';
  readonly title: string;
  readonly view: ExampleView;
  readonly opening?: (viewer: SeqViewer) => void;
}

const tracks = P53Views.find((v) => v.id === 'tracks')!;
const bigwig = P53Views.find((v) => v.id === 'signal-from-a-bigwig')!;
const mid = GeneOffset + Math.floor(GeneLength / 2);

export const Cells: readonly CellEntry[] = [
  { id: 'conservation-out', title: 'Conservation, zoomed out (the whole protein)', view: tracks },
  {
    id: 'conservation-in',
    title: 'Conservation, zoomed in (residues 94 to 140)',
    view: tracks,
    opening: (v) => v.view.setRange('/views/main/protein', 94, 140),
  },
  { id: 'bigwig-out', title: 'bigWig signal, zoomed out (the region in bins)', view: bigwig },
  {
    id: 'bigwig-in',
    title: 'bigWig signal, zoomed in (50 bases)',
    view: bigwig,
    opening: (v) => v.view.setRange('/views/main/region', mid, mid + 50),
  },
];

/** The viewer config item that picks a look. */
export const lookConfig = (look: ValuesLook) => [SeqCoreConfig.render, { values: { look } }] as const;
