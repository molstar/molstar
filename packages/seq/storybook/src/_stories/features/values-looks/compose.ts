/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { ExampleView } from '@molstar/seqx-examples';
import { GeneLength } from '@molstar/seqx-examples/p53/data';
import { valuesColors } from '@molstar/seqx-examples/p53/layers';
import { GeneOffset } from '@molstar/seqx-examples/p53/region';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// The values look story's composition, DOM-free: the story and its headless twin call these. Each cell is the values
// layer of one of p53's views, zoomed out (a profile) or in (bars); the color cells show one confidence as a uniform
// series against one colored by pLDDT's bands (full strength: the legend's colors), as bars, a line and points.

/** A cell: a values layer in a view, and where the view opens (a function of the viewer; none: at its fit). */
export interface CellEntry {
  readonly id: 'conservation-out' | 'conservation-in' | 'bigwig-out' | 'bigwig-in' | 'colors-out' | 'colors-in';
  readonly title: string;
  readonly view: ExampleView;
  /** The range the title says, when the cell opens zoomed: its section and `[start, end)`. */
  readonly zoom?: { readonly section: string; readonly start: number; readonly end: number };
}

/** Moves a zoomed cell to its range; true once the view holds it (the story retries until it does). */
export function zoomed(cell: CellEntry, viewer: SeqViewer): boolean {
  const z = cell.zoom;
  if (!z) return true;
  viewer.view.setRange(z.section, z.start, z.end);
  const r = viewer.view.get(z.section);
  return r?.start === z.start && r?.end === z.end;
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
    zoom: { section: '/views/main/protein', start: 94, end: 140 },
  },
  { id: 'bigwig-out', title: 'bigWig signal, zoomed out (the region in bins)', view: bigwig },
  {
    id: 'bigwig-in',
    title: 'bigWig signal, zoomed in (50 bases)',
    view: bigwig,
    zoom: { section: '/views/main/region', start: mid, end: mid + 50 },
  },
];

const colors: ExampleView = {
  dataset: 'p53',
  id: 'values-colors',
  title: 'Uniform against by value',
  summary: 'One confidence as bars, a line and points: uniform, and colored by pLDDT’s bands.',
  build: valuesColors,
};

/** Uniform against by-value color, at both densities: zoomed out (a profile, a column's mark) and in (bars, dots). */
export const ColorCells: readonly CellEntry[] = [
  { id: 'colors-out', title: 'Uniform and by value, zoomed out (the whole protein)', view: colors },
  {
    id: 'colors-in',
    title: 'Uniform and by value, zoomed in (residues 94 to 130: 6 px or more a value, dots on the line)',
    view: colors,
    zoom: { section: '/views/main/protein', start: 94, end: 130 },
  },
];

/** Every cell of the story. */
export const AllCells: readonly CellEntry[] = [...Cells, ...ColorCells];
