/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { TrnaScaleView, TrnaViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = (id: string) => TrnaViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Examples/tRNA/Views' };
export default meta;

/**
 * tRNA-Phe (PDB 1EHZ, chain A): its 27 base pairs as arcs over the letters, the four stems nested and the six tertiary
 * pairs among them (four cross the stems' pairs), colored by kind (canonical, wobble, non-canonical); its 14 modified
 * residues (11 components) as their parent letters, each with a dot at its cell's top right (hover: `U · PSU` and its
 * name), their codes also in a track; the bookmarks list's “Modified residues” runs `(elements :modified true)` and
 * highlights the 14. Hover a base: the arc to its partner lights. Click an arc: both bases are selected. Pan: an arc whose end leaves the window moves continuously. Editing: hide the tertiary pairs (a filter by
 * that field in the arcs' `data`, a patch).
 */
export const Arcs: StoryObj = { render: () => exampleStory(view('arcs')) };

/** The largest nested set of pairs above the letters, the pairs that cross it below (`crossing: split`). */
export const SplitCrossing: StoryObj = { render: () => exampleStory(view('split')) };

/** The tRNA on a ring, its pairs as chords through the hole (`side: below`). */
export const Chords: StoryObj = { render: () => exampleStory(view('ring')) };

/**
 * The base-pair dot plot under the arcs, its axes following their window: above the diagonal the family's covariation
 * on 1EHZ's positions (computed from the Rfam seed), below it the structure's pairs mirrored. Hover a cell: a
 * crosshair and both bases lit on the letters. Drag a rectangle: both ranges selected; Cmd/Ctrl+drag adds one,
 * Alt+drag takes one out. Editing: change the covariation's palette (a theme patch).
 */
export const DotPlot: StoryObj = { render: () => exampleStory(view('dot-plot')) };

/** The seed row 1EHZ is placed on, with the seed's `SS_cons` pairs projected onto it: nested only. */
export const SeedRow: StoryObj = { render: () => exampleStory(view('seed-row')) };

/**
 * Rfam RF00005's seed, 954 rows of 118 columns, under its conservation logo. The `SS_cons` arc and the link from
 * 1EHZ's residues need a pairs node and a map on an alignment's columns, which the language does not have yet.
 */
export const Seed: StoryObj = { render: () => exampleStory(view('seed')) };

/** The dot plot at 4,000 × 4,000: a generated grid, the plot alone, at the frame target. */
export const DotPlotScale: StoryObj = { render: () => exampleStory(TrnaScaleView) };
