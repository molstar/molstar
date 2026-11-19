/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { TrnaViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = (id: string) => TrnaViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Examples/tRNA' };
export default meta;

/**
 * tRNA-Phe (1EHZ): its base pairs as arcs over the letters, the four stems nested and the tertiary pairs crossing them,
 * colored by kind (canonical, wobble, non-canonical). Hover a base: the arc to its partner lights. Click an arc: both
 * bases are selected. Pan: an arc whose end leaves the window moves continuously.
 */
export const Arcs: StoryObj = { render: () => exampleStory(view('arcs')) };

/** The largest nested set of pairs above the letters, the pairs that cross it below (`crossing: split`). */
export const SplitCrossing: StoryObj = { render: () => exampleStory(view('split')) };

/** The tRNA on a ring, its pairs as chords through the hole (`side: below`). */
export const Chords: StoryObj = { render: () => exampleStory(view('ring')) };

/**
 * The base-pair dot plot under the arcs, its axes following their window: probabilities above the diagonal, the
 * structure's pairs mirrored below. Hover a cell: a crosshair and both bases lit on the letters. Drag a rectangle:
 * both ranges selected; Cmd/Ctrl+drag adds one, Alt+drag takes one out.
 */
export const DotPlot: StoryObj = { render: () => exampleStory(view('dot-plot')) };
