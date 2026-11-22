/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RegionReadsViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';

const view = (id: string) => RegionReadsViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Features/Region reads/Pairs' };
export default meta;

/**
 * A generated synteny alignment between two toy genomes, on one section of both assemblies kept to `chrA2` and `chrB2`
 * (the section’s `include`): pairs over the two assemblies (in the `seq` format), read a chromosome at a time, each
 * found from either end, about 80 in view. Under 4,096 the pairs layer reads the items and draws each alignment block
 * as one arc between the genomes; B reads a 6 Mb stretch reversed, its arcs twisted. The data waits a moment, as a
 * remote source’s would (EC-208, EC-210).
 */
export const ZoomedIn: StoryObj = { render: () => exampleStory(view('pairs-zoomed-in'), { stacked: true }) };

/**
 * 12,000 generated links on a 2 Mb chromosome, twice: held in the spec (the first track) and read from a BEDPE file’s
 * text (the second, by a host adapter serving the file’s pairs by region). Past the 4,096 arcs a view draws, the BEDPE
 * track reads bin pairs and strokes each arc by its count (the pairs level, made by the store’s summaries); the held
 * track, not read through the store, bins arcs by their ends’ device columns, a density strip under. Hover a bin pair:
 * its card names both bins, the count and the bins’ width, never a pair; zoom in until fewer than 4,096 remain and the
 * items themselves are drawn (EC-171, EC-208).
 */
export const ZoomedOut: StoryObj = { render: () => exampleStory(view('pairs-zoomed-out'), { stacked: true }) };

/**
 * The genome ring at the fit: 12,000 generated links over a toy assembly’s six chromosomes (a `seq` source with three
 * levels), read by region a chromosome at a time, as chords between the chromosomes and within each (core’s pairs
 * layer, no extension). Past the 4,096 arcs a view draws, the arcs are binned by their ends’ device columns (EC-208).
 */
export const OnARing: StoryObj = { render: () => exampleStory(view('pairs-ring'), { stacked: true }) };
