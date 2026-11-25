/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ExampleView, GlobinsViews, P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { liveStory } from './live.js';
import * as smallSpec from './small-spec.js';
import smallSpecSource from './small-spec.ts?raw';

// The Docs pages' live examples, each a page's `<Canvas of>`: not listed in the sidebar (`!dev`), the stories they
// illustrate are (each page links them). Light and dark by the toolbar's theme.

const meta: Meta = { title: 'Docs/Live', tags: ['!dev'] };
export default meta;

const view = (views: readonly ExampleView[], id: string) => views.find((v) => v.id === id)!;

/** p53's letters, its domains in lanes and a conservation signal: the Introduction's view. */
export const Tracks: StoryObj = { render: () => liveStory(view(P53Views, 'tracks')) };

/** The gene's section following the protein's through the codon map: Concepts' coordinate systems and maps. */
export const ProteinOnItsGene: StoryObj = { render: () => liveStory(view(P53Views, 'protein-on-gene')) };

/** A tracks overview and the protein as text under it, one view of two layouts: Concepts' sections and layouts. */
export const TextUnderTracks: StoryObj = { render: () => liveStory(view(P53Views, 'text-under-tracks')) };

/** The language page's spec: the code, the view, and the spec printed as an S-expression and as JSON. */
export const ASmallSpec: StoryObj = {
  render: () =>
    liveStory(
      {
        dataset: 'docs',
        id: 'small-spec',
        title: 'A small spec',
        summary: 'A protein, its domains and a signal, one section of three tracks.',
        build: smallSpec.build,
      },
      { code: smallSpecSource.replace(/^\/\*\*[\s\S]*?\*\/\s*/, ''), printed: ['sexpr', 'json'] },
    ),
};

/** The gene's region read region by region, transcripts arriving as nodes: the Data page's pending and arriving. */
export const GeneByRegion: StoryObj = { render: () => liveStory(view(P53Views, 'gene-by-region')) };

/** The canvas and the SVG export of one spec, side by side: the Embedding page's two backends. */
export const CanvasAndSvg: StoryObj = { render: () => liveStory(view(P53Views, 'tracks'), { svg: true }) };

/** Ten million alignment rows, only those in view ever made: the Performance page's cost by pixels. */
export const TenMillionRows: StoryObj = { render: () => liveStory(view(GlobinsViews, 'lazy')) };
