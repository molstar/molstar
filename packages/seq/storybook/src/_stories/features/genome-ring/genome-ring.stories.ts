/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { rondoViewer } from '@molstar/seqx-ext-rondo/viewer';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';
import {
  partnerKeys,
  squareClicked,
  squareDragged,
  squareResized,
  squareScattered,
  stemMenuOpen,
} from './square-gestures.js';
import { hoverStrongestStem, ringStoryView } from './views.js';

// The genome ring (Design › Genome ring), `@molstar/seqx-ext-rondo`'s showcase: GRCh38's 23 chromosomes, 12,000
// generated links, UCSC's recorded cytobands, through `rondoRing`; a viewer composed with `rondoViewer` (the fallback's
// with none). A story a default of the design, each the canvas and the SVG of one view, in the toolbar's theme (light or
// dark; the last story dark).

const meta: Meta = { title: 'Features/Genome Ring' };
export default meta;

const Rondo = [rondoViewer];
const ring = (id: string, options: Parameters<typeof exampleStory>[1] = {}) =>
  exampleStory(ringStoryView(id), { stacked: true, extensions: Rondo, ...options });

/**
 * The genome at the fit: 23 arcs clockwise from the top, the cytobands outside a thin hue ring, stems as trumpets
 * inside, ribbons stem to stem wide by count, the strongest at the bottom of the stack. Hover a stem for its links,
 * double-click a ribbon to open its pair.
 */
export const TheGenomeAtTheFit: StoryObj = { render: () => ring('ring') };

/**
 * The same spec in a viewer without the extension (EC-235): the arcs and bands as with it, the links as core's bin-pair
 * chords stroked by count, no stems; the issues say once "drawn as `pairs`", with the version the spec needs.
 */
export const WithoutTheExtension: StoryObj = {
  render: () => exampleStory(ringStoryView('ring'), { stacked: true, sourceIssues: true }),
};

/**
 * One chromosome (`include: chr4`): its bands round the whole ring, its own links as ribbons (its 646 links as items
 * once their bins have landed), its outward links as stubs inward in their partners' hues, labelled at the tips.
 */
export const OneChromosome: StoryObj = { render: () => ring('ring-chromosome') };

/** A pair (chr4 × chr7, the square's cell): two arcs, their links between them as ribbons, the rest stubs. */
export const APair: StoryObj = { render: () => ring('ring-pair') };

/**
 * The square beside the ring: 23 × 23 cells by count on a square-root ramp, the diagonal the links within; its window
 * the chromosomes the ring shows. One click on a cell shows that pair at once, on a name that chromosome; drag the box
 * along the diagonal or by its corner; Cmd or Ctrl+click adds or removes; Home shows all (each one undo step).
 */
export const TheSquareBeside: StoryObj = { render: () => ring('ring', { navigators: 'include' }) };

/** The square's box dragged: the ring showed chr4 to chr6, its box moved three chromosomes on (chr7 to chr9). */
export const TheSquareDragged: StoryObj = {
  render: () => ring('ring', { navigators: 'include', opening: squareDragged }),
};

/** The square's box resized by its corner: chr4 and chr5 grown to chr4 to chr9. */
export const TheSquareResized: StoryObj = {
  render: () => ring('ring', { navigators: 'include', opening: squareResized }),
};

/** One click on the cell chr4 × chr7: the pair at once, the window that cell's run. */
export const TheSquareClicked: StoryObj = {
  render: () => ring('ring', { navigators: 'include', opening: squareClicked }),
};

/** chr4, then Cmd+click on chr9's cell: a scattered set, outlined cell by cell (its diagonal and cross cells). Dark. */
export const TheSquareScatteredDark: StoryObj = {
  globals: { theme: 'dark' },
  render: () => ring('ring', { navigators: 'include', opening: squareScattered }),
};

/**
 * The stem menu: a secondary click on the strongest stem (opened here): its partners by count shown as a pair, the
 * stem and its partners selected, its chromosome alone.
 */
export const TheStemMenu: StoryObj = { render: () => ring('ring', { opening: stemMenuOpen }) };

/** The partner keys: the cursor on the strongest stem, then `.` twice (`,` back): on to its second partner, said. */
export const ThePartnerKeys: StoryObj = { render: () => ring('ring', { opening: partnerKeys }) };

/** Hovering a stem (the strongest, once placed): its ribbons lit over the washed hole, the card with its counts. */
export const HoveringAStem: StoryObj = { render: () => ring('ring', { opening: hoverStrongestStem }) };

/**
 * Two assemblies on one ring (EC-210, peach and grape): a group gap between them, synteny ribbons from a peach stem to a
 * grape stem, an inverted run crossed.
 */
export const TwoAssembliesOnOneRing: StoryObj = { render: () => ring('ring-two-assemblies') };

/**
 * Crowded (3,000 contigs): the contigs under 2 px a run; stems only on blocks a device pixel wide; links to unread
 * contigs end at their places in the run.
 */
export const Crowded: StoryObj = { name: 'Crowded: 3,000 Contigs', render: () => ring('ring-crowded') };

/** The links pending (their answer never comes): the bands drawn, the hole's pending hatch, no stems. */
export const LinksPending: StoryObj = { render: () => ring('ring-pending') };

/** The links’ file missing: the issue said once, no ribbon, no hatch waiting (“links failed”). */
export const LinksFailed: StoryObj = { render: () => ring('ring-failed', { sourceIssues: true }) };

/** No links at all: an empty hole, its ARIA description "no links". */
export const NoLinks: StoryObj = { render: () => ring('ring-no-links') };

/**
 * Dark: tokens only (the hues `rondo:hues`' derived dark, the cytobands' stain from the text token, the wash and the
 * hover the theme's); the toolbar's theme set dark.
 */
export const Dark: StoryObj = { globals: { theme: 'dark' }, render: () => ring('ring') };
