/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FeatureViews, GlobinsViews, P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { viewerItem } from '@molstar/seqx-viewer/config';
import { exampleStory } from '../../examples/example-root.js';
import { CoverageView, openStrip } from './compose.js';
import { hostNavigatorStory } from './host-navigator-root.js';

// Navigators as such (the navigators design): the ruler rail every linear section has by default, and a declared strip
// that replaces it; on p53, whose tracks show both best; a strip with a follower, and the ideogram of a chromosome over
// a gene's region; the vertical forms over the Globins' rows (the rows rail, the vertical overview, both); a host's
// own picture (a coverage strip), and a navigator the host draws itself in DOM. The other datasets' strips are drawn in
// their own examples.

const meta: Meta = { title: 'Features/Navigators' };
export default meta;

const view = (id: string) => P53Views.find((v) => v.id === id)!;

/** The rows rail replaced by a declared overview (the rule before `navigator.rowsRail`). */
const OverviewAlone = [[viewerItem('navigator.rowsRail'), 'unlessOverview']] as const;

/**
 * The ruler rail, every linear section's navigator by default: the whole axis on the ruler's bottom border, the window a
 * segment with a grip each end. Drag on the ruler to pan, drag a grip to zoom, click beside the window to centre there,
 * drag beside it to sweep a region, double-click to fit (fitted, to zoom back in). It grows while used, and is a quiet
 * hairline when the view shows everything at its zoom limit. An export leaves it out.
 */
export const Rail: StoryObj = { render: () => exampleStory(view('tracks')) };

/**
 * A declared strip over the same tracks replaces the rail: its picture the tracks small, its window the view's. Its
 * window is a slider: focus it and use the arrows, Shift+arrows, Home and End; its value is said when the view rests.
 * The SVG includes it, its window an outline titled with its positions. It opens zoomed in on the DNA-binding domain, so
 * the window is smaller than the strip and can be dragged.
 */
export const Strip: StoryObj = {
  render: () =>
    exampleStory(view('tracks-strip'), {
      navigators: 'include',
      opening: openStrip,
    }),
};

/**
 * The protein on its gene with a strip of the protein: its window and its gene follower's, carried along the codon map,
 * drawn together; a drag on the gene's window moves the protein, which leads.
 */
export const ProteinStripWithItsGene: StoryObj = {
  render: () => exampleStory(view('protein-gene-strip'), { navigators: 'include' }),
};

/**
 * The gene's region with chromosome 17's ideogram above it (GRCh38's bands, the centromere two triangles): the region,
 * a sliver of the chromosome, marked on it, and a zoom strip under it over the region alone, joined to the mark by two
 * lines. Hover a band for its name, the mark for the region; drag the zoom strip's window, its grips, click outside it,
 * double-click. The SVG beside it includes the ideogram and its zoom strip, the window an outline titled with its
 * positions.
 */
export const GeneRegionOnChr17: StoryObj = {
  render: () => exampleStory(view('gene-region-ideogram'), { navigators: 'include' }),
};

/**
 * The rows rail, an alignment's rows' navigator when none is declared: the rail upright at the right edge of the rows'
 * band, beside their numbers, its window the rows in view. Drag it to scroll the rows, a grip to zoom them, click
 * beside the window to centre there; `]` and `[` zoom them. Idle while every row shows.
 */
export const RowsRail: StoryObj = {
  render: () => exampleStory(FeatureViews.find((v) => v.id === 'rows-rail')!),
};

/**
 * The vertical overview of the same rows, declared, alone (`navigator.rowsRail: unlessOverview`, the rail replaced):
 * every row summarized down a column at the right, its window the rows and the columns in view, a grip on each of its
 * four sides. Drag the window to move both, its top and bottom grips to zoom the rows, its left and right grips to zoom
 * the columns; click outside it to centre there. The SVG includes it.
 */
export const VerticalOverview: StoryObj = {
  render: () =>
    exampleStory(GlobinsViews.find((v) => v.id === 'family')!, {
      navigators: 'include',
      config: OverviewAlone,
    }),
};

/**
 * Both together, the default (`navigator.rowsRail: always`): the rows rail at the band's right edge beside the numbers,
 * the precise control, and the overview beside the canvas, the picture. They draw one window: a drag on either moves
 * the other in the same frame.
 */
export const RowsRailAndOverview: StoryObj = {
  render: () => exampleStory(GlobinsViews.find((v) => v.id === 'family')!, { navigators: 'include' }),
};

/**
 * A host's picture: `host:coverage`, registered with the viewer (a `navigator-picture` entry) and drawn as spec tracks
 * by the navigator's own renderer: p53's conservation filling the strip, colored as the view's. The window, its grips,
 * gestures and keys are the navigator's, whatever the picture; hover it for the picture's readout (the residue and its
 * conservation). A reader without the picture draws the axis and says so. The SVG includes it.
 */
export const CoverageStrip: StoryObj = {
  render: () => exampleStory(CoverageView, { navigators: 'include' }),
};

/**
 * A navigator the host draws itself in DOM elements, no canvas: p53's domains and the window placed from the viewer's
 * `NavigatorModel` (`viewer.navigator(def)`), which holds every rule. Drag the window to pan, an end to zoom, click
 * beside it to centre there, double-click to fit, the wheel to zoom; focus it and use the arrows, Shift+arrows, Home and
 * End. Each of the host's events is one of the model's methods.
 */
export const HostDrawn: StoryObj = {
  render: () =>
    hostNavigatorStory({
      hint: 'A navigator the host draws in DOM elements: its gestures, keys, readout and slider value are the viewer’s NavigatorModel’s.',
    }),
};
