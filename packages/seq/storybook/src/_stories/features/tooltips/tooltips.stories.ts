/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';
import { TooltipViews } from './tooltip-views.js';

// The spec's tooltips: words attached to data, shown on the marks that draw it by the card's `seq:tooltips` provider
// (the application composes the card). Each story the canvas and the SVG of one view, light and dark by the toolbar's
// theme; the headless twin is `_test/tooltips-stories.test.ts`.

const meta: Meta = { title: 'Features/Tooltips' };
export default meta;

const story = (id: string): StoryObj => {
  const view = TooltipViews.find((v) => v.id === id)!;
  return { render: () => exampleStory(view, { summary: view.hint, stacked: id === 'svg-titles' }) };
};

/** By data: TP53’s transcripts from Ensembl’s GFF3 speak their own attributes (a tooltip `on` the data node). */
export const ByData: StoryObj = story('by-data');

/** By range: a curator’s notes over p53, inline in the tooltip, on every mark they cover. */
export const ByRangeInline: StoryObj = story('by-range-inline');

/** By range from a source: the same notes as JSON records, read through the store beside the marks they decorate. */
export const ByRangeSource: StoryObj = story('by-range-source');

/** Over an assembly: notes whose rows name their chromosome, on the features they cover. */
export const OverAssembly: StoryObj = story('over-assembly');

/** Markdown in the template with an in-view link (`!seq-focus`): the values stay text. */
export const MarkdownLink: StoryObj = story('markdown-link');

/** The SVG’s titles: `export.titles: card`, each exported item titled by its card (the browser’s tooltip). */
export const SvgTitles: StoryObj = story('svg-titles');
