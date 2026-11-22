/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { BuilderCases } from './builder-compose.js';
import { builderStory } from './builder-view.js';

// The builder (`createSeqSpec`), each story the code beside the view it builds: the code is the source of the module
// the view is built from (Vite's `?raw`), so it cannot drift. Light and dark by the toolbar's theme.

const meta: Meta = { title: 'Features/Builder' };
export default meta;

const story = (id: string): StoryObj => {
  const c = BuilderCases.find((x) => x.id === id)!;
  return { render: () => builderStory(c) };
};

/** The smallest spec that draws: a sequence from letters, one track, one layer. */
export const FirstView: StoryObj = story('first-view');

/**
 * Every call returns the address of what it adds, typed by its kind: the sequence's handle is the ranges' target, the
 * ranges' handle the layer's data. The last line is a type error, kept under `@ts-expect-error` so the typecheck proves it.
 */
export const HandlesAndReferences: StoryObj = story('handles');

/** `theme(...)` and `measure(...)`: a scheme, a categorical field, rules over a value. */
export const LayersAndThemes: StoryObj = story('themes');

/**
 * An assembly with the genome builder, a genome section of it with a link between two chromosomes (`a_name`,
 * `b_name`), and a section of one chromosome of it.
 */
export const AGenomeView: StoryObj = story('genome');

/** A spec with a mistake the types cannot catch: `check()`'s issues are listed beside the view; the SVG export draws what resolves, as the canvas does, and says the issues in its `<desc>`. */
export const Check: StoryObj = story('check');

/** `build()` as data: the normalized tree printed as JSON and as the S-expression, beside the view; the SVG export draws what resolves, as the canvas does, and says the issues in its `<desc>`. */
export const BuildAsData: StoryObj = story('as-data');

/** `views.tooltip` and `view.tooltip`: words on rendered marks, by data and by range, Markdown in the template, an in-view link. Rest on a domain or residue 175. */
export const Tooltips: StoryObj = story('tooltips');
