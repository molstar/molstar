/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BuiltInSchemes } from '@molstar/seqx-core/theme/index';
import { FeatureViews, P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';
import { SchemesGalleryRoot, ThemesEditingRoot } from './themes-roots.js';

// p53's themes: the view with every scheme that fits its letters, its legends and their export; the view edited (a scheme
// previewed while hovered, picked in one step; the shared scale's domain set for both signals in one step; the reader's
// color-blind substitution); and the gallery of the built-in schemes in the light and dark tables, each seen through the
// three simulated color-vision deficiencies. Strings from the spec go in as text, never markup.

const view = P53Views.find((v) => v.id === 'themes')!;
const editingView = FeatureViews.find((v) => v.id === 'themes-editing')!;

const meta: Meta = { title: 'Features/Themes' };
export default meta;

/** p53 under every theme form: schemes on the letters, categories with their legend, pLDDT bands and a ramp, a shared scale, a log-sized variant track; the export draws the legends below. */
export const Themes: StoryObj = { render: () => exampleStory(view) };

/** The themes edited: a scheme previewed on hover and picked in one step, the shared scale set for both signals, the reader's substitution. */
export const EditingThemes: StoryObj = {
  render: () => {
    const el = document.createElement(ThemesEditingRoot.tag) as ThemesEditingRoot;
    el.model = editingView;
    return el;
  },
};

/** The built-in schemes in both tables, each under the three simulated color-vision deficiencies. */
export const SchemesGallery: StoryObj = {
  render: () => {
    const el = document.createElement(SchemesGalleryRoot.tag) as SchemesGalleryRoot;
    el.model = BuiltInSchemes;
    return el;
  },
};
