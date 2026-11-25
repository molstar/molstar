/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FeatureViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';

// Track groups (Design/Track groups): p53's protein section in three groups, a story a fold state. Each is a viewer and
// its SVG export side by side, in the light and the dark theme.

const meta: Meta = { title: 'Features/Track groups' };
export default meta;

const view = (id: string) => FeatureViews.find((v) => v.id === id)!;

/**
 * The opening state the spec writes: Features open (a header line, its tracks indented), Variants folded to its
 * authored density track (shown only while folded), Signals folded to a strip a signal. Click a group's header to fold
 * or open it (not undone: a fold is the reader's); the section's menu has Collapse all and Expand all groups (Alt+Left
 * after one brings the folds back).
 */
export const Opening: StoryObj = { render: () => exampleStory(view('groups')) };

/** Every group open: each group's header line above its tracks, `trackGap` within a group and `groupGap` around it. */
export const AllOpen: StoryObj = { render: () => exampleStory(view('groups-open')) };

/**
 * Every group folded to one row: its label column the header with a quiet count, its band the summary (Features' ranges
 * and secondary structure in one lane, Variants' density, Signals' three strips).
 */
export const AllFolded: StoryObj = { render: () => exampleStory(view('groups-folded')) };
