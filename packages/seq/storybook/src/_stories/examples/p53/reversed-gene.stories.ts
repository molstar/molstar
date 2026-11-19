/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = P53Views.find((v) => v.id === 'reversed-gene')!;

const meta: Meta = { title: 'Examples/p53/Reversed gene' };
export default meta;

/**
 * The gene alone, reversed: positions count down, its first label says so, the letters are the other strand's (the
 * complement). Shift+C collapses the introns to hairlines; Shift+R turns it back.
 */
export const ReversedGene: StoryObj = { render: () => exampleStory(view) };
