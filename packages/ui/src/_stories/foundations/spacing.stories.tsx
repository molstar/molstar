/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Heading, Note, page, story, TokenTable } from '../support.js';

const meta: Meta = { title: 'Foundations/Spacing' };
export default meta;

const Scale = [
  ['--msui-unit', 'The grid: every padding, margin, gap and size is a multiple (hairlines apart)'],
  ['--msui-space-1', 'Inside a row (icon to text); a label to its hint'],
  ['--msui-space-2', "Between controls in a group; a control's text inset"],
  ['--msui-space-3', 'The inset; between groups of a form'],
  ['--msui-space-4', "A card's padding; between sections"],
  ['--msui-space-5', 'Between blocks of a page'],
  ['--msui-space-6', 'Empty states; large separations'],
] as const;

const Sizes = [
  ['--msui-row', 'A row of a list or tree'],
  ['--msui-control-height', 'Fields, selects, buttons: a unit less than a row, centered in it'],
  ['--msui-header', "A surface's header"],
  ['--msui-icon', 'An icon'],
] as const;

/**
 * The spacing scale on one 4px unit (styles §3): every step is on the grid, so anything placed by the scale lands where
 * something else placed by it does.
 */
export const Scale_: StoryObj = {
  name: 'Scale',
  render: () =>
    story(() => (
      <div style={page}>
        <TokenTable rows={Scale} bars />
      </div>
    )),
};

/** The sizes: rows, controls and headers are grid multiples, and a control sits in a row (styles §4). */
export const Sizes_: StoryObj = {
  name: 'Sizes',
  render: () =>
    story(() => (
      <div style={page}>
        <TokenTable rows={Sizes} bars />
        <Heading>A list beside a form</Heading>
        <Note>A form row is a list row (a control sits in it): the two columns line up row by row.</Note>
        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: var(--msui-space-4); max-width: 640px;">
          <div class="msui-surface">
            {['Polymer', 'Ligand', 'Water', 'Ions'].map((l) => (
              <div class="msui-row msui-row-plain">
                <span>{l}</span>
                <span class="msui-num" style="color: var(--msui-text-subtle)">
                  {l.length * 37}
                </span>
              </div>
            ))}
          </div>
          <div class="msui-surface">
            <div class="msui-form">
              {['Color', 'Size', 'Opacity', 'Quality'].map((l) => [
                <span class="msui-label">{l}</span>,
                <span class="msui-control">{l === 'Opacity' ? '0.8' : 'Auto'}</span>,
              ])}
            </div>
          </div>
        </div>
      </div>
    )),
};
