/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Heading, Note, page, story } from '../support.js';

const meta: Meta = { title: 'Foundations/Type' };
export default meta;

const Sizes: readonly [size: string, use: string][] = [
  ['xs', 'Hints, counts, captions'],
  ['sm', 'The UI: rows, labels, values, controls'],
  ['md', 'Headings in a panel; reading text'],
  ['lg', 'Page and dialog titles'],
];

/**
 * The sizes, each with its line height on the 4px grid (11/16, 12/16, 14/20, 16/24): text of different sizes in
 * neighboring columns still lands on the grid, so rows of a list and lines of a paragraph next to it agree.
 */
export const Sizes_: StoryObj = {
  name: 'Sizes',
  render: () =>
    story(() => (
      <div style={page}>
        {Sizes.map(([s, use]) => (
          <div style="display: grid; grid-template-columns: 280px minmax(0, 1fr); column-gap: var(--msui-space-4); align-items: baseline;">
            <code style="font-size: var(--msui-text-xs); color: var(--msui-text-muted);">
              --msui-text-{s} / --msui-leading-{s}
            </code>
            <div style={{ fontSize: `var(--msui-text-${s})`, lineHeight: `var(--msui-leading-${s})` }}>
              The quick brown fox jumps over the lazy dog. <span style="color: var(--msui-text-muted)">{use}</span>
            </div>
          </div>
        ))}
        <Heading>Weights</Heading>
        <Note>Two: regular for everything, strong (600) for headers and the one action that matters.</Note>
        <div style="display: flex; gap: var(--msui-space-5);">
          <span style="font-weight: var(--msui-weight-regular)">Regular 400</span>
          <span style="font-weight: var(--msui-weight-strong)">Strong 600</span>
        </div>
        <Heading>Monospace</Heading>
        <Note>Sequences, code, identifiers.</Note>
        <code style="font-family: var(--msui-font-mono); font-size: var(--msui-text-sm);">
          MALYGTHSHGLFKKLGIPGPTPLP · 1TQN · auth_seq_id
        </code>
        <Heading>Numbers</Heading>
        <Note>
          Tabular figures wherever numbers stack (`.msui-num`): every digit the same width, so places line up.
        </Note>
        <div style="display: grid; grid-template-columns: 120px 120px; column-gap: var(--msui-space-5);">
          <div style="display: grid;">
            <span style="color: var(--msui-text-muted)">Proportional</span>
            {['1,111.11', '88.08', '7.10'].map((n) => (
              <span style="text-align: end; font-variant-numeric: proportional-nums;">{n}</span>
            ))}
          </div>
          <div style="display: grid;">
            <span style="color: var(--msui-text-muted)">Tabular</span>
            {['1,111.11', '88.08', '7.10'].map((n) => (
              <span class="msui-num">{n}</span>
            ))}
          </div>
        </div>
      </div>
    )),
};
