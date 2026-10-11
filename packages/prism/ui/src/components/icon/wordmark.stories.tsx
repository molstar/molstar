/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../testing/story-root.js';
import { MolstarWordmark } from './wordmark.js';

const meta: Meta = { title: 'Components/Wordmark' };
export default meta;

const surface =
  'display: flex; gap: var(--msp-space-6); align-items: center; flex-wrap: wrap; padding: var(--msp-space-6); background: var(--msp-color-bg); color: var(--msp-color-text);';
const item = 'display: grid; gap: var(--msp-space-3); justify-items: start;';
const caption = 'color: var(--msp-color-text-muted); font-size: var(--msp-font-size-sm);';

function Labelled(props: { label: string; children: Node }) {
  return (
    <div style={item}>
      {props.children}
      <span style={caption}>{props.label}</span>
    </div>
  );
}

/** The flash of the logo's colors: now and then (rare, a heading), often, or wiped in once and kept (reveal). */
export const Flash: StoryObj = {
  render: () =>
    story(() => (
      <div style={surface}>
        <Labelled label="rare">
          <MolstarWordmark flash="rare" />
        </Labelled>
        <Labelled label="often">
          <MolstarWordmark flash="often" />
        </Labelled>
        <Labelled label="reveal">
          <MolstarWordmark flash="reveal" />
        </Labelled>
      </div>
    )),
};

/** Whole multiples of its drawn size keep the lines on the pixel grid. */
export const Scale: StoryObj = {
  render: () =>
    story(() => (
      <div style={surface}>
        <Labelled label="scale 1">
          <MolstarWordmark scale={1} />
        </Labelled>
        <Labelled label="scale 3">
          <MolstarWordmark scale={3} />
        </Labelled>
      </div>
    )),
};

/** Alive under the pointer (the empty canvas): move across it, on the background and on a dark surface. */
export const Live: StoryObj = {
  render: () =>
    story(() => (
      <div>
        <div style={surface}>
          <Labelled label="on the background">
            <MolstarWordmark scale={3} flash="reveal" live />
          </Labelled>
        </div>
        <div style={surface} data-msp-theme="dark">
          <Labelled label="on a dark surface">
            <MolstarWordmark scale={3} flash="reveal" live />
          </Labelled>
        </div>
      </div>
    )),
};
