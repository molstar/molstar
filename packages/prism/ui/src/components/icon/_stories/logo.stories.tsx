/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../../testing/story-root.js';
import { Icon } from '../icon.js';
import { MolstarWordmark } from '../wordmark.js';

const meta: Meta = { title: 'Components/Logo' };
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

function LogoPanel(props: { theme: 'light' | 'dark' }) {
  return (
    <div
      data-msp-theme={props.theme}
      style="display: grid; gap: 20px; padding: 20px; border-radius: var(--msp-radius-lg); background: var(--msp-color-panel); color: var(--msp-color-text);"
    >
      <div style="display: grid; gap: var(--msp-space-3)">
        <span style={caption}>
          The wordmark (Home's heading; the empty canvas at 3×): hover it for the logo's colors
        </span>
        <div style="display: flex; gap: var(--msp-space-6); align-items: end; flex-wrap: wrap">
          <MolstarWordmark />
          <MolstarWordmark scale={2} />
          <MolstarWordmark scale={3} live />
        </div>
      </div>
      <div style="display: grid; gap: var(--msp-space-3)">
        <span style={caption}>The Home icon (16, 20, 24, 48 px)</span>
        <div style="display: flex; gap: var(--msp-space-5); align-items: center">
          {[16, 20, 24, 48].map((size) => (
            <Icon name="prism" size={size} />
          ))}
        </div>
      </div>
    </div>
  );
}

/**
 * The logo: the Mol* wordmark with Prism's mark after it, a triangle as tall as the asterisk pointing away from it (a
 * prism, a play button) with a thicker core inside (a tetrahedron from above); Home's icon is the mark as a delta.
 */
export const Logo: StoryObj = {
  render: () =>
    story(() => (
      <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: var(--msp-space-5); padding: var(--msp-space-5)">
        <LogoPanel theme="light" />
        <LogoPanel theme="dark" />
      </div>
    )),
};

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
