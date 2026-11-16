/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../../_stories/support.js';
import type { Child } from '../../../dom/index.js';
import { MolstarWordmark } from '../wordmark.js';

const meta: Meta = { title: 'Components/Logo' };
export default meta;

const surface =
  'display: flex; gap: var(--msui-space-6); align-items: center; flex-wrap: wrap; padding: var(--msui-space-6); background: var(--msui-bg); color: var(--msui-text);';
const item = 'display: grid; gap: var(--msui-space-3); justify-items: start;';
const caption = 'color: var(--msui-text-muted);';

function Labelled(props: { label: string; children: Child }) {
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
      data-msui-theme={props.theme}
      class="msui-surface"
      style="gap: var(--msui-space-3); padding: var(--msui-space-5);"
    >
      <span style={caption}>The wordmark at 1×, 2× and 3× (live): hover it for the logo's colors</span>
      <div style="display: flex; gap: var(--msui-space-6); align-items: end; flex-wrap: wrap">
        <MolstarWordmark />
        <MolstarWordmark scale={2} />
        <MolstarWordmark scale={3} live />
      </div>
    </div>
  );
}

/** The Mol* wordmark, outlined in the text color, on a light and a dark theme. */
export const Logo: StoryObj = {
  render: () =>
    story(() => (
      <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: var(--msui-space-5); padding: var(--msui-space-5)">
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

/** Alive under the pointer (an empty canvas): move across it, on the background and on a dark theme. */
export const Live: StoryObj = {
  render: () =>
    story(() => (
      <div>
        <div style={surface}>
          <Labelled label="on the background">
            <MolstarWordmark scale={3} flash="reveal" live />
          </Labelled>
        </div>
        <div style={surface} data-msui-theme="dark">
          <Labelled label="on a dark theme">
            <MolstarWordmark scale={3} flash="reveal" live />
          </Labelled>
        </div>
      </div>
    )),
};
