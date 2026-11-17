/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../../_stories/support.js';
import { State } from '@molstar/core/reactive';
import { type Child, Show } from '../../../dom/index.js';
import { Button } from '../../button/button.js';
import { GlassWordmark } from '../glass-wordmark.js';
import { MolstarWordmark } from '../wordmark.js';

const meta: Meta = { title: 'Brand' };
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

function GlassPanel(props: { theme: 'light' | 'dark' }) {
  // Remounted on demand (a new key): it assembles each time it shows.
  const shown = new State(1);
  return (
    <div
      data-msui-theme={props.theme}
      class="msui-surface"
      style="gap: var(--msui-space-4); padding: var(--msui-space-6);"
    >
      <span style={caption}>
        The wordmark as glass at 1×, 2× and 3×: it assembles as it shows; the pieces lean under the pointer; a press
        shatters it and it snaps back. At rest it draws nothing.
      </span>
      <div style="display: flex; gap: var(--msui-space-6); align-items: center; flex-wrap: wrap; min-height: 96px;">
        <Show when={shown}>
          {() => (
            <>
              <GlassWordmark />
              <GlassWordmark scale={2} />
              <GlassWordmark scale={3} />
            </>
          )}
        </Show>
      </div>
      <div>
        <Button size="sm" onClick={() => shown.set(shown.value + 1)}>
          Show again
        </Button>
      </div>
    </div>
  );
}

/** The glass wordmark (a view's heading, the shell's Home): on a light and a dark theme. */
export const Glass: StoryObj = {
  render: () =>
    story(() => (
      <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: var(--msui-space-5); padding: var(--msui-space-5)">
        <GlassPanel theme="light" />
        <GlassPanel theme="dark" />
      </div>
    )),
};
