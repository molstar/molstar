/** @jsxImportSource @molstar/plugin-prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../_stories/story-root.js';
import { ScrollArea } from './scroll-area.js';

const meta: Meta = { title: 'Components/Scroll area' };
export default meta;

const panel =
  'display: grid; grid-template-columns: repeat(auto-fill, 260px); gap: var(--msp-space-4); padding: var(--msp-space-4); background: var(--msp-color-panel); color: var(--msp-color-text); font-size: var(--msp-font-size-sm);';
const caption =
  'color: var(--msp-color-text-muted); font-size: var(--msp-font-size-xs); margin: 0 0 var(--msp-space-2);';
const row =
  'display: flex; justify-content: space-between; height: 26px; align-items: center; padding: 0 var(--msp-space-3);';
const rows = (n: number) =>
  Array.from({ length: n }, (_, i) => (
    <div style={row}>
      <span>Snapshot {i + 1}</span>
      <span style="color: var(--msp-color-text-subtle)">#{i + 1}</span>
    </div>
  ));

/**
 * Scrolling areas: overlay bars (shown on hover and while scrolling, hidden when nothing overflows; drag the thumb or
 * press the track), and edges that fade where there is more to scroll.
 */
export const Areas: StoryObj = {
  render: () =>
    story(() => (
      <div style={panel}>
        <section>
          <p style={caption}>A long list, at a fixed height</p>
          <div style="display: grid; height: 240px">
            <ScrollArea aria-label="Long list">{rows(40)}</ScrollArea>
          </div>
        </section>
        <section>
          <p style={caption}>Short: nothing overflows, no bar, no fade</p>
          <ScrollArea aria-label="Short list">{rows(4)}</ScrollArea>
        </section>
        <section>
          <p style={caption}>On its own surface (the surface does not fade, its content does)</p>
          <div style="display: flex; flex-direction: column; max-height: 200px; border-radius: var(--msp-radius-md); background: var(--msp-color-well); box-shadow: inset 0 0 0 1px var(--msp-glass-edge);">
            <ScrollArea aria-label="On a surface">{rows(20)}</ScrollArea>
          </div>
        </section>
        <section>
          <p style={caption}>Across (wide content)</p>
          <ScrollArea axis="x" aria-label="Wide">
            <div style="display: flex; gap: var(--msp-space-2); padding: var(--msp-space-2)">
              {Array.from({ length: 30 }, (_, i) => (
                <span style="flex: none; width: 48px; height: 48px; display: grid; place-items: center; border-radius: var(--msp-radius-md); background: var(--msp-color-control)">
                  {i + 1}
                </span>
              ))}
            </div>
          </ScrollArea>
        </section>
      </div>
    )),
};
