/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../../_stories/support.js';
import { ScrollArea } from '../scroll-area.js';

const meta: Meta = { title: 'Components/Scroll area' };
export default meta;

const panel =
  'display: grid; grid-template-columns: repeat(auto-fill, 260px); gap: var(--msui-space-3); padding: var(--msui-space-3); background: var(--msui-surface); color: var(--msui-text); font-size: var(--msui-text-sm);';
const caption = 'color: var(--msui-text-muted); font-size: var(--msui-text-xs); margin: 0 0 var(--msui-space-1);';
const row =
  'display: flex; justify-content: space-between; height: 26px; align-items: center; padding: 0 var(--msui-space-2);';
const rows = (n: number) =>
  Array.from({ length: n }, (_, i) => (
    <div style={row}>
      <span>Snapshot {i + 1}</span>
      <span style="color: var(--msui-text-subtle)">#{i + 1}</span>
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
          <div style="display: flex; flex-direction: column; max-height: 200px; border-radius: var(--msui-radius-md); background: var(--msui-surface-sunken); box-shadow: inset 0 0 0 1px var(--msui-line);">
            <ScrollArea aria-label="On a surface">{rows(20)}</ScrollArea>
          </div>
        </section>
        <section>
          <p style={caption}>Across (wide content)</p>
          <ScrollArea axis="x" aria-label="Wide">
            <div style="display: flex; gap: var(--msui-space-1); padding: var(--msui-space-1)">
              {Array.from({ length: 30 }, (_, i) => (
                <span style="flex: none; width: 48px; height: 48px; display: grid; place-items: center; border-radius: var(--msui-radius-md); background: var(--msui-control)">
                  {i + 1}
                </span>
              ))}
            </div>
          </ScrollArea>
        </section>
      </div>
    )),
};
