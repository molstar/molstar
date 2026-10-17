/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../testing/story-root.js';

const meta: Meta = { title: 'Internals/Performance/Input rate' };
export default meta;

// How many input events a frame brings (rendering.md §2.3, the frame lane): scroll the list with a trackpad or a
// wheel, move the pointer over the box, drag in it. Per kind: events a frame (mean, the most), frames with more than
// one, over the last two seconds of each. A trackpad on macOS sends wheel and pointer events at 120 Hz or more, several
// a frame: input that sets state on every event would render several times a frame; a frame state publishes once.

type Kind = 'wheel' | 'pointermove' | 'scroll';

export const InputRate: StoryObj = {
  name: 'Input rate',
  parameters: { a11y: { disable: true, test: 'off' } },
  render: () =>
    story(() => {
      const counts: Record<Kind, number> = { wheel: 0, pointermove: 0, scroll: 0 };
      // Per kind: events in each frame (the last 2 s of frames with any), the time it was first seen.
      const frames: Record<Kind, number[]> = { wheel: [], pointermove: [], scroll: [] };
      const readout = document.createElement('pre');
      readout.style.cssText = 'margin: 0; font: 12px/1.5 var(--msp-font-mono); color: var(--msp-color-text);';
      let on = true;
      const tick = () => {
        if (!on) return;
        for (const k of Object.keys(counts) as Kind[]) {
          if (counts[k] > 0) {
            frames[k].push(counts[k]);
            if (frames[k].length > 240) frames[k].shift();
          }
          counts[k] = 0;
        }
        const lines = (Object.keys(frames) as Kind[]).map((k) => {
          const f = frames[k];
          if (!f.length) return `${k.padEnd(12)} no events yet`;
          const mean = f.reduce((a, b) => a + b, 0) / f.length;
          const max = Math.max(...f);
          const several = f.filter((n) => n > 1).length;
          return `${k.padEnd(12)} ${mean.toFixed(2)} a frame (most ${max}), ${several} of ${f.length} frames with more than one`;
        });
        readout.textContent = lines.join('\n');
        requestAnimationFrame(tick);
      };
      requestAnimationFrame(tick);
      const list = document.createElement('div');
      list.style.cssText =
        'height: 260px; width: 280px; overflow: auto; border: 1px solid var(--msp-color-border, #8884); font-size: 12px;';
      for (let i = 0; i < 2000; i++) {
        const row = document.createElement('div');
        row.textContent = `Row ${i}`;
        row.style.cssText = 'height: 22px; padding: 0 8px;';
        list.append(row);
      }
      list.addEventListener('wheel', () => counts.wheel++, { passive: true });
      list.addEventListener('scroll', () => counts.scroll++, { passive: true });
      const box = document.createElement('div');
      box.style.cssText =
        'height: 260px; width: 280px; border: 1px solid var(--msp-color-border, #8884); display: grid; place-items: center; color: var(--msp-color-text-subtle); touch-action: none;';
      box.textContent = 'Move or drag the pointer here';
      box.addEventListener('pointermove', () => counts.pointermove++);
      const root = (
        <div style="display: grid; gap: var(--msp-space-4); padding: var(--msp-space-5); color: var(--msp-color-text);">
          <p style="margin: 0; max-width: 60ch;">
            Scroll the list with the trackpad, move the pointer over the box. Each line: events of that kind a frame
            (the last 240 frames that had any).
          </p>
          {readout}
          <div style="display: flex; gap: var(--msp-space-5);">
            {list}
            {box}
          </div>
        </div>
      ) as HTMLElement;
      (window as unknown as { __mspInputRate?: () => unknown }).__mspInputRate = () => ({ ...frames });
      queueMicrotask(() => {
        const stop = new MutationObserver(() => {
          if (!root.isConnected) {
            on = false;
            stop.disconnect();
          }
        });
        stop.observe(document.body, { childList: true, subtree: true });
      });
      return root;
    }),
};
