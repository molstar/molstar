// @vitest-environment jsdom
/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { viewerItem } from '../config.js';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { legendRows } from '../legends.js';
import { mountLegend } from '../mount-legend.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { stubCanvasContext } from './canvas-stub.js';
import { type Frames, installFrames } from './frames.js';

let frames: Frames | undefined;
let viewer: SeqViewer | undefined;
beforeEach(stubCanvasContext);
afterEach(() => {
  viewer?.retire();
  frames?.restore();
  viewer = undefined;
  frames = undefined;
});

const Chromosomes = '/views/main/p/chroms';

function spec(n: number) {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'dna', letters: 'A'.repeat(n * 10) });
  const d = b.data.ranges('d', {
    target: p,
    rows: {
      start: Array.from({ length: n }, (_, i) => i * 10),
      end: Array.from({ length: n }, (_, i) => i * 10 + 8),
      id: Array.from({ length: n }, (_, i) => `i${i}`),
      category: Array.from({ length: n }, (_, i) => `chr${i + 1}`),
    },
  });
  b.views.view('main', (v) => {
    const s = v.section('p', { of: p });
    s.track('letters').layer('l', { type: 'letters' });
    s.track('chroms', { label: 'Chromosomes' }).layer('d', {
      type: 'ranges',
      data: d,
      color: theme.field('category', { list: false }),
    });
  });
  return b.build();
}

async function rest() {
  for (let k = 0; k < 3; k++) {
    await viewer!.settle();
    frames!.timers.fire();
    frames!.flush();
  }
}

async function load(n: number) {
  frames = installFrames();
  viewer = createSeqViewer({
    ...DefaultSeqViewerSpec(),
    ui: undefined,
    config: [
      [viewerItem('viewport.maxPxPerPosition'), 1000],
      [viewerItem('view.motion'), 'cut'],
    ],
  });
  await viewer.load(spec(n));
  viewer.resize(n > 1000 ? 30000 : 500);
  await viewer.settle();
  frames.flush();
}

const more = (el: HTMLElement) => el.querySelector<HTMLButtonElement>('button[data-legend-more]');
const labels = (el: HTMLElement) => [...el.querySelectorAll('[role="option"]')].map((o) => o.textContent);

async function mounted(n: number, options: Parameters<typeof mountLegend>[2] = {}) {
  await load(n);
  await rest();
  const el = document.createElement('div');
  document.body.append(el);
  return { el, dispose: mountLegend(el, viewer!, { of: Chromosomes, ...options }) };
}

describe('a legend’s "+n"', () => {
  it('is a button with a label and aria-expanded; open lists every entry, close returns to 12 and "+n"', async () => {
    const { el } = await mounted(30);
    expect(labels(el)).toHaveLength(12);
    const button = more(el)!;
    expect(button.tagName).toBe('BUTTON');
    expect(button.textContent).toBe('+18');
    expect(button.getAttribute('aria-expanded')).toBe('false');
    expect(button.getAttribute('aria-label')).toBe('+18, show 18 more: Chromosomes');
    button.click();
    expect(labels(el)).toHaveLength(30);
    expect(labels(el)[29]).toBe('chr30');
    expect(more(el)!.getAttribute('aria-expanded')).toBe('true');
    expect(more(el)!.getAttribute('aria-label')).toBe('Show fewer: Chromosomes');
    more(el)!.click();
    expect(labels(el)).toHaveLength(12);
    expect(more(el)!.textContent).toBe('+18');
    expect(more(el)!.getAttribute('aria-expanded')).toBe('false');
  });

  it('keeps the focus on the button across the redraw, and Enter does not select an entry', async () => {
    const { el } = await mounted(30);
    const first = more(el)!;
    first.focus();
    first.dispatchEvent(new KeyboardEvent('keydown', { key: 'Enter', bubbles: true }));
    expect(viewer!.selection.value.value.elements.keys).toEqual([]);
    first.click();
    expect(document.activeElement).toBe(more(el));
    expect(document.activeElement).not.toBe(first);
  });

  it('the state is the session’s: it survives a redraw and a pan that keeps the legend, and is no spec’s', async () => {
    const { el } = await mounted(30);
    more(el)!.click();
    viewer!.legends.setCvdSafe(true);
    expect(labels(el)).toHaveLength(30);
    // a pan to 20 of the 30 items: the legend (its values changed, its key with them) is the same one, still open
    const before = viewer!.legends.legends.value[0].key;
    viewer!.runCommand('view.focus', { section: '/views/main/p', focus: { start: 0, end: 200 } });
    await rest();
    expect(viewer!.legends.legends.value[0].key).not.toBe(before);
    expect(more(el)!.getAttribute('aria-expanded')).toBe('true');
    expect(labels(el).length).toBeGreaterThan(12);
    expect(labels(el).length).toBeLessThan(30);
    // the rows say so without the DOM
    expect(legendRows(viewer!, { of: Chromosomes })[0].open).toBe(true);
    // a second part over the viewer shows the same
    const other = document.createElement('div');
    document.body.append(other);
    mountLegend(other, viewer!, { of: Chromosomes });
    expect(more(other)!.getAttribute('aria-expanded')).toBe('true');
  });

  it('a legend over 500 values opens to a bounded, scrollable list', async () => {
    const { el } = await mounted(500);
    expect(labels(el)).toHaveLength(12);
    expect(more(el)!.textContent).toBe('+488');
    more(el)!.click();
    // the DOM stays bounded: the first 200, the rest counted inside the box
    expect(labels(el)).toHaveLength(200);
    expect(el.querySelector<HTMLElement>('[data-beyond]')!.dataset.beyond).toBe('300');
    const group = el.querySelector<HTMLElement>('[role="group"]')!;
    expect(group.style.maxHeight).not.toBe('');
    expect(group.style.overflowY).toBe('auto');
    more(el)!.click();
    expect(el.querySelector<HTMLElement>('[role="group"]')!.style.maxHeight).toBe('');
  });

  it('builds rows only for what is shown: 12 closed, a bounded number open, over 3,000 values', async () => {
    await mounted(3000);
    const names = vi.spyOn(viewer!.legends, 'names');
    const closed = legendRows(viewer!, { of: Chromosomes })[0];
    expect(closed.entries).toHaveLength(12);
    expect(closed.more).toBe(2988);
    expect(names).toHaveBeenCalledTimes(12);
    names.mockClear();
    viewer!.legends.setOpen(closed.id, true);
    names.mockClear();
    const open = legendRows(viewer!, { of: Chromosomes })[0];
    expect(open.entries).toHaveLength(200);
    expect(open.beyond).toBe(2800);
    expect(names).toHaveBeenCalledTimes(200);
  });
});
