/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Wrapped text scrolls to a focus (seqspec:interaction#focus; the `focusViewport` behavior's `reveal`): the lines of
// the focus come into view, as the tracks' zoom brings the focus into view; nothing moves where they already are.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import { DefaultPresets, KernelRenderers, type LayerRenderer } from '../src/index.js';
import { opened } from './layout-fixtures.js';

function text(n: number): SeqViewSpecDoc {
  const b = createSeqViewSpec();
  const value = 'MSEEKLKELAEKAIRLAKEGNPEEVLRALQ'.repeat(Math.ceil(n / 30)).slice(0, n);
  b.data.sequence({ ref: 'a', alphabet: 'protein', value });
  b.view({ ref: 'v', label: 'V', layout: 'text' }).track({ ref: 'v.letters' }).letters();
  return b.getDocument();
}

const BUILT_IN = {
  presets: false as const,
  renderers: [...KernelRenderers] as LayerRenderer[],
  registry: [{ presets: [...DefaultPresets] }],
};

/** 2,000 positions in a 400 × 120 px strip: about 60 lines, a few visible. */
async function open() {
  const t = await opened(text(2000), {
    ...BUILT_IN,
    size: { 'v.main': [400, 120] },
    render: { wrappedGroupLabels: 'none', wrappedGutter: 0 },
  });
  t.layout('v.main');
  return t;
}
const scrollY = (t: Awaited<ReturnType<typeof open>>) => t.ctx.viewports.get('v.main')?.live.current.scrollY ?? -1;
const focus = (t: Awaited<ReturnType<typeof open>>, start: number, stop: number) => {
  const loci = t.ctx.refs.toLoci('a', { start, stop });
  expect(loci).not.toBeNull();
  return t.ctx.focus.fromLoci(loci as never, { section: 'v.main' });
};

/** The line (top, bottom) holding shown position k, in the scrolled content's px. */
function lineOf(t: Awaited<ReturnType<typeof open>>, k: number) {
  const l = t.layout('v.main');
  for (let n = 0; n < l.lines.count; n++)
    if (l.lines.start[n] <= k && k < l.lines.stop[n]) {
      return { top: l.lines.top[n], bottom: n + 1 < l.lines.count ? l.lines.top[n + 1] : l.scrolledHeight };
    }
  throw new Error(`no line holds ${k}`);
}

describe('wrapped text scrolls to a focus', () => {
  it('a focus on a far residue brings its line into view', async () => {
    const t = await open();
    expect(t.ctx.layout.isWrapped('v.main')).toBe(true);
    expect(scrollY(t)).toBe(0);
    expect(focus(t, 1500, 1510)).toBe(true);
    const view = t.layout('v.main').height - t.layout('v.main').scrolledTop;
    const line = lineOf(t, 1500);
    const y = scrollY(t);
    expect(y).toBeGreaterThan(0);
    expect(line.top).toBeGreaterThanOrEqual(y);
    expect(line.bottom).toBeLessThanOrEqual(y + view + 0.5);
  });

  it('no scroll where the lines already are visible', async () => {
    const t = await open();
    expect(focus(t, 3, 8)).toBe(true);
    expect(scrollY(t)).toBe(0);
    focus(t, 1500, 1510);
    const there = scrollY(t);
    // Again, and a residue on the same line: nothing moves.
    focus(t, 1500, 1510);
    focus(t, 1502, 1504);
    expect(scrollY(t)).toBe(there);
  });

  it('scrolls back up to a focus above the view', async () => {
    const t = await open();
    focus(t, 1500, 1510);
    expect(scrollY(t)).toBeGreaterThan(0);
    focus(t, 0, 5);
    expect(scrollY(t)).toBe(0);
  });

  it('is a behavior param: `reveal: false` leaves the text where it is', async () => {
    const t = await open();
    t.ctx.behaviors.get('focusViewport')?.setParams({ reveal: false });
    focus(t, 1500, 1510);
    expect(scrollY(t)).toBe(0);
  });
});
