/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The track menu's "Focus" in both default UIs on real surfaces and the real clock (rAF, timers): the menu opens, the
// item is chosen, the menu closes and gives its trigger focus back, and the section's viewport animates onto what the
// track covers. Pan the view away first (arrow keys on the surface, as a person would), then Focus.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import { afterEach, describe, expect, it } from 'vitest';
import { userEvent } from 'vitest/browser';
import { windowScheduler } from '../../src/canvas/index.js';
import * as D from '../../src/dom/index.js';
import '../../src/skin.css';
import { createSeqViewer, DefaultSeqViewerSpec, SeqViewerConfig, type SeqViewerContext } from '../../src/index.js';
import { SeqViewerLayout, SeqViewerProvider } from '../../src/react/index.js';
import { at } from '../dom/fixtures.js';

(globalThis as { IS_REACT_ACT_ENVIRONMENT?: boolean }).IS_REACT_ACT_ENVIRONMENT = true;

const cleanups: (() => void | Promise<void>)[] = [];
afterEach(async () => {
  for (const c of cleanups.splice(0).reverse()) await c();
});

const SEQ = Array.from({ length: 360 }, (_, i) => 'ACDEFGHIKLMNPQRSTVWY'[(i * 7) % 20]).join('');
/**
 * A 360-residue sequence and one feature set with a Domain at [0, 360) and two Motifs near its middle, shown as two
 * tracks that each filter the set by category (as UniProt-style views do): Focus on "Motif" zooms onto [150, 210), not
 * onto the whole set.
 */
function doc() {
  const b = createSeqViewSpec({ title: 'Focus' });
  const s = b.data.sequence({ ref: 's', value: SEQ });
  const a = b.data.annotation({
    ref: 'a',
    geometry: 'ranges',
    target: s,
    value: [
      { start: 0, stop: 360, label: 'Whole', category: 'Domain' },
      { start: 150, stop: 170, label: 'One', category: 'Motif' },
      { start: 190, stop: 210, label: 'Two', category: 'Motif' },
    ],
  });
  const v = b.view({ ref: 'v', label: 'V' }).block({ sequence: s });
  v.track({ ref: 'v.ruler', label: '' }).ruler();
  v.track({ ref: 'v.letters', label: 'Residues' }).letters();
  v.track({ ref: 'v.domain', label: 'Domain' }).ranges({ data: a, filter: { field: 'category', values: ['Domain'] } });
  v.track({ ref: 'v.motif', label: 'Motif' }).ranges({ data: a, filter: { field: 'category', values: ['Motif'] } });
  return b.getDocument();
}

const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
const frames = async (n: number) => {
  for (let i = 0; i < n; i++) await new Promise((r) => requestAnimationFrame(r));
};

async function open(ui: 'react' | 'dom'): Promise<{ ctx: SeqViewerContext; container: HTMLElement }> {
  // `full` chrome: the track menu the test opens (the default is clean).
  const spec = DefaultSeqViewerSpec();
  const ctx = createSeqViewer(
    { ...spec, config: [...(spec.config ?? []), [SeqViewerConfig.UI.Chrome, 'full']] },
    { scheduler: windowScheduler() },
  );
  const r = await ctx.setDocument(doc());
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  const container = document.createElement('div');
  container.style.cssText = "width:900px;--seq-font-ui:Arial;--seq-font-mono:'Courier New'";
  document.body.appendChild(container);
  cleanups.push(
    () => ctx.dispose(),
    () => container.remove(),
  );
  if (ui === 'react') {
    const root = createRoot(container);
    await act(async () =>
      root.render(
        <SeqViewerProvider context={ctx} ui={{}}>
          <SeqViewerLayout />
        </SeqViewerProvider>,
      ),
    );
    cleanups.push(() => act(async () => root.unmount()));
  } else {
    const m = D.mountSeqViewer(container, ctx, {});
    cleanups.push(() => m.dispose());
  }
  const sec = ctx.viewports.sections.value[0];
  for (let k = 0; k < 100; k++) {
    await frames(1);
    const w = ctx.layout.get(sec)?.current.width ?? 0;
    if (w > 100) break;
  }
  await frames(3);
  return { ctx, container };
}

describe.each(['react', 'dom'] as const)("the track menu's Focus (%s)", (ui) => {
  it('zooms the viewport onto what the track covers, after arrow keys moved the view', async () => {
    const { ctx, container } = await open(ui);
    const sec = ctx.viewports.sections.value[0];
    const vp = ctx.viewports.get(sec);
    if (!vp) throw new Error('no viewport');
    const issues: string[] = [];
    ctx.events.issue.subscribe((i) => issues.push(`${i.code} ${i.message}`));
    const surface = container.querySelector('.seqv-surface [tabindex]') as HTMLElement;
    // The person works on the surface first: it holds the keyboard focus, the cursor is on a track.
    surface.focus();
    for (const key of ['ArrowDown', 'ArrowRight', 'ArrowRight'])
      surface.dispatchEvent(new KeyboardEvent('keydown', { key, bubbles: true, cancelable: true }));
    await frames(3);
    // The pointer over the Features header shows its actions box; its trigger opens the menu.
    const h = at(ctx, sec, 'v.motif');
    await userEvent.hover(surface, { position: { x: h.header, y: h.y } });
    await frames(3);
    const trigger = container.querySelector('.seqv-track-actions button') as HTMLElement;
    expect(trigger).toBeTruthy();
    const before = vp.live.current.cell;
    await userEvent.click(trigger);
    await frames(2);
    const item = [...container.querySelectorAll<HTMLElement>('[role^=menuitem]')].find(
      (e) => e.textContent === 'Focus',
    );
    expect(item).toBeTruthy();
    await userEvent.click(item as HTMLElement);
    for (let k = 0; k < 100 && (k < 3 || vp.live.current.animating); k++) await wait(50);
    await frames(2);
    const s = vp.live.current;
    const f = ctx.focus.current.value.get(sec);
    // The Motif items cover [150, 210) of 360, not the Domain's whole sequence the set also holds: the view fits them.
    expect({ focus: f && [f.block, f.start, f.stop], issues }).toEqual({ focus: [0, 150, 210], issues: [] });
    const w = s.width - (ctx.layout.get(sec)?.current.headerWidth ?? 0);
    expect(s.cell).toBeGreaterThan(before * 4);
    // The 60 residues fill most of the data area (a tenth of padding each side).
    expect(s.cell * 60).toBeGreaterThan(w * 0.7);
    expect(s.cell * 60).toBeLessThan(w);
  });
});
