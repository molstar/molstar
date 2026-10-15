/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover follows the content under a resting pointer (seqspec:interaction#input-events), in both default UIs on real
// surfaces and the manual clock: the pointer moves once over the rows, then a vertical wheel, a horizontal wheel, a
// Ctrl-wheel zoom, a programmatic zoom and scroll move the view under it; at the next frame the hover is what is under
// the pointer and the tooltip says so. A pointer that left gets nothing; a `follow`er moved by its leader updates under
// the pointer.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import { afterEach, describe, expect, it } from 'vitest';
import * as D from '../../src/dom/index.js';
import '../../src/skin.css';
import {
  createSeqViewer,
  DefaultSeqViewerSpec,
  Hit,
  HoverKind,
  type ManualScheduler,
  manualScheduler,
  type SeqViewerContext,
} from '../../src/index.js';
import { SeqViewerLayout, SeqViewerProvider } from '../../src/react/index.js';
import { fixtureMeasurer } from '../../src/testing/index.js';

(globalThis as { IS_REACT_ACT_ENVIRONMENT?: boolean }).IS_REACT_ACT_ENVIRONMENT = true;

const cleanups: (() => void | Promise<void>)[] = [];
afterEach(async () => {
  for (const c of cleanups.splice(0).reverse()) await c();
});

/** 80 rows × 120 columns by value (`ROW<k>`), in a box of 8 rows; with `follow`, a second section following it. */
function doc(follow: boolean) {
  const AA = 'ACDEFGHIKLMNPQRSTVWY';
  const lines = ['# STOCKHOLM 1.0'];
  for (let k = 0; k < 80; k++) {
    let r = '';
    for (let c = 0; c < 120; c++) r += (c + k) % 9 === 0 ? '-' : AA[(c * 7 + k * 3) % 20];
    lines.push(`ROW${k} ${r}`);
  }
  lines.push('//', '');
  const b = createSeqViewSpec({ title: 'Resting' });
  const aln = b.data.alignment({ ref: 'aln', source: { value: lines.join('\n'), format: 'stockholm' } });
  const v = b.view({ ref: 'v', label: 'V', style: 'msa' });
  const lead = v.section({ ref: 'v.lead', height: 'box', height_rows: 8 } as never).block({ alignment: aln });
  lead.track({ label: '' }).ruler();
  lead.track({ label: 'Rows' }).rows();
  if (follow) {
    const f = v.section({ ref: 'v.follow', follow: 'v.lead', height: 'box', height_rows: 8 } as never);
    f.block({ alignment: aln }).track({ label: 'Rows' }).rows();
  }
  return b.getDocument();
}

interface Opened {
  ctx: SeqViewerContext;
  container: HTMLElement;
  scheduler: ManualScheduler;
  surface(ref: string): HTMLElement;
  /** `input` (an event, a move), then the next frame, in React's `act`. */
  frame(input?: () => void): Promise<void>;
}

/** The default UI (`react` or `dom`) on real surfaces, the manual clock; resolved once every surface is sized. */
async function open(ui: 'react' | 'dom', follow = false): Promise<Opened> {
  const scheduler = manualScheduler();
  const ctx = createSeqViewer(DefaultSeqViewerSpec(), { scheduler, text: fixtureMeasurer() });
  const r = await ctx.setDocument(doc(follow));
  if (!r.ok) throw new Error(JSON.stringify(r.issues));
  const container = document.createElement('div');
  container.style.cssText = "width:900px;--seq-font-ui:Arial;--seq-font-mono:'Courier New'";
  document.body.appendChild(container);
  cleanups.push(() => ctx.dispose());
  cleanups.push(() => container.remove());
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
  const frame = (input?: () => void) =>
    act(async () => {
      input?.();
      scheduler.advance(16);
      await new Promise((res) => setTimeout(res, 0));
    });
  const refs = ctx.viewports.sections.value;
  const surface = (ref: string) => {
    const sections = [...container.querySelectorAll('.seqv-surface')];
    return sections[refs.indexOf(ref)]?.querySelector('[tabindex]') as HTMLElement;
  };
  // Settled: every layout as wide as its surface for three browser frames (the ResizeObserver has run).
  const sized = () =>
    refs.every((ref) => {
      const w = ctx.layout.get(ref)?.current.width ?? 0;
      return w > 0 && Math.abs(w - (surface(ref)?.clientWidth ?? -1)) < 1;
    });
  for (let k = 0, stable = 0; k < 200 && stable < 3; k++) {
    await frame();
    await new Promise((res) => requestAnimationFrame(res));
    stable = sized() ? stable + 1 : 0;
  }
  await frame();
  return { ctx, container, scheduler, surface, frame };
}

/** The point (x, y in the surface's px) as the browser sends it: a primary pointer's move, a wheel. */
function moveTo(el: HTMLElement, x: number, y: number) {
  const r = el.getBoundingClientRect();
  el.dispatchEvent(
    new PointerEvent('pointermove', {
      clientX: r.left + x,
      clientY: r.top + y,
      isPrimary: true,
      pointerId: 1,
      pointerType: 'mouse',
      bubbles: true,
    }),
  );
}
function wheel(el: HTMLElement, x: number, y: number, deltaX: number, deltaY: number, ctrlKey = false) {
  const r = el.getBoundingClientRect();
  const e = new WheelEvent('wheel', {
    clientX: r.left + x,
    clientY: r.top + y,
    deltaX,
    deltaY,
    ctrlKey,
    bubbles: true,
    cancelable: true,
  });
  el.dispatchEvent(e);
  return e.defaultPrevented;
}

/** The hovered cell (row, column) and what is under the pointer now, by the layout's own hit test. */
function cells(t: Opened, ref: string) {
  const h = t.ctx.hover.live.current;
  const p = t.ctx.pointer.current;
  const g = Hit.create();
  t.ctx.layout.hitTest(ref, p.x, p.y, g);
  return {
    hover: h.kind === HoverKind.aligned ? [h.aligned.row, h.aligned.column] : [h.kind],
    under: g.kind === HoverKind.aligned ? [g.aligned.row, g.aligned.column] : [g.kind],
  };
}
const tip = (t: Opened) => t.container.querySelector('.seqv-tooltip')?.textContent ?? '';

/** The pointer over the rows' third row, at 40% of the data area. */
async function rest(t: Opened, ref: string) {
  const l = t.ctx.layout.get(ref)?.current;
  const ax = t.ctx.layout.axis(ref);
  if (!l || !ax) throw new Error('no layout');
  const x = ax.left + 0.4 * (l.width - ax.left) + 0.25;
  const y = l.scrolledTop + 2.5 * t.ctx.layout.props(ref).rowHeight;
  await t.frame(() => moveTo(t.surface(ref), x, y));
  const c = cells(t, ref);
  expect(c.hover).toEqual(c.under);
  expect(c.hover.length).toBe(2);
  return { x, y, el: t.surface(ref), at: c.hover };
}

for (const ui of ['react', 'dom'] as const) {
  describe(`the hover follows a resting pointer (${ui} UI)`, () => {
    it('a vertical wheel: the row under the pointer at the next frame, and the tooltip names it', async () => {
      const t = await open(ui);
      const { x, y, el, at } = await rest(t, 'v.lead');
      expect(tip(t)).toContain(`ROW${at[0]}`);
      let handled = false;
      await t.frame(() => {
        handled = wheel(el, x, y, 0, 3 * t.ctx.layout.props('v.lead').rowHeight);
      });
      expect(handled).toBe(true);
      const c = cells(t, 'v.lead');
      expect(c.hover).toEqual(c.under);
      expect(c.hover[0]).toBe(at[0] + 3);
      expect(c.hover[1]).toBe(at[1]);
      expect(tip(t)).toContain(`ROW${at[0] + 3}residue`);
    });

    it('a horizontal wheel, a Ctrl-wheel zoom, a programmatic zoom and scroll: the cell under the pointer each time', async () => {
      const t = await open(ui);
      const vp = t.ctx.viewports.get('v.lead');
      await t.frame(() => vp?.zoomTo(0, 10, 40));
      const { x, y, el, at } = await rest(t, 'v.lead');
      await t.frame(() => wheel(el, x, y, 80, 0));
      let c = cells(t, 'v.lead');
      expect(c.hover).toEqual(c.under);
      expect(c.hover[1]).toBeGreaterThan(at[1]);
      expect(tip(t)).toContain(`column ${(c.hover[1] as number) + 1}`);
      // Anchored at the pointer: the cell stays what is under it.
      await t.frame(() => wheel(el, x, y, 0, -120, true));
      c = cells(t, 'v.lead');
      expect(c.hover).toEqual(c.under);
      const before = c.hover;
      // Out, anchored at the data area's left: the pointer's column moves.
      await t.frame(() => vp?.zoomAt(0.6, 0, 'api'));
      c = cells(t, 'v.lead');
      expect(c.hover).toEqual(c.under);
      expect(c.hover[1]).not.toBe(before[1]);
      const s = vp?.live.current;
      const dy = 5 * t.ctx.layout.props('v.lead').rowHeight;
      await t.frame(() => vp?.scrollTo(s?.scrollX ?? 0, (s?.scrollY ?? 0) + dy, 'api'));
      c = cells(t, 'v.lead');
      expect(c.hover).toEqual(c.under);
      expect(c.hover[0]).toBe((before[0] as number) + 5);
      expect(tip(t)).toContain(`ROW${c.hover[0]}`);
      expect(tip(t)).toContain(`column ${(c.hover[1] as number) + 1}`);
    });

    it('a pointer that left gets no hover: a move after the leave hits nothing', async () => {
      const t = await open(ui);
      const { el } = await rest(t, 'v.lead');
      await t.frame(() => el.dispatchEvent(new PointerEvent('pointerleave', { isPrimary: true, pointerId: 1 })));
      expect(t.ctx.hover.live.current.kind).toBe(HoverKind.none);
      const v = t.ctx.hover.live.version.value;
      await t.frame(() => t.ctx.viewports.get('v.lead')?.panBy(0, 40, 'api'));
      expect([t.ctx.hover.live.current.kind, t.ctx.hover.live.version.value, tip(t)]).toEqual([HoverKind.none, v, '']);
    });

    it('a `follow`er moved by its leader: the hover in the follower, under the pointer, follows', async () => {
      const t = await open(ui, true);
      await t.frame(() => t.ctx.viewports.get('v.lead')?.zoomTo(0, 10, 40));
      const { at } = await rest(t, 'v.follow');
      await t.frame(() => t.ctx.viewports.get('v.lead')?.zoomTo(0, 50, 80));
      const c = cells(t, 'v.follow');
      expect(c.hover).toEqual(c.under);
      expect(c.hover[1]).not.toBe(at[1]);
      expect(tip(t)).toContain(`column ${(c.hover[1] as number) + 1}`);
    });
  });
}
