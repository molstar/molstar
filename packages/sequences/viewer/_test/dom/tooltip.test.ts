/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The tooltip's placement (UX2; seqspec:interaction#accessibility): from offsets up to its offset parent, right of and
// below the anchor, flipped left and kept as low as the viewer holds it where it would leave the viewer, read again at
// every call so it follows a surface that moved, and writing only a change; the pointer's anchor is the pointer, the
// keyboard's the cursor's cell.

import { afterEach, describe, expect, it } from 'vitest';
import * as D from '../../src/dom/index.js';
import '../../src/skin.css';
import { placeTooltip, tipSpot } from '../../src/dom/shared.js';
import type { SeqViewerContext } from '../../src/index.js';
import { cursorRect } from '../../src/keyboard.js';
import { at, flush, keys, openedContext, pointer, stubSurface, WIDTH, xOf } from './fixtures.js';

const cleanups: (() => void)[] = [];
afterEach(() => {
  for (const c of cleanups.splice(0)) c();
});

/** Where the tooltip was put, as the `transform` the skin makes of its two custom properties. */
const where = (el: HTMLElement) =>
  `translate(${el.style.getPropertyValue('--seqv-x')}, ${el.style.getPropertyValue('--seqv-y')})`;

/** A viewer's parent (relative, 400 × 300), a surface in it (left 20, top 30) and a tooltip (100 × 40). */
function stage() {
  const parent = document.createElement('div');
  parent.style.cssText = 'position:relative;width:400px;height:260px';
  const surface = document.createElement('div');
  surface.style.cssText = 'position:absolute;left:20px;top:30px;width:300px;height:200px';
  const tip = document.createElement('div');
  tip.style.cssText = 'position:absolute;left:0;top:0;width:100px;height:40px';
  parent.append(surface, tip);
  document.body.append(parent);
  cleanups.push(() => parent.remove());
  const state = { owner: 'pointer', p: { inside: true, section: 0, x: 0, y: 0 } };
  const ctx = {
    hover: { live: { current: state } },
    pointer: { current: state.p },
    layout: { sectionRef: () => 's' },
    ui: { keyCursor: { value: null } },
  } as unknown as SeqViewerContext;
  const spot = tipSpot();
  const place = (x: number, y: number) => {
    Object.assign(state.p, { x, y });
    placeTooltip(ctx, new Map([['s', surface]]), tip, spot);
    return where(tip);
  };
  return { parent, surface, tip, state, place, spot };
}

describe('placeTooltip from offsets', () => {
  it("goes right of and below the pointer, offset by the surface's place in the viewer", () => {
    const t = stage();
    expect(t.place(10, 10)).toBe('translate(44px, 58px)');
    expect(t.tip.style.visibility).toBe('');
  });

  it("flips left at the viewer's right edge, and at its bottom edge stays beside the pointer as low as it fits", () => {
    const t = stage();
    // 20 + 290 + 14 + 100 > 400: left of the pointer; 30 + 10 + 18 + 40 fits.
    expect(t.place(290, 10)).toBe('translate(198px, 58px)');
    // 30 + 190 + 18 + 40 > 260: beside the pointer, its bottom at the viewer's, never up over a ruler above the tracks.
    expect(t.place(10, 190)).toBe('translate(44px, 220px)');
    expect(t.place(290, 190)).toBe('translate(198px, 220px)');
  });

  it('goes above the pointer at the bottom of a viewer too narrow to hold it beside the pointer', () => {
    const t = stage();
    t.parent.style.width = '120px';
    // Neither right (20 + 10 + 14 + 100 > 120) nor left of the pointer fits: kept in at 0, over the pointer across, so
    // it goes above it (30 + 190 - 12 - 40), not down over what it labels.
    expect(t.place(10, 190)).toBe('translate(0px, 168px)');
  });

  it('stays inside the viewer when it is larger than it', () => {
    const t = stage();
    t.parent.style.width = '60px';
    expect(t.place(10, 10)).toBe('translate(0px, 58px)');
  });

  it('writes only when the place changes: a move that rounds to the same place and a repeated call write nothing', async () => {
    const t = stage();
    t.place(10, 10);
    const writes: number[] = [];
    const mo = new MutationObserver((r) => writes.push(r.length));
    mo.observe(t.tip, { attributes: true });
    t.place(10, 10);
    t.place(10.2, 10.1);
    await Promise.resolve();
    expect(mo.takeRecords().length).toBe(0);
    t.place(11, 10);
    expect(mo.takeRecords().length).toBeGreaterThan(0);
    mo.disconnect();
  });

  it('is hidden when no pointer is over a surface, and shown again, writing the visibility only on the change', async () => {
    const t = stage();
    t.place(10, 10);
    const mo = new MutationObserver(() => {});
    mo.observe(t.tip, { attributes: true });
    t.state.p.inside = false;
    placeTooltip(
      {
        hover: { live: { current: t.state } },
        pointer: { current: t.state.p },
        layout: { sectionRef: () => 's' },
        ui: { keyCursor: { value: null } },
      } as unknown as SeqViewerContext,
      new Map([['s', t.surface]]),
      t.tip,
      t.spot,
    );
    expect(t.tip.style.visibility).toBe('hidden');
    expect(mo.takeRecords().length).toBe(1);
    t.state.p.inside = true;
    expect(t.place(10, 10)).toBe('translate(44px, 58px)');
    expect(t.tip.style.visibility).toBe('');
    mo.disconnect();
  });

  it('follows the surface when a part above it grows, with no cached geometry', () => {
    const t = stage();
    expect(t.place(10, 10)).toBe('translate(44px, 58px)');
    t.surface.style.top = '80px';
    expect(t.place(10, 10)).toBe('translate(44px, 108px)');
  });
});

async function mounted() {
  const t = await openedContext();
  cleanups.push(() => t.ctx.dispose());
  const host = document.createElement('div');
  host.style.width = `${WIDTH}px`;
  document.body.appendChild(host);
  cleanups.push(() => host.remove());
  const m = D.mountSeqViewer(host, t.ctx, { surface: stubSurface });
  cleanups.unshift(() => m.dispose());
  await flush();
  return { ...t, host };
}
const move = (ctx: SeqViewerContext, x: number, y: number) => pointer(ctx)('move', 'chains.main', x, y);
const translate = (el: HTMLElement) => (where(el).match(/-?[\d.]+/g) ?? []).map(Number);

describe('the tooltip in a viewer', () => {
  it('follows the surface when a part above the sections grows between two hovers', async () => {
    const { host, ctx } = await mounted();
    const res = at(ctx, 'chains.main', 't.res');
    move(ctx, xOf(ctx, 'chains.main', 12), res.y);
    const tip = host.querySelector('.seqv-tooltip') as HTMLElement;
    const [x, y] = translate(tip);
    const surface = host.querySelector('.seqv-surface') as HTMLElement;
    const top = surface.getBoundingClientRect().top;
    // A part above the sections grows: a banner before the view (and the viewer's gap after it).
    const banner = document.createElement('div');
    banner.style.height = '50px';
    (host.querySelector('.seqv-root') as HTMLElement).prepend(banner);
    const grown = surface.getBoundingClientRect().top - top;
    expect(grown).toBeGreaterThanOrEqual(50);
    move(ctx, xOf(ctx, 'chains.main', 12) + 2, res.y);
    expect(translate(tip)).toEqual([x + 2, y + grown]);
  });

  it("is anchored at the keyboard cursor's cell, and follows the cursor along the track", async () => {
    const { host, ctx } = await mounted();
    const key = keys(ctx);
    ctx.input.focus(ctx.layout.sectionIndex('chains.main'), true);
    expect(ctx.hover.live.current.owner).toBe('keyboard');
    const tip = host.querySelector('.seqv-tooltip') as HTMLElement;
    expect(tip).not.toBeNull();
    const rect = new Float64Array(6);
    const root = host.querySelector('.seqv-root') as HTMLElement;
    const surface = host.querySelector('.seqv-surface [tabindex]') as HTMLElement;
    const expected = () => {
      expect(cursorRect(ctx, ctx.ui.keyCursor.value as never, rect)).toBe(true);
      const [sx, sy] = [
        surface.getBoundingClientRect().left - root.getBoundingClientRect().left,
        surface.getBoundingClientRect().top - root.getBoundingClientRect().top,
      ];
      return [Math.round(sx + rect[0] + rect[2] + 6), Math.round(sy + rect[1] + rect[3] + 4)];
    };
    expect(tip.style.visibility).toBe('');
    expect(translate(tip)).toEqual(expected());
    const before = translate(tip);
    key('chains.main', 'ArrowRight');
    key('chains.main', 'ArrowRight');
    expect(translate(tip)).toEqual(expected());
    expect(translate(tip)[0]).toBeGreaterThan(before[0]);
    // The pointer's tooltip is not drawn at the pointer while the keyboard owns the hover.
    ctx.input.focus(ctx.layout.sectionIndex('chains.main'), false);
    await flush();
    expect(host.querySelector('.seqv-tooltip')).toBeNull();
  });

  it("never blocks the canvas's own hover: it takes no pointer events, the surface under it gets them", async () => {
    const { host, ctx } = await mounted();
    const res = at(ctx, 'chains.main', 't.res');
    move(ctx, xOf(ctx, 'chains.main', 12), res.y);
    const tip = host.querySelector('.seqv-tooltip') as HTMLElement;
    expect(getComputedStyle(tip).pointerEvents).toBe('none');
    const r = tip.getBoundingClientRect();
    const under = document.elementFromPoint(r.left + r.width / 2, r.top + r.height / 2);
    expect(under && !tip.contains(under) && host.contains(under)).toBe(true);
  });

  it('is hidden over a pending region or off the surfaces: nothing is drawn without a label', async () => {
    const { host, ctx } = await mounted();
    move(ctx, 3, -5);
    expect(host.querySelector('.seqv-tooltip')).toBeNull();
    ctx.input.leave(ctx.layout.sectionIndex('chains.main'));
    expect(host.querySelector('.seqv-tooltip')).toBeNull();
  });
});
