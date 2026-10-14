/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// CG3 in Chromium (seqspec:interaction#linked-sections): a ring's surface follows the linear section its window shows:
// a move of that section (a pan, a zoom) asks the ring for a frame and redraws its overlay, never its tracks; the
// window's sector is on the overlay canvas where the linear view is, and leaves it when the view moves away.

import { afterEach, describe, expect, it } from 'vitest';
import { SectionSurface } from '../../src/canvas/index.js';
import { Pass } from '../../src/draw/backend.js';
import { createSeqViewer, DefaultSeqViewerSpec, manualScheduler, type SeqViewerContext } from '../../src/index.js';
import { Ring } from '../../src/ring.js';
import { fixtureMeasurer } from '../../src/testing/index.js';
import type { Viewport } from '../../src/viewport.js';
import { ringWindowDoc } from '../ring-window-fixtures.js';

const open: { ctx: SeqViewerContext; host: HTMLElement; surface: SectionSurface }[] = [];
afterEach(() => {
  for (const t of open.splice(0)) {
    t.surface.dispose();
    t.ctx.dispose();
    t.host.remove();
  }
});

describe("a ring's surface and its window", () => {
  it('redraws its overlay when the linked linear section moves, the sector where its view is', async () => {
    const scheduler = manualScheduler();
    const ctx = createSeqViewer(DefaultSeqViewerSpec(), { scheduler, text: fixtureMeasurer() });
    const r = await ctx.setDocument(ringWindowDoc());
    if (!r.ok) throw new Error(JSON.stringify(r.issues));
    ctx.viewports.get('lin')?.setSize(1000, 200, 1);
    await ctx.graph.update();
    const host = document.createElement('div');
    host.style.cssText = 'position:absolute;left:0;top:0;width:600px;height:520px';
    document.body.appendChild(host);
    const surface = new SectionSurface(ctx, 'ring', host, { dpr: 1, readback: true });
    open.push({ ctx, host, surface });
    scheduler.advance(0);
    const lin = ctx.viewports.get('lin') as Viewport | null;
    const l = ctx.layout.get('ring')?.current;
    const av = ctx.layout.axis('ring');
    if (!lin || !l || !av) throw new Error('no sections');
    // The overlay's alpha at ring position `pos`, halfway down the features' ring.
    const overlay = surface.element.querySelectorAll('canvas')[Pass.overlay];
    const g = overlay.getContext('2d') as CanvasRenderingContext2D;
    const ring = new Ring().from(l, av, 1);
    const out = new Float64Array(2);
    const alpha = (pos: number) => {
      ring.toSurface(av.blockX[0] + pos * av.span[0] * av.cell, l.scrolledTop + l.entryTop[1] + 4, out);
      return g.getImageData(Math.round(out[0]), Math.round(out[1]), 1, 1).data[3];
    };
    const frames = surface.stats.frames;
    const tracks = surface.composer.stats.fullDraws;
    // The linear view on [100, 160): asked once, drawn as an overlay frame only.
    const area = lin.live.current.width - (lin.axis?.insetX ?? 0);
    lin.setView(area / 60, 100 * (area / 60), 'api');
    expect(scheduler.pendingFrames).toBe(1);
    scheduler.advance(0);
    expect(surface.stats.frames).toBe(frames + 1);
    expect(surface.composer.stats.fullDraws).toBe(tracks);
    expect([alpha(130) > 0, alpha(400)]).toEqual([true, 0]);
    // Panned to [380, 440): the sector follows.
    lin.setView(area / 60, 380 * (area / 60), 'api');
    scheduler.advance(0);
    expect([alpha(130), alpha(410) > 0]).toEqual([0, true]);
  });
});
