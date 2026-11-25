/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Anchor, DrawLists, drawListsToSvg, lineTo, moveTo, Passes, Round } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import { CanvasBackend } from '../backend.js';
import { canvasTextMetrics } from '../metrics.js';
import { cssFont } from '../font.js';
import { FakeCanvas } from './recording.js';

const RED = 0xff0000ff;
const BLUE = 0x0000ff80;

function setup() {
  const base = new FakeCanvas();
  const overlay = new FakeCanvas();
  const labels = new FakeCanvas();
  const images = new FakeCanvas();
  const backend = new CanvasBackend(base, overlay, labels, { createCanvas: () => images });
  backend.resize(200, 100, 2);
  base.context.take();
  overlay.context.take();
  labels.context.take();
  return { base, overlay, labels, images, backend };
}

function lists(): DrawLists {
  const l = new DrawLists();
  l.width = 200;
  l.height = 100;
  l.writer('background').rect(0, 0, 200, 100, RED);
  return l;
}

describe('CanvasBackend', () => {
  it('sizes both canvases for the device pixel ratio and sets the CSS size', () => {
    const base = new FakeCanvas();
    const overlay = new FakeCanvas();
    const b = new CanvasBackend(base, overlay, new FakeCanvas());
    expect(b.resize(200, 100, 1.5)).toBe(true);
    expect([base.width, base.height, overlay.width, overlay.height]).toEqual([300, 150, 300, 150]);
    expect(base.style).toEqual({ width: '200px', height: '100px' });
    expect(b.resize(200, 100, 1.5)).toBe(false);
  });

  it('draws rects with their fill and stroke, in CSS px under the device transform', () => {
    const { base, backend } = setup();
    const l = lists();
    l.writer('data').rect(10, 20, 30, 40, RED, BLUE, 2);
    backend.draw(l, true);
    expect(base.context.take()).toEqual([
      'setTransform(2, 0, 0, 2, 0, 0)',
      'textBaseline=alphabetic',
      'fillStyle=rgba(255, 0, 0, 1)',
      'fillRect(0, 0, 200, 100)',
      'fillStyle=rgba(255, 0, 0, 1)',
      'fillRect(10, 20, 30, 40)',
      'strokeStyle=rgba(0, 0, 255, 0.502)',
      'lineWidth=2',
      'strokeRect(10, 20, 30, 40)',
    ]);
  });

  it('draws a squircle through the core’s curve: a cubic a rounded corner', () => {
    const { base, backend } = setup();
    const l = new DrawLists();
    l.writer('data').rect(0, 0, 40, 20, RED, 0, 0, 6, Round.both);
    backend.draw(l, true);
    expect(base.context.calls('bezierCurveTo').length).toBe(4);
    const l2 = new DrawLists();
    l2.writer('data').rect(0, 0, 40, 20, RED, 0, 0, 6, Round.left);
    base.context.take();
    backend.draw(l2, true);
    expect(base.context.calls('bezierCurveTo').length).toBe(2);
  });

  it('draws clips as a saved, clipped state restored at the unclip, and a pass left open is closed', () => {
    const { base, backend } = setup();
    const l = new DrawLists();
    const out = l.writer('data');
    out.clip(5, 5, 50, 20);
    out.rect(0, 0, 100, 100, RED);
    out.unclip();
    out.clip(1, 1, 2, 2);
    backend.draw(l, true);
    const log = base.context.take();
    expect(log.filter((c) => c === 'save()').length).toBe(2);
    expect(log.filter((c) => c === 'restore()').length).toBe(2);
    expect(log).toContain('rect(5, 5, 50, 20)');
  });

  it('a clip inside a clip intersects it: each is its own saved state over the one around it, restored innermost first', () => {
    const { base, backend } = setup();
    const l = new DrawLists();
    const out = l.writer('data');
    out.clip(0, 0, 100, 50);
    out.clip(20, 30, 100, 50);
    out.rect(0, 0, 200, 200, RED);
    out.unclip();
    out.unclip();
    backend.draw(l, true);
    const log = base.context.take();
    // the second clip is made with the first still in force (the canvas's clip() intersects with the current one), and
    // the first's state is restored only after the second's
    const first = log.indexOf('rect(0, 0, 100, 50)');
    const second = log.indexOf('rect(20, 30, 100, 50)');
    expect(first).toBeGreaterThan(-1);
    expect(second).toBeGreaterThan(first);
    expect(log.slice(first, second).includes('restore()')).toBe(false);
    expect(log.indexOf('restore()')).toBeGreaterThan(log.indexOf('fillRect(0, 0, 200, 200)'));
    expect(log.filter((c) => c === 'restore()').length).toBe(2);
  });

  it('draws a squircle clip (a track’s well) through the core’s curve, as the SVG’s clip path; a box clip a rect', () => {
    const { base, backend } = setup();
    const l = new DrawLists();
    const out = l.writer('data');
    out.clip(5, 5, 60, 20, 6, Round.both);
    out.rect(0, 0, 100, 100, RED);
    out.unclip();
    backend.draw(l, true);
    const log = base.context.take();
    // the clip's path: four corners, then the clip; no box
    const at = log.indexOf('clip()');
    expect(at).toBeGreaterThan(0);
    expect(log.slice(0, at).filter((c) => c.startsWith('bezierCurveTo')).length).toBe(4);
    expect(log).not.toContain('rect(5, 5, 60, 20)');
    const svg = drawListsToSvg(l);
    expect(svg).toMatch(/<clipPath id="[^"]*"><path d="M 11 5 L [^"]*C[^"]*"\/><\/clipPath>/);
    // a clip with no radius stays a box in both
    const l2 = new DrawLists();
    l2.writer('data').clip(5, 5, 60, 20);
    backend.draw(l2, true);
    expect(base.context.take()).toContain('rect(5, 5, 60, 20)');
    expect(drawListsToSvg(l2)).toContain('<clipPath id="c0"><rect x="5" y="5" width="60" height="20"/></clipPath>');
  });

  it('draws an emphasis squircle as the SVG does: an outline concentric, square where cut', () => {
    const { overlay, backend } = setup();
    const l = new DrawLists();
    // a hover's outline round at its left end, cut at its right (the view's edge)
    l.writer('hover').rect(10, 10, 30, 24, 0, RED, 1.5, 8.42, Round.left);
    backend.draw(l, false);
    expect(overlay.context.calls('bezierCurveTo').length).toBe(2);
    expect(overlay.context.calls('stroke').length).toBe(1);
    expect(
      (
        drawListsToSvg(l)
          .match(/<path d="[^"]*"/)?.[0]
          .match(/C/g) ?? []
      ).length,
    ).toBe(2);
  });

  it('draws glyph runs a letter at a time, each in its own color or one for the run', () => {
    const { base, backend } = setup();
    const l = new DrawLists();
    const out = l.writer('data');
    out.beginGlyphs('14px mono', 30, 0);
    out.letter(77, 10, RED);
    out.letter(69, 28, BLUE);
    out.endGlyphs();
    out.beginGlyphs('bold 10px sans', 50, RED);
    out.letter(65, 10);
    out.endGlyphs();
    backend.draw(l, true);
    const log = base.context.take();
    expect(log).toContain(`font=${cssFont('14px mono')}`);
    expect(log).toContain('fillText("M", 10, 30)');
    expect(log).toContain('fillText("E", 28, 30)');
    expect(log.indexOf('fillStyle=rgba(255, 0, 0, 1)')).toBeLessThan(log.indexOf('fillText("M", 10, 30)'));
    expect(log).toContain(`font=${cssFont('bold 10px sans')}`);
    expect(log.filter((c) => c === 'textAlign=center').length).toBe(2);
  });

  it('fades a text past its max width through a gradient clip, and draws one that fits plainly', () => {
    const { labels: base, backend } = setup();
    const l = new DrawLists();
    const out = l.writer('labels');
    out.text('short', 8, 20, Anchor.start, '12px sans', RED, 30, 100);
    out.text('a very long label', 8, 40, Anchor.start, '12px sans', RED, 150, 100);
    backend.draw(l, true);
    const log = base.context.take();
    expect(log.filter((c) => c.startsWith('createLinearGradient'))).toEqual(['createLinearGradient(96, 0, 108, 0)']);
    expect(log.filter((c) => c.startsWith('fillText'))).toEqual([
      'fillText("short", 8, 20)',
      'fillText("a very long label", 8, 40)',
    ]);
  });

  it('draws a path from its commands, dashed when asked and the dash undone after', () => {
    const { base, backend } = setup();
    const l = new DrawLists();
    l.writer('data').path(
      (p) => {
        moveTo(p, 1, 2);
        lineTo(p, 3, 4);
      },
      0,
      RED,
      1,
      4,
    );
    backend.draw(l, true);
    const log = base.context.take();
    expect(log).toContain('moveTo(1, 2)');
    expect(log).toContain('lineTo(3, 4)');
    expect(log).toContain('setLineDash(4/4)');
    expect(log.at(-1)).toBe('setLineDash()');
  });

  it('draws an image through one reusable pixel buffer, copied at each draw', () => {
    const { base, images, backend } = setup();
    const bytes = new Uint8Array(2 * 2 * 4).fill(9);
    const l = new DrawLists();
    l.writer('data').image(10, 10, 20, 20, { width: 2, height: 2, data: bytes });
    backend.draw(l, true);
    backend.draw(l, true);
    expect(images.context.calls('createImageData').length).toBe(1);
    expect(images.context.calls('putImageData').length).toBe(2);
    expect(base.context.calls('drawImage')).toEqual(['drawImage(10, 10, 20, 20)', 'drawImage(10, 10, 20, 20)']);
  });

  it('draws the background and data on the base canvas, the emphasis on the overlay, the labels above both', () => {
    const { base, overlay, labels, backend } = setup();
    const l = lists();
    l.writer('hover').rect(1, 1, 5, 5, 0, BLUE, 1);
    l.writer('labels').text('x', 0, 0, Anchor.start, '12px sans', RED, 5);
    backend.draw(l, true);
    expect(labels.context.calls('fillText').length).toBe(1);
    expect(labels.context.calls('clearRect').length).toBe(1);
    expect(base.context.calls('fillText').length).toBe(0);
    expect(base.context.calls('strokeRect').length).toBe(0);
    expect(overlay.context.calls('strokeRect')).toEqual(['strokeRect(1, 1, 5, 5)']);
  });

  it('a hover redraws the overlay only: the base canvas is not touched', () => {
    const { base, overlay, labels, backend } = setup();
    const l = lists();
    backend.draw(l, true);
    expect(backend.counts).toEqual({ base: 1, overlay: 0, labels: 1, tip: 0 });
    base.context.take();
    labels.context.take();
    overlay.context.take();
    // the hover appears: only the emphasis passes are rebuilt (the renderer resets them alone)
    l.reset(['highlight', 'selection', 'hover', 'cursor']);
    l.writer('hover').rect(1, 1, 5, 5, 0, BLUE, 1);
    backend.draw(l, false);
    expect(base.context.take()).toEqual([]);
    expect(labels.context.take()).toEqual([]);
    expect(overlay.context.calls('clearRect')).toEqual(['clearRect(0, 0, 200, 100)']);
    expect(overlay.context.calls('strokeRect').length).toBe(1);
    // it moves: the overlay is cleared and drawn again
    overlay.context.take();
    l.reset(['highlight', 'selection', 'hover', 'cursor']);
    l.writer('hover').rect(9, 9, 5, 5, 0, BLUE, 1);
    backend.draw(l, false);
    expect(base.context.take()).toEqual([]);
    expect(overlay.context.calls('strokeRect')).toEqual(['strokeRect(9, 9, 5, 5)']);
    // it goes: cleared once, then nothing at all
    overlay.context.take();
    l.reset(['highlight', 'selection', 'hover', 'cursor']);
    backend.draw(l, false);
    expect(overlay.context.calls('clearRect').length).toBe(1);
    overlay.context.take();
    backend.draw(l, false);
    expect(overlay.context.take()).toEqual([]);
    expect(base.context.take()).toEqual([]);
    expect(backend.counts).toEqual({ base: 1, overlay: 3, labels: 1, tip: 0 });
  });

  it('draws with the backing size over the CSS size: a rounded canvas is not drawn at the ratio asked for', () => {
    const base = new FakeCanvas();
    const b = new CanvasBackend(base, new FakeCanvas(), new FakeCanvas());
    b.resize(101, 50, 1.5);
    expect([base.width, base.height]).toEqual([152, 75]);
    base.context.take();
    const l = new DrawLists();
    l.writer('background').rect(0, 0, 101, 50, RED);
    b.draw(l, true);
    expect(base.context.take()[0]).toBe(`setTransform(${Math.round((152 / 101) * 100) / 100}, 0, 0, 1.5, 0, 0)`);
  });

  it('keeps a slot per image strip: a strip’s buffer is kept while its size stays, made again in place when it changes', () => {
    const made: FakeCanvas[] = [];
    const base = new FakeCanvas();
    const b = new CanvasBackend(base, new FakeCanvas(), new FakeCanvas(), {
      createCanvas: () => {
        const c = new FakeCanvas();
        made.push(c);
        return c;
      },
    });
    b.resize(100, 100, 1);
    const strip = (w: number, h: number) => ({ width: w, height: h, data: new Uint8Array(w * h * 4) });
    const draw = (sizes: number[]) => {
      const l = new DrawLists();
      sizes.forEach((n, i) => l.writer('data').image(0, i * 10, 50, 8, strip(n, 1)));
      b.draw(l, true);
    };
    draw([4, 8]);
    draw([4, 8]);
    expect(made.length).toBe(2);
    expect(made.map((c) => c.context.calls('createImageData').length)).toEqual([1, 1]);
    // the second strip changes its size: its canvas is resized, no new canvas
    draw([4, 16]);
    expect(made.length).toBe(2);
    expect(made[1].width).toBe(16);
    expect(made[1].context.calls('createImageData').length).toBe(2);
    // more strips than slots: the extra share the last one, nothing grows without bound
    draw(Array.from({ length: 40 }, () => 4));
    expect(made.length).toBeLessThanOrEqual(16);
  });

  it('the canvases stack as the SVG’s passes do: the same order, pass by pass', () => {
    const base = new FakeCanvas();
    const overlay = new FakeCanvas();
    const labels = new FakeCanvas();
    const tip = new FakeCanvas();
    const backend = new CanvasBackend(base, overlay, labels, { tip });
    backend.resize(200, 100, 2);
    // one rect a pass, each a color of its own whose red channel is the pass's place
    const l = new DrawLists();
    Passes.forEach((name, i) => l.writer(name).rect(1, 1, 5, 5, ((10 + i) << 24) | 0xff));
    backend.draw(l, true);
    const red = (log: string[]) =>
      log
        .map((c) => /fillStyle=rgba\((\d+),/.exec(c)?.[1])
        .filter((x) => x !== undefined)
        .map(Number);
    // bottom to top: the base canvas, the overlay, the labels, the tip
    const canvasOrder = [
      ...red(base.context.log),
      ...red(overlay.context.log),
      ...red(labels.context.log),
      ...red(tip.context.log),
    ];
    const svg = drawListsToSvg(l);
    const svgOrder = [...svg.matchAll(/data-pass="(\w+)"/g)].map(
      (m) => Passes.indexOf(m[1] as (typeof Passes)[number]) + 10,
    );
    expect(svgOrder).toEqual(Passes.map((_, i) => i + 10));
    expect(canvasOrder).toEqual(svgOrder);
  });

  it('keeps colors and glyph strings: drawing the same lists twice makes the same calls', () => {
    const { base, backend } = setup();
    const l = lists();
    const out = l.writer('data');
    out.beginGlyphs('14px mono', 30, 0);
    for (let i = 0; i < 20; i++) out.letter(65 + i, 10 + i * 10, RED + i * 256);
    out.endGlyphs();
    backend.draw(l, true);
    const first = base.context.take();
    backend.draw(l, true);
    expect(base.context.take()).toEqual(first);
  });
});

describe('canvasTextMetrics', () => {
  it('measures in the CSS font the backend draws with, setting the font only when it changes', () => {
    const sets: string[] = [];
    const ctx = {
      set font(v: string) {
        sets.push(v);
      },
      measureText: (t: string) => ({ width: t.length * 6 }),
    };
    const port = canvasTextMetrics(ctx as never);
    expect(port.measure('12px sans', 'abcd')).toBe(24);
    expect(port.measure('12px sans', 'ab')).toBe(12);
    expect(port.measure('bold 10px mono', 'a')).toBe(6);
    expect(sets).toEqual([cssFont('12px sans'), cssFont('bold 10px mono')]);
  });

  it('reads a font\u2019s cap height from the ascent of H, in the font it is set in', () => {
    let font = '';
    const ctx = {
      set font(v: string) {
        font = v;
      },
      measureText: (t: string) => ({
        width: 0,
        actualBoundingBoxAscent: t === 'H' ? (font.includes('mono') ? 70 : 72) : 0,
      }),
    };
    const port = canvasTextMetrics(ctx as never, { sans: 'sans-serif', mono: 'monospace' } as never);
    expect(port.capHeight!('100px sans')).toBe(72);
    expect(port.capHeight!('100px mono')).toBe(70);
  });

  it('the tip pass is drawn above the labels: on a canvas of its own, redrawn alone with the hover; else on the overlay', () => {
    const base = new FakeCanvas();
    const overlay = new FakeCanvas();
    const labels = new FakeCanvas();
    const tip = new FakeCanvas();
    const backend = new CanvasBackend(base, overlay, labels, { tip });
    backend.resize(200, 100, 2);
    // no tip yet: no backing store of the view's size
    expect([tip.width, tip.height]).toEqual([1, 1]);
    const l = lists();
    l.writer('tip').rect(1, 1, 5, 5, 0, BLUE, 1);
    backend.draw(l, true);
    expect([tip.width, tip.height]).toEqual([400, 200]);
    expect(tip.context.calls('strokeRect')).toEqual(['strokeRect(1, 1, 5, 5)']);
    expect(overlay.context.calls('strokeRect')).toEqual([]);
    expect(labels.context.calls('strokeRect')).toEqual([]);
    // the hover goes: the tip canvas shrinks back (cleared), the labels are not drawn again
    const clears = labels.context.calls('clearRect').length;
    l.reset(['highlight', 'selection', 'hover', 'cursor', 'tip']);
    backend.draw(l, false);
    expect(labels.context.calls('clearRect').length).toBe(clears);
    expect([tip.width, tip.height]).toEqual([1, 1]);
    expect(backend.counts.tip).toBe(2);
    // and nothing more to do while there is no tip
    backend.draw(l, false);
    expect(backend.counts.tip).toBe(2);
    // without a tip canvas the pass is the overlay's last
    const o2 = new FakeCanvas();
    const b2 = new CanvasBackend(new FakeCanvas(), o2, new FakeCanvas());
    b2.resize(200, 100, 2);
    b2.draw(l, true);
    l.writer('tip').rect(1, 1, 5, 5, 0, BLUE, 1);
    b2.draw(l, false);
    expect(o2.context.calls('strokeRect')).toEqual(['strokeRect(1, 1, 5, 5)']);
  });
});
