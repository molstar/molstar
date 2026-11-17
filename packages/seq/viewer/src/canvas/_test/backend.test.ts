/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Anchor, DrawLists, lineTo, moveTo, Round } from '@molstar/seqx-core/render';
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
  const images = new FakeCanvas();
  const backend = new CanvasBackend(base, overlay, { createCanvas: () => images });
  backend.resize(200, 100, 2);
  base.context.take();
  overlay.context.take();
  return { base, overlay, images, backend };
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
    const b = new CanvasBackend(base, overlay);
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
    const { base, backend } = setup();
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

  it('draws the base passes on the base canvas and the emphasis passes on the overlay', () => {
    const { base, overlay, backend } = setup();
    const l = lists();
    l.writer('hover').rect(1, 1, 5, 5, 0, BLUE, 1);
    l.writer('labels').text('x', 0, 0, Anchor.start, '12px sans', RED, 5);
    backend.draw(l, true);
    expect(base.context.calls('fillText').length).toBe(1);
    expect(base.context.calls('strokeRect').length).toBe(0);
    expect(overlay.context.calls('strokeRect')).toEqual(['strokeRect(1, 1, 5, 5)']);
  });

  it('a hover redraws the overlay only: the base canvas is not touched', () => {
    const { base, overlay, backend } = setup();
    const l = lists();
    backend.draw(l, true);
    expect(backend.counts).toEqual({ base: 1, overlay: 0 });
    base.context.take();
    overlay.context.take();
    // the hover appears: only the emphasis passes are rebuilt (the renderer resets them alone)
    l.reset(['highlight', 'selection', 'hover', 'cursor']);
    l.writer('hover').rect(1, 1, 5, 5, 0, BLUE, 1);
    backend.draw(l, false);
    expect(base.context.take()).toEqual([]);
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
    expect(backend.counts).toEqual({ base: 1, overlay: 3 });
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
});
