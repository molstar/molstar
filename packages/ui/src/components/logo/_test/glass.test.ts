/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { beforeEach, describe, expect, it } from 'vitest';
import { Glass, type GlassColors } from '../glass.js';
import { fakeContext, stubPath2D } from './fake-canvas.js';

const colors: GlassColors = { logo: ['#111', '#222', '#333', '#444'], edge: '#999' };
const Quick = { duration: 0.9, stagger: 0.3 };
const run = (g: Glass, seconds: number) => {
  for (let t = 0; t < seconds; t += 1 / 60) g.step(1 / 60);
};

describe('the glass wordmark', () => {
  beforeEach(() => stubPath2D());

  it('starts whole and resting; scattered, it assembles and rests again (no frames needed)', () => {
    const g = new Glass(23, 0.35);
    expect(g.resting).toBe(true);
    g.scatter(Quick);
    expect(g.phase).toBe('assembling');
    expect(g.resting).toBe(false);
    run(g, 1.3);
    expect(g.phase).toBe('whole');
    expect(g.resting).toBe(true);
  });

  it('a press shatters it, it holds apart, then snaps back whole', () => {
    const g = new Glass(23, 0.35);
    g.shatter(20, 10);
    expect(g.phase).toBe('broken');
    run(g, 0.5);
    expect(g.phase).toBe('broken');
    run(g, 1);
    expect(g.phase).toBe('whole');
  });

  it('a pointer resting over it: the pieces lean and settle, then it rests (no frames) while drawn leaned', () => {
    const g = new Glass(23, 0.35);
    const { ctx, counts } = fakeContext();
    g.setPointer(20, 10);
    g.step(1 / 60);
    expect(g.resting).toBe(false);
    run(g, 2);
    // Settled under the pointer: nothing moves, nothing asks for a frame...
    expect(g.resting).toBe(true);
    // ...yet it is drawn leaned, in pieces (not the word whole).
    counts.stroke = 0;
    g.draw(ctx, 0, 0, 1, colors);
    expect(counts.drawImage + counts.stroke).toBeGreaterThan(4);
    // The pointer leaves: it springs back and is drawn whole again.
    g.clearPointer();
    run(g, 2);
    counts.stroke = 0;
    counts.drawImage = 0;
    g.draw(ctx, 0, 0, 1, colors);
    expect(counts.drawImage).toBe(0);
    expect(counts.stroke).toBe(4);
  });

  it('makes the logo gradient once for the same colors, again when they change', () => {
    const g = new Glass(23, 0.35);
    const { ctx, counts } = fakeContext();
    g.draw(ctx, 0, 0, 1, colors);
    g.draw(ctx, 0, 0, 1, colors);
    g.draw(ctx, 0, 0, 1, { ...colors });
    expect(counts.gradients).toBe(1);
    g.draw(ctx, 0, 0, 1, { ...colors, logo: ['#000', '#222', '#333', '#444'] });
    expect(counts.gradients).toBe(2);
  });
});
