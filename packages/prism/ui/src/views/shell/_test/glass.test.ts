/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { Glass } from '../glass.js';

// The splash's wordmark as glass: it assembles, shatters when pressed, and comes back together.

/** Node has no Path2D: a shape that records nothing will do for the states. */
class NoPath {
  moveTo() {}
  lineTo() {}
  closePath() {}
}

beforeEach(() => vi.stubGlobal('Path2D', NoPath));
afterEach(() => vi.unstubAllGlobals());

function run(glass: Glass, seconds: number) {
  for (let t = 0; t < seconds; t += 1 / 60) glass.step(1 / 60);
}

describe('the wordmark as glass', () => {
  it('assembles from pieces, then is whole and says so', () => {
    const glass = new Glass();
    const whole: boolean[] = [];
    let snaps = 0;
    glass.onWhole = (w) => whole.push(w);
    glass.onSnap = () => snaps++;
    glass.scatter();
    expect(glass.phase).toBe('assembling');
    expect(glass.whole).toBe(false);
    // Smooth: a while (its pieces' starts spread).
    run(glass, 1.5);
    expect(glass.whole).toBe(false);
    run(glass, 2);
    expect(glass.whole).toBe(true);
    // In pieces first (the page's wordmark hides), then whole.
    expect(whole).toEqual([false, true]);
    expect(snaps).toBe(1);
  });

  it('shatters when pressed on the word, comes back together after a moment', () => {
    const glass = new Glass();
    const whole: boolean[] = [];
    glass.onWhole = (w) => whole.push(w);
    glass.scatter();
    run(glass, 3.5);
    expect(glass.hits(20, 10)).toBe(true);
    expect(glass.hits(200, 10)).toBe(false);
    glass.shatter(20, 10);
    expect(glass.phase).toBe('broken');
    expect(whole).toEqual([false, true, false]);
    run(glass, 0.5);
    expect(glass.phase).toBe('broken');
    run(glass, 0.5);
    expect(glass.phase).toBe('assembling');
    // A snap back: quick.
    run(glass, 0.6);
    expect(glass.whole).toBe(true);
    expect(whole).toEqual([false, true, false, true]);
  });
});
