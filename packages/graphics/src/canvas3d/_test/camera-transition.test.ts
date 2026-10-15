/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Camera } from '../camera.js';

// A camera transition runs on the time its ticks bring. The plugin's animation loop restarts its clock at 0 (after an
// image is drawn, for one), so a transition may see the time go back: it keeps going from where it was rather than
// extrapolating far outside its source and target.

function camera() {
  const c = new Camera({
    position: Vec3.create(0, 0, 100),
    target: Vec3.create(0, 0, 0),
    radius: 30,
    radiusMax: 50,
  });
  return c;
}

describe('camera transition', () => {
  it('interpolates between its source and target', () => {
    const c = camera();
    c.transition.tick(1000);
    c.setState({ radius: 10, position: Vec3.create(0, 0, 40) }, 200);
    c.transition.tick(1100);
    expect(c.state.radius).toBeGreaterThan(10);
    expect(c.state.radius).toBeLessThan(30);
    c.transition.tick(1200);
    expect(c.state.radius).toBe(10);
    expect(c.transition.inTransition).toBe(false);
  });

  it('keeps going from where it was when the clock restarts during it', () => {
    const c = camera();
    c.transition.tick(26000);
    c.setState({ radius: 10, position: Vec3.create(0, 0, 40) }, 250);
    c.transition.tick(26050);
    const before = c.state.radius;
    // The loop's clock starts again at 0.
    c.transition.tick(0);
    expect(c.state.radius).toBeGreaterThanOrEqual(10);
    expect(c.state.radius).toBeLessThanOrEqual(before);
    c.transition.tick(100);
    expect(c.state.radius).toBeGreaterThanOrEqual(10);
    expect(c.state.radius).toBeLessThan(before);
    // The rest of its 250 ms on the new clock: done.
    c.transition.tick(250);
    expect(c.state.radius).toBe(10);
    expect(c.transition.inTransition).toBe(false);
  });
});

describe('camera settled snapshot', () => {
  it('is the transition target while one runs, else the camera as it is', () => {
    const c = camera();
    expect(c.getSettledSnapshot().radius).toBe(30);
    c.transition.tick(0);
    c.setState({ radius: 10, position: Vec3.create(0, 0, 40) }, 200);
    c.transition.tick(100);
    expect(c.getSnapshot().radius).not.toBe(10);
    expect(c.getSettledSnapshot().radius).toBe(10);
    expect(Vec3.equals(c.getSettledSnapshot().position, Vec3.create(0, 0, 40))).toBe(true);
  });
});
