/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Fluid } from '../fluid.js';

function run(fluid: Fluid, seconds: number) {
  for (let t = 0; t < seconds; t += 1 / 60) fluid.step(1 / 60);
}

function at(fluid: Fluid, x: number, y: number) {
  fluid.sample(x, y);
  return { u: fluid.su, v: fluid.sv };
}

describe('Fluid', () => {
  it('is off over nothing, and still air stays still', () => {
    const none = new Fluid();
    none.layout(0, 0, 0, 0);
    expect(none.on).toBe(false);
    none.step(1 / 60);
    expect(at(none, 10, 10)).toEqual({ u: 0, v: 0 });

    const fluid = new Fluid();
    fluid.layout(0, 0, 400, 300);
    expect(fluid.on).toBe(true);
    run(fluid, 0.5);
    expect(at(fluid, 200, 150)).toEqual({ u: 0, v: 0 });
  });

  it('carries a stroke on, and lets it die down', () => {
    const fluid = new Fluid();
    fluid.layout(0, 0, 600, 400);
    for (let k = 0; k < 10; k++) {
      fluid.drive(150 + k * 20, 200, 400, 0, 60, 0.3);
      fluid.step(1 / 60);
    }
    run(fluid, 0.2);
    // Ahead of where the stroke ended, the air moves its way.
    const ahead = at(fluid, 380, 200);
    expect(ahead.u).toBeGreaterThan(20);
    expect(Math.abs(ahead.v)).toBeLessThan(ahead.u);
    run(fluid, 8);
    expect(Math.abs(at(fluid, 380, 200).u)).toBeLessThan(ahead.u * 0.1);
  });

  it('spins round a swirl (clockwise on the screen for a positive one)', () => {
    const fluid = new Fluid();
    fluid.layout(0, 0, 400, 400);
    fluid.swirl(200, 200, 100, 120);
    // Right of it the air goes down, left of it up, above it right.
    expect(at(fluid, 240, 200).v).toBeGreaterThan(30);
    expect(at(fluid, 160, 200).v).toBeLessThan(-30);
    expect(at(fluid, 200, 160).u).toBeGreaterThan(30);
    // A swirl has no divergence to take out: the projection keeps it.
    const before = at(fluid, 240, 200).v;
    fluid.step(1 / 60);
    expect(at(fluid, 240, 200).v).toBeGreaterThan(before * 0.7);
  });

  it('takes out a push that only piles air up (it has nowhere to go)', () => {
    const fluid = new Fluid();
    fluid.layout(0, 0, 400, 400);
    // Air pushed in from both sides toward the middle.
    fluid.drive(170, 200, 200, 0, 40, 1);
    fluid.drive(230, 200, -200, 0, 40, 1);
    const before = at(fluid, 170, 200).u;
    fluid.step(1 / 60);
    expect(at(fluid, 170, 200).u).toBeLessThan(before * 0.6);
  });
});
