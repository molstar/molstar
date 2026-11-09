/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Camera } from '../camera.js';
import { formatLength, niceValue, worldPerPixel } from '../helper/scale-line.js';

describe('scale line', () => {
  it('is the largest of 1, 2, 5 × 10ⁿ that fits', () => {
    expect(niceValue(1)).toBe(1);
    expect(niceValue(1.99)).toBe(1);
    expect(niceValue(2)).toBe(2);
    expect(niceValue(4.9)).toBe(2);
    expect(niceValue(5)).toBe(5);
    expect(niceValue(99)).toBe(50);
    expect(niceValue(120)).toBe(100);
    expect(niceValue(0.37)).toBeCloseTo(0.2, 12);
  });

  it('reads in Å, nm or µm', () => {
    expect(formatLength(0.2)).toBe('0.2 Å');
    expect(formatLength(50)).toBe('50 Å');
    expect(formatLength(500)).toBe('500 Å');
    expect(formatLength(1000)).toBe('100 nm');
    expect(formatLength(5000)).toBe('500 nm');
    expect(formatLength(20000)).toBe('2 µm');
    // as written, without float noise
    expect(formatLength(0.1 + 0.2)).toBe('0.3 Å');
  });

  // a 200 CSS pixel high view (400 device pixels at a ratio of 2), 100 Å from its target at 90°: 200 Å high
  it('measures a CSS pixel at the camera target', () => {
    const camera = new Camera(
      { position: Vec3.create(0, 0, 100), target: Vec3(), fov: Math.PI / 2 },
      { x: 0, y: 0, width: 400, height: 400 },
    );
    expect(worldPerPixel(camera, 2)).toBeCloseTo(1, 6);
    // twice the distance: twice the size
    camera.setState({ position: Vec3.create(0, 0, 200) }, 0);
    expect(worldPerPixel(camera, 2)).toBeCloseTo(2, 6);
    // a cropped part of the frame: the field of view still spans the frame
    Camera.setCrop(camera.crop, camera.viewport.width, camera.viewport.height, 0, 0);
    Object.assign(camera.viewport, { height: camera.viewport.height / 2 });
    expect(worldPerPixel(camera, 2)).toBeCloseTo(2, 6);
  });
});
