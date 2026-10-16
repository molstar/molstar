/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { fft2 } from '../quantum-field.js';

const N = 128;

describe('the quantum splash’s FFT', () => {
  it('comes back to where it started (forward, then inverse over N²)', () => {
    const re = new Float32Array(N * N);
    const im = new Float32Array(N * N);
    let s = 7;
    for (let i = 0; i < re.length; i++) {
      s = (s * 16807) % 2147483647;
      re[i] = (s / 2147483647) * 2 - 1;
      im[i] = Math.sin(i * 0.37);
    }
    const re0 = re.slice();
    const im0 = im.slice();
    fft2(re, im, false);
    fft2(re, im, true);
    let worst = 0;
    for (let i = 0; i < re.length; i++) {
      worst = Math.max(worst, Math.abs(re[i] / (N * N) - re0[i]), Math.abs(im[i] / (N * N) - im0[i]));
    }
    expect(worst).toBeLessThan(1e-4);
  });

  it('puts a plane wave in one bin: its wavenumber', () => {
    const re = new Float32Array(N * N);
    const im = new Float32Array(N * N);
    const kx = 5;
    const ky = 3;
    for (let y = 0; y < N; y++) {
      for (let x = 0; x < N; x++) {
        const phase = (2 * Math.PI * (kx * x + ky * y)) / N;
        re[y * N + x] = Math.cos(phase);
        im[y * N + x] = Math.sin(phase);
      }
    }
    fft2(re, im, false);
    const at = ky * N + kx;
    expect(re[at]).toBeCloseTo(N * N, -1);
    let rest = 0;
    for (let i = 0; i < re.length; i++) if (i !== at) rest = Math.max(rest, Math.hypot(re[i], im[i]));
    expect(rest).toBeLessThan(1);
  });
});
