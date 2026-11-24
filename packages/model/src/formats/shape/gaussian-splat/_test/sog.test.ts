/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Quat } from '@molstar/core/math/linear-algebra';
import type { SogFile, SogImage, SogMeta } from '@molstar/io/reader/sog/schema';
import { SH_C0 } from '../data.js';
import { gaussianSplatsFromSog } from '../sog.js';
import { decode, expectRotation } from './decode.js';

interface TestSplat {
  position: number[];
  /** log scales, codebook entries */
  logScale: number[];
  rotation: Quat;
  opacity: number;
  /** the zeroth harmonic, codebook entries */
  dc: number[];
}

const Splats: TestSplat[] = [
  {
    position: [1.5, -2, 0.25],
    logScale: [-1, -2, -3],
    rotation: Quat.create(0.2, -0.5, 0.1, 0.8),
    opacity: 0.75,
    dc: [1, 0, -1],
  },
  {
    position: [-3, 0, 4],
    logScale: [0, -2, 1],
    rotation: Quat.create(-0.9, 0.2, 0.1, 0.3),
    opacity: 0.25,
    dc: [0.5, -0.5, 0],
  },
  {
    position: [0, 2.5, -1],
    logScale: [1, 0, -1],
    rotation: Quat.create(0.1, 0.8, -0.2, -0.3),
    opacity: 1,
    dc: [0, 1, 0.5],
  },
];

/** the codebook of the values used: a byte each */
const Codebook = [-3, -2, -1, -0.5, 0, 0.5, 1];
const code = (v: number) => Codebook.indexOf(v);

function image(width: number, height: number, fill: (i: number, out: Uint8Array) => void): SogImage {
  const data = new Uint8Array(width * height * 4);
  const px = new Uint8Array(4);
  for (let i = 0; i < width * height; ++i) {
    px.fill(0);
    fill(i, px);
    data.set(px, i * 4);
  }
  return { width, height, data };
}

/** the log domain of positions: `sign(x) * ln(|x| + 1)` */
const toLog = (x: number) => Math.sign(x) * Math.log(Math.abs(x) + 1);

/** splat-transform's packing of version 3's exponential map (`w, x, y, z`, w made positive) */
function expMapBytes(q: Quat) {
  let [x, y, z, w] = Quat.normalize(Quat(), q);
  if (w < 0) [x, y, z, w] = [-x, -y, -z, -w];
  const sv = Math.sqrt(x * x + y * y + z * z);
  const k = sv > 1e-12 ? Math.atan2(sv, w) / sv / (Math.PI / 2) : 1 / (Math.PI / 2);
  return [x, y, z].map((v) => Math.min(255, Math.max(0, Math.round(v * k * 127 + 128))));
}

/** version 2's smallest three (`w, x, y, z`, the largest made positive, its index + 252 in alpha) */
function smallestThreeBytes(q: Quat) {
  const [x, y, z, w] = Quat.normalize(Quat(), q);
  const v = [w, x, y, z];
  let largest = 0;
  for (let k = 1; k < 4; ++k) if (Math.abs(v[k]) > Math.abs(v[largest])) largest = k;
  const sign = v[largest] < 0 ? -1 : 1;
  const others = [0, 1, 2, 3].filter((k) => k !== largest);
  return [...others.map((k) => Math.round(((v[k] * sign * Math.SQRT2 + 1) / 2) * 255)), 252 + largest];
}

/** a SOG of `Splats` in images of 2x2 (one pixel unused), version 2 or 3 */
function sogFile(version: 2 | 3, sh?: { bands: 1; palette: number[][][]; labels: number[] }): SogFile {
  const mins = [-3, -3, -3].map(toLog);
  const maxs = [4, 4, 4].map(toLog);
  const q16 = (i: number, c: number) =>
    Math.round(((toLog(Splats[i].position[c]) - mins[c]) / (maxs[c] - mins[c])) * 65535);
  const count = Splats.length;
  const meta: SogMeta = {
    version,
    count,
    means: { mins, maxs, files: [] },
    scales: { codebook: Codebook, files: [] },
    quats: { files: [] },
    sh0: { codebook: Codebook, files: [] },
    shN: sh
      ? { codebook: Codebook, files: [], ...(version === 3 ? { count: sh.palette.length, bands: 1 } : {}) }
      : undefined,
  };
  const per = (fill: (s: TestSplat, i: number, out: Uint8Array) => void) =>
    image(2, 2, (i, out) => {
      if (i < count) fill(Splats[i], i, out);
    });
  return {
    meta,
    version,
    count,
    meansLow: per((_, i, out) => out.set([0, 1, 2].map((c) => q16(i, c) & 0xff))),
    meansHigh: per((_, i, out) => out.set([0, 1, 2].map((c) => q16(i, c) >> 8))),
    scales: per((s, _, out) => out.set(s.logScale.map(code))),
    quats: per((s, _, out) => out.set(version === 3 ? expMapBytes(s.rotation) : smallestThreeBytes(s.rotation))),
    sh0: per((s, _, out) => {
      out.set(s.dc.map(code));
      out[3] = version === 3 ? Math.round(Math.sqrt(s.opacity) * 204) : Math.round(s.opacity * 255);
    }),
    shN: sh
      ? {
          // 64 entries a row, 3 coefficients (pixels) each
          centroids: image(64 * 3, 1, (p, out) => {
            const entry = sh.palette[Math.floor(p / 3)];
            if (entry) out.set(entry[p % 3].map(code));
          }),
          labels: per((_, i, out) => out.set([sh.labels[i] & 0xff, sh.labels[i] >> 8, 0])),
        }
      : undefined,
  };
}

describe('gaussian splats of SOG files', () => {
  for (const version of [2, 3] as const) {
    it(`reads version ${version}: positions in the log domain, codebook scales and colors, its rotations`, () => {
      const data = gaussianSplatsFromSog(sogFile(version), { shDegree: 0 });
      expect(data.count).toBe(3);
      expect(data.upright).toBe(false);
      decode(data).forEach((d, i) => {
        const s = Splats[i];
        expect(d.x).toBeCloseTo(s.position[0], 3);
        expect(d.y).toBeCloseTo(s.position[1], 3);
        expect(d.z).toBeCloseTo(s.position[2], 3);
        expect([d.logScaleX, d.logScaleY, d.logScaleZ]).toEqual(s.logScale);
        expectRotation(d, s.rotation, 2);
        expect(d.opacity).toBeCloseTo(s.opacity, 2);
        expect([d.r, d.g, d.b]).toEqual(s.dc.map((v) => 0.5 + SH_C0 * v));
      });
    });
  }

  it("clamps version 3's opacities above 1 (to 1.5625)", () => {
    const file = sogFile(3);
    file.sh0.data[3] = 255;
    expect(decode(gaussianSplatsFromSog(file, { shDegree: 0 }))[0].opacity).toBe(1);
  });

  it('reads the view-dependent color from the palette by label, its bands from the palette in version 2', () => {
    // 2 entries of 3 coefficients (rgb codebook values each), the splats' labels
    const palette = [
      [
        [1, 0, -1],
        [0.5, -0.5, 0],
        [-3, -2, 1],
      ],
      [
        [0, 0, 0],
        [-1, 1, 0.5],
        [0.5, 0.5, -0.5],
      ],
    ];
    const labels = [1, 0, 1];
    for (const version of [2, 3] as const) {
      const file = sogFile(version, { bands: 1, palette, labels });
      expect(gaussianSplatsFromSog(file, { shDegree: 0 }).sh).toBeUndefined();
      const { sh } = gaussianSplatsFromSog(file, { shDegree: 3 });
      expect(sh?.degree).toBe(1);
      for (let i = 0; i < 3; ++i) {
        for (let k = 0; k < 3; ++k) {
          for (let c = 0; c < 3; ++c) expect(sh!.get(i, k, c)).toBe(palette[labels[i]][k][c]);
        }
      }
    }
  });
});
