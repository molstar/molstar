/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { SogFile, SogMeta } from '@molstar/io/reader/sog/schema';
import {
  createShAccessor,
  type GaussianSplat,
  type GaussianSplatData,
  type GaussianSplatReadOptions,
  SH_C0,
  setCenterRotation,
} from './data.js';

/** Coefficients past the zeroth a channel, by bands (0 to 4) */
const SogShCoefficients = [0, 3, 8, 15, 24];

/** Version 2's smallest three: the indices (in `w, x, y, z`) of the stored three, by the largest component's index */
const QuatIndices = [1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2];

/** A position from the log domain: `sign(n) * (exp(|n|) - 1)` */
function fromLogDomain(n: number) {
  return Math.sign(n) * (Math.exp(Math.abs(n)) - 1);
}

/**
 * Splat `i`'s center coordinate `c` from SOG's means (`meansLow`, `meansHigh`: RGBA bytes, a pixel a splat): 16 bits an
 * axis within `meta.means`' range, out of the log domain
 */
export function createSogMeansReader(meta: SogMeta, meansLow: Uint8Array, meansHigh: Uint8Array) {
  const { mins, maxs } = meta.means;
  return (i: number, c: number) => {
    const o = i * 4 + c;
    const v = meansLow[o] | (meansHigh[o] << 8);
    return fromLogDomain(mins[c] + (maxs[c] - mins[c]) * (v / 65535));
  };
}

/**
 * The layout of SOG's palette of view-dependent color (`meta.shN`) with a centroids image of `width` × `height`: the
 * coefficients an entry (pixels), its bands (0 to 3: none or unknown 0) and its entries (version 2's `meta.json` gives
 * neither `count` nor `bands`: from the image, 64 entries a row)
 */
export function getSogShLayout(meta: SogMeta, width: number, height: number) {
  const shN = meta.shN!;
  const coefficients = shN.bands !== undefined ? SogShCoefficients[shN.bands] : width / 64;
  const bands = Math.max(0, SogShCoefficients.indexOf(coefficients));
  return { coefficients, bands: Math.min(3, bands), count: shN.count ?? height * 64 };
}

/**
 * The gaussian splats of PlayCanvas's SOG (versions 2 and 3, see `SogFile`): in the frame of the 3DGS PLY (y-down).
 * Version 3's opacities above 1 (to 1.5625) are clamped to 1
 */
export function gaussianSplatsFromSog(file: SogFile, options: GaussianSplatReadOptions): GaussianSplatData {
  const { meta, version } = file;
  const scaleCodes = meta.scales.codebook;
  const colorCodes = meta.sh0.codebook;
  const scales = file.scales.data;
  const quats = file.quats.data;
  const sh0 = file.sh0.data;
  const position = createSogMeansReader(meta, file.meansLow.data, file.meansHigh.data);
  // `w, x, y, z`
  const q = [1, 0, 0, 0];
  const rotation =
    version === 3
      ? (o: number) => {
          // the exponential map: the axis times the angle over pi
          const px = (quats[o] - 128) / 127;
          const py = (quats[o + 1] - 128) / 127;
          const pz = (quats[o + 2] - 128) / 127;
          const r = Math.sqrt(px * px + py * py + pz * pz);
          const t = (Math.PI / 2) * r;
          const k = r > 1e-12 ? Math.sin(t) / r : Math.PI / 2;
          q[0] = Math.cos(t);
          q[1] = px * k;
          q[2] = py * k;
          q[3] = pz * k;
        }
      : (o: number) => {
          const largest = quats[o + 3] - 252;
          if (largest < 0 || largest > 3) {
            q[0] = 1;
            q[1] = q[2] = q[3] = 0;
            return;
          }
          const a = ((quats[o] / 255) * 2 - 1) / Math.SQRT2;
          const b = ((quats[o + 1] / 255) * 2 - 1) / Math.SQRT2;
          const c = ((quats[o + 2] / 255) * 2 - 1) / Math.SQRT2;
          const base = largest * 3;
          q[QuatIndices[base]] = a;
          q[QuatIndices[base + 1]] = b;
          q[QuatIndices[base + 2]] = c;
          q[largest] = Math.sqrt(Math.max(0, 1 - (a * a + b * b + c * c)));
        };
  const opacity =
    version === 3
      ? (a: number) => {
          const s = a / 204;
          // to 1.5625 in the file: clamped (the opacity's encoding is 0 to 1)
          return Math.min(1, s * s);
        }
      : (a: number) => a / 255;

  let sh: GaussianSplatData['sh'];
  if (meta.shN && file.shN && options.shDegree > 0) {
    const { centroids, labels } = file.shN;
    const { coefficients, bands, count: paletteCount } = getSogShLayout(meta, centroids.width, centroids.height);
    const codes = meta.shN.codebook;
    const width = centroids.width;
    const c = centroids.data;
    const l = labels.data;
    sh = createShAccessor(bands > 0 ? Math.min(options.shDegree, bands) : 0, (i, k, channel) => {
      const o = i * 4;
      const label = l[o] | (l[o + 1] << 8) | (version === 3 ? l[o + 2] << 16 : 0);
      // a label past the palette: none (as PlayCanvas reads it)
      if (label >= paletteCount) return 0;
      const u = (label % 64) * coefficients + k;
      const v = Math.floor(label / 64);
      return codes[c[(v * width + u) * 4 + channel]];
    });
  }

  return {
    count: file.count,
    upright: false,
    getSplat: (i: number, out: GaussianSplat) => {
      const o = i * 4;
      rotation(o);
      setCenterRotation(out, position(i, 0), position(i, 1), position(i, 2), q[1], q[2], q[3], q[0]);
      out.logScaleX = scaleCodes[scales[o]];
      out.logScaleY = scaleCodes[scales[o + 1]];
      out.logScaleZ = scaleCodes[scales[o + 2]];
      out.r = 0.5 + SH_C0 * colorCodes[sh0[o]];
      out.g = 0.5 + SH_C0 * colorCodes[sh0[o + 1]];
      out.b = 0.5 + SH_C0 * colorCodes[sh0[o + 2]];
      out.opacity = opacity(sh0[o + 3]);
    },
    sh,
  };
}
