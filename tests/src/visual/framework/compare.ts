/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export interface Tolerance {
  /** Per-pixel perceptual (YIQ) color delta in [0, 1] above which a pixel counts as different */
  threshold: number;
  /** Maximum allowed ratio of differing pixels */
  maxDiffRatio: number;
  /** Maximum allowed number of differing pixels */
  maxDiffPixels: number;
}

export const DefaultTolerance: Tolerance = {
  threshold: 0.1,
  maxDiffRatio: 0.001,
  maxDiffPixels: Infinity,
};

export interface CompareMetrics {
  width: number;
  height: number;
  sizeMismatch: boolean;
  diffPixels: number;
  diffRatio: number;
  /** Largest per-pixel YIQ delta in [0, 1] */
  maxDelta: number;
  meanDelta: number;
  /** Over RGB channels; Infinity for identical images */
  psnr: number;
}

// max squared YIQ distance between black and white
const MaxYiqDelta = 35215;

function blend(c: number, a: number) {
  return 255 + (c - 255) * a;
}

function yiqDelta(a: Uint8ClampedArray, b: Uint8ClampedArray, i: number) {
  const aa = a[i + 3] / 255,
    ab = b[i + 3] / 255;
  const r1 = blend(a[i], aa),
    g1 = blend(a[i + 1], aa),
    b1 = blend(a[i + 2], aa);
  const r2 = blend(b[i], ab),
    g2 = blend(b[i + 1], ab),
    b2 = blend(b[i + 2], ab);

  const y = (r1 - r2) * 0.29889531 + (g1 - g2) * 0.58662247 + (b1 - b2) * 0.11448223;
  const iq = (r1 - r2) * 0.59597799 - (g1 - g2) * 0.2741761 - (b1 - b2) * 0.32180189;
  const q = (r1 - r2) * 0.21147017 - (g1 - g2) * 0.52261711 + (b1 - b2) * 0.31114694;
  return Math.sqrt((0.5053 * y * y + 0.299 * iq * iq + 0.1957 * q * q) / MaxYiqDelta);
}

/**
 * Compares two RGBA images (not premultiplied). The diff image shows the
 * reference faded, sub-threshold differences in yellow and failing pixels in red.
 */
export function compareImages(
  reference: ImageData,
  target: ImageData,
  threshold: number,
): { metrics: CompareMetrics; diff: ImageData } {
  const width = Math.max(reference.width, target.width);
  const height = Math.max(reference.height, target.height);
  const diff = new ImageData(width, height);

  if (reference.width !== target.width || reference.height !== target.height) {
    diff.data.fill(255);
    return {
      metrics: {
        width,
        height,
        sizeMismatch: true,
        diffPixels: width * height,
        diffRatio: 1,
        maxDelta: 1,
        meanDelta: 1,
        psnr: 0,
      },
      diff,
    };
  }

  const a = reference.data,
    b = target.data,
    d = diff.data;
  const n = width * height;
  let diffPixels = 0,
    maxDelta = 0,
    sumDelta = 0,
    sumSq = 0;

  for (let p = 0; p < n; ++p) {
    const i = p * 4;
    const delta = yiqDelta(a, b, i);
    sumDelta += delta;
    if (delta > maxDelta) maxDelta = delta;

    for (let c = 0; c < 3; ++c) {
      const e = a[i + c] - b[i + c];
      sumSq += e * e;
    }

    if (delta > threshold) {
      diffPixels += 1;
      d[i] = 255;
      d[i + 1] = 0;
      d[i + 2] = 0;
    } else if (delta > 0) {
      d[i] = 255;
      d[i + 1] = 200;
      d[i + 2] = 0;
    } else {
      const l = 0.299 * a[i] + 0.587 * a[i + 1] + 0.114 * a[i + 2];
      const v = 255 - (255 - l) * 0.1;
      d[i] = v;
      d[i + 1] = v;
      d[i + 2] = v;
    }
    d[i + 3] = 255;
  }

  const mse = sumSq / (n * 3);
  return {
    metrics: {
      width,
      height,
      sizeMismatch: false,
      diffPixels,
      diffRatio: diffPixels / n,
      maxDelta,
      meanDelta: sumDelta / n,
      psnr: mse === 0 ? Infinity : 10 * Math.log10((255 * 255) / mse),
    },
    diff,
  };
}

export function isWithinTolerance(metrics: CompareMetrics, tolerance: Tolerance) {
  return (
    !metrics.sizeMismatch &&
    metrics.diffRatio <= tolerance.maxDiffRatio &&
    metrics.diffPixels <= tolerance.maxDiffPixels
  );
}

/** True if all pixels are identical, e.g. nothing was rendered */
export function isUniform(image: ImageData) {
  const d = new Uint32Array(image.data.buffer, image.data.byteOffset, image.width * image.height);
  const first = d[0];
  for (let i = 1, il = d.length; i < il; ++i) {
    if (d[i] !== first) return false;
  }
  return true;
}
