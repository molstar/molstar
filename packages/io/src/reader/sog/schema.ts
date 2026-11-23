/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { DecodedImage } from '../../common/image.js';

/**
 * PlayCanvas's SOG gaussian splats, versions 2 and 3 (`meta.json`, see `SogFile`). Unknown fields ignored
 */
export interface SogMeta {
  /** 2 or 3 (version 1 has none) */
  readonly version?: number;
  readonly count: number;
  /** positions by axis within the ranges, in the log domain (`sign(n) * (exp(|n|) - 1)`) */
  readonly means: { readonly mins: number[]; readonly maxs: number[]; readonly files: string[] };
  /** natural logs of the scales by byte */
  readonly scales: { readonly codebook: number[]; readonly files: string[] };
  readonly quats: { readonly files: string[] };
  /** the zeroth spherical harmonic by byte (rgb), the opacity in alpha */
  readonly sh0: { readonly codebook: number[]; readonly files: string[] };
  /**
   * the spherical harmonics past the zeroth: a palette of `count` entries, a label a splat (version 2 gives neither
   * `count` nor `bands`: from the palette image, 64 entries a row of `coefficients` pixels each)
   */
  readonly shN?: {
    readonly count?: number;
    readonly bands?: number;
    readonly codebook: number[];
    readonly files: string[];
  };
  /** version 3, where an opacity exceeds 1 */
  readonly opacity?: { readonly max: number };
  /** version 3: 'antialiased' or '2dgs', none for ordinary gaussians */
  readonly model?: string;
}

/** A decoded image of a SOG file: RGBA bytes, row-major from the top-left, not premultiplied (`decodeImage`) */
export type SogImage = DecodedImage;

/** An image of a SOG file as stored: its name and lossless WebP bytes (a view of the archive's entry) */
export interface SogImageBytes {
  readonly name: string;
  readonly bytes: Uint8Array;
}

/**
 * A SOG file's archive: its `meta.json` and its images as stored (not decoded: the graphics may upload them as they
 * are, `SogFile` for the decoded images)
 */
export interface SogArchive {
  readonly meta: SogMeta;
  readonly version: 2 | 3;
  readonly count: number;
  readonly meansLow: SogImageBytes;
  readonly meansHigh: SogImageBytes;
  readonly scales: SogImageBytes;
  readonly quats: SogImageBytes;
  readonly sh0: SogImageBytes;
  readonly shN?: { readonly centroids: SogImageBytes; readonly labels: SogImageBytes };
}

/** Decodes the lossless WebP `bytes` of the file `name` */
export type SogImageDecoder = (bytes: Uint8Array, name: string) => Promise<SogImage>;

/**
 * A SOG file: its `meta.json` and its images decoded (the only arrays made: 4 bytes a pixel), a splat a pixel in
 * row-major order (splat `i` at `4 i` of each but the palette):
 * - `means`: the low and high bytes of 16 bits an axis (rgb), a fraction of `meta.means`' range
 * - `scales`: rgb bytes into `meta.scales.codebook`
 * - `quats`: version 3 an exponential map (rgb `(b - 128) / 127`, `w, x, y, z`); version 2 the smallest three (rgb
 *   `(b / 255 * 2 - 1) / sqrt 2`, alpha 252 + the largest component's index in `w, x, y, z`)
 * - `sh0`: rgb bytes into `meta.sh0.codebook`, the opacity in alpha (version 3 `(a / 204)²`, version 2 `a / 255`)
 * - `shN`: the palette (`centroids`: 64 entries a row, an entry `coefficients` pixels, rgb bytes into
 *   `meta.shN.codebook`) and a label a splat (`labels`: rg, version 3 rgb, little-endian)
 *
 * In the frame of the 3DGS PLY (y-down)
 */
export interface SogFile {
  readonly meta: SogMeta;
  readonly version: 2 | 3;
  readonly count: number;
  readonly meansLow: SogImage;
  readonly meansHigh: SogImage;
  readonly scales: SogImage;
  readonly quats: SogImage;
  readonly sh0: SogImage;
  readonly shN?: { readonly centroids: SogImage; readonly labels: SogImage };
}
