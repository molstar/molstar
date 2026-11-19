/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import type { RenderSvgOptions, RenderSvgResult } from './export.js';
import { renderSvg } from './render-svg.js';

// Headless PNG (_spec/README.md "Render", Backends): the SVG export rasterized, so a PNG is the SVG's picture. The
// rasterizer is a port (`Rasterizer`): this package stays free of native code, the application composes one (the
// render CLI's, over resvg).

/** What a rasterizer is asked for: an SVG string at a size in CSS px and the scale it is drawn at. */
export interface RasterizeRequest {
  /** The SVG's own size, CSS px. */
  readonly width: number;
  readonly height: number;
  /** Output pixels per CSS px (`scale` × `dpr`). */
  readonly scale: number;
}

/** An SVG string as a PNG file's bytes, `round(width × scale)` by `round(height × scale)` pixels. */
export type Rasterizer = (svg: string, request: RasterizeRequest) => Uint8Array | Promise<Uint8Array>;

export interface RenderPngOptions extends RenderSvgOptions {
  readonly rasterize: Rasterizer;
  /** Output pixels per CSS px, 1 by default; the picture's pixels are `scale` × `dpr` per CSS px. */
  readonly scale?: number;
}

export interface RenderPngResult extends RenderSvgResult {
  readonly png: Uint8Array;
  /** The PNG's size in pixels. */
  readonly pixels: { readonly width: number; readonly height: number };
}

/** The size an SVG string's root declares (its CSS px). */
export function svgSize(svg: string): { width: number; height: number } {
  const root = /^<svg [^>]*>/.exec(svg)?.[0] ?? '';
  const w = /\swidth="([\d.]+)"/.exec(root)?.[1];
  const h = /\sheight="([\d.]+)"/.exec(root)?.[1];
  return { width: w ? +w : 0, height: h ? +h : 0 };
}

/**
 * A spec drawn to a PNG: `renderSvg`'s SVG (the same options, so the same picture) rasterized by the given port at
 * `scale` × `dpr` pixels per CSS px (the SVG paints its own background, so the pixels are opaque). A render with
 * errors (a blank SVG) is still rasterized: its issues say why.
 */
export async function renderPng(spec: Node, options: RenderPngOptions): Promise<RenderPngResult> {
  const r = await renderSvg(spec, options);
  const size = svgSize(r.svg);
  const scale = (options.scale ?? 1) * (options.dpr ?? 1);
  const png = await options.rasterize(r.svg, { width: size.width, height: size.height, scale });
  return {
    ...r,
    png,
    pixels: { width: Math.round(size.width * scale), height: Math.round(size.height * scale) },
  };
}
