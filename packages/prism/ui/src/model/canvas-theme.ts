/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, map, type Observable } from 'rxjs';

/**
 * The theme of what is drawn over the canvas (toolbar, cards, labels), chosen by the canvas background, not the UI
 * theme, so it stays visible: `dark` chrome (light text) on a dark background, `light` chrome on a light one.
 */
export type CanvasTheme = 'light' | 'dark';

/** The relative luminance (WCAG, 0..1) of a `0xRRGGBB` color. */
export function luminance(color: number): number {
  const channel = (c: number) => {
    const v = c / 255;
    return v <= 0.04045 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4;
  };
  return (
    0.2126 * channel((color >> 16) & 0xff) + 0.7152 * channel((color >> 8) & 0xff) + 0.0722 * channel(color & 0xff)
  );
}

/** 0.179: where light and dark text have the same contrast against the background. */
const Threshold = 0.179;

export function themeForLuminance(l: number): CanvasTheme {
  return l < Threshold ? 'dark' : 'light';
}

export function themeForColor(color: number): CanvasTheme {
  return themeForLuminance(luminance(color));
}

/** The parts of the Canvas3D props that decide what the background looks like. */
export interface CanvasBackgroundProps {
  renderer: { backgroundColor: number };
  transparentBackground?: boolean;
  postprocessing?: {
    background?: {
      variant:
        | { name: 'off' }
        | { name: 'horizontalGradient'; params: { topColor: number; bottomColor: number; opacity?: number } }
        | { name: 'radialGradient'; params: { centerColor: number; edgeColor: number; opacity?: number } }
        | { name: 'image'; params?: unknown }
        | { name: 'skybox'; params?: unknown };
    };
  };
}

/**
 * The background's luminance as the chrome sees it: the background color, with a gradient's two colors averaged
 * over it by the gradient's opacity. Images and skyboxes count as the background color for now (sampling the
 * rendered frame would be more exact). `undefined` for a transparent canvas: the host page decides.
 */
export function backgroundLuminance(props: CanvasBackgroundProps): number | undefined {
  if (props.transparentBackground) return undefined;
  const base = luminance(props.renderer.backgroundColor);
  const variant = props.postprocessing?.background?.variant;
  if (!variant || variant.name === 'off' || variant.name === 'image' || variant.name === 'skybox') return base;
  const [a, b] =
    variant.name === 'horizontalGradient'
      ? [variant.params.topColor, variant.params.bottomColor]
      : [variant.params.centerColor, variant.params.edgeColor];
  const opacity = variant.params.opacity ?? 1;
  return opacity * ((luminance(a) + luminance(b)) / 2) + (1 - opacity) * base;
}

/** The canvas theme for the given props; `fallback` when the background is transparent. */
export function canvasThemeOf(props: CanvasBackgroundProps, fallback: CanvasTheme = 'light'): CanvasTheme {
  const l = backgroundLuminance(props);
  return l === undefined ? fallback : themeForLuminance(l);
}

/** Follows the canvas props: a new theme only when it changes. */
export function canvasTheme(props: Observable<CanvasBackgroundProps>, fallback?: CanvasTheme): Observable<CanvasTheme> {
  return props.pipe(
    map((p) => canvasThemeOf(p, fallback)),
    distinctUntilChanged(),
  );
}

/** Parses `#rgb`, `#rrggbb` or `rgb()/rgba()` into `0xRRGGBB` (alpha ignored); undefined otherwise. */
export function parseColor(text: string): number | undefined {
  const t = text.trim();
  const hex = /^#?([0-9a-f]{3}|[0-9a-f]{6})$/i.exec(t)?.[1];
  if (hex) return Number.parseInt(hex.length === 3 ? [...hex].map((c) => c + c).join('') : hex, 16);
  const rgb = /^rgba?\(\s*(\d+)[\s,]+(\d+)[\s,]+(\d+)/i.exec(t);
  if (rgb) return (Number(rgb[1]) << 16) | (Number(rgb[2]) << 8) | Number(rgb[3]);
  return undefined;
}
