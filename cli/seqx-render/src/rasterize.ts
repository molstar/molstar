/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { fileURLToPath } from 'node:url';
import type { Rasterizer } from '@molstar/seqx-core/render';
import { Resvg } from '@resvg/resvg-js';

// The render CLI's rasterizer port (the core's `Rasterizer`): resvg over the bundled fonts alone. The text metrics the
// layouts measure with are Arial's advances for the sans role and a fixed 600 advance for the mono role; Liberation
// Sans and Liberation Mono (SIL OFL 1.1, `fonts/OFL.txt`) are metric-compatible with Arial and Courier New, so the
// labels the layout fitted fit as drawn. System fonts are never loaded: the same SVG draws the same on every machine.

const Dir = fileURLToPath(new URL('../fonts/', import.meta.url));
const Families = { sans: 'Liberation Sans', mono: 'Liberation Mono' };

/** The four bundled font files (regular and bold of each family). */
export const FontFiles: readonly string[] = [
  'LiberationSans-Regular.ttf',
  'LiberationSans-Bold.ttf',
  'LiberationMono-Regular.ttf',
  'LiberationMono-Bold.ttf',
].map((f) => `${Dir}${f}`);

const MonoRole = /mono|menlo|consolas|courier/i;

/**
 * The SVG with each `font-family` resolved to a bundled family: a stack naming a mono face (the token table's mono role)
 * to Liberation Mono, any other to Liberation Sans. Why: resvg 2.6.2 drops a stack to its default family when an
 * unknown family comes before a bundled one, and `monospace` ignores `monospaceFamily`, so the stacks cannot be left
 * for it to resolve (the letters would draw in the sans).
 */
export function resolveFonts(svg: string): string {
  return svg.replace(
    /font-family="([^"]*)"/g,
    (_, stack: string) => `font-family="${MonoRole.test(stack) ? Families.mono : Families.sans}"`,
  );
}

/** An SVG string as a PNG through resvg: zoomed by the request's scale, bundled fonts only. */
export const resvgRasterizer: Rasterizer = (svg, request) => {
  const resvg = new Resvg(resolveFonts(svg), {
    fitTo: { mode: 'zoom', value: request.scale },
    font: {
      loadSystemFonts: false,
      fontFiles: [...FontFiles],
      defaultFontFamily: Families.sans,
      sansSerifFamily: Families.sans,
      serifFamily: Families.sans,
      monospaceFamily: Families.mono,
    },
  });
  return new Uint8Array(resvg.render().asPng());
};
