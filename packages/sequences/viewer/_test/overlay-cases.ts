/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The overlay's hover cases, for their heap guard and its quick twin (overlay.heap.ts, overlay.quick.ts): the pointer
// moving over 2 × 200,000 positions zoomed out (`scale`) and over an alignment's rows (`rows`), each with a selection
// and a highlight, an overlay frame a move (`hoverFrames`, overlay-fixtures.ts), a fresh viewer for each
// (`vi.resetModules`).

import { vi } from 'vitest';

/** The cases. */
export const KINDS = ['scale', 'rows'] as const;

/** The hover frames over `kind`: the frame, and how many moves changed the hovered element. */
export async function hovering(kind: (typeof KINDS)[number]) {
  vi.resetModules();
  const { hoverFrames } = await import('./overlay-fixtures.js');
  return hoverFrames(kind);
}
