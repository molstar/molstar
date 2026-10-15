/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The streamed genome's cases, for their heap guard and its quick twin (streams-port.heap.ts, streams-port.quick.ts):
// over the streamed hg38-sized genome with its regions held, a pan's pass, and a pan's pass and composer frame (the
// values drawn from the regions, `StreamView.values` per device column), at the whole view and at 10 kb; and the same
// frame with dynamic nodes held (the manifest's genes drawn, every data region's peak). Each call loads the viewer and
// its frame afresh (`vi.resetModules`).

import { vi } from 'vitest';

/** The cases: name, the view, whether the frame is drawn (else the pass only) and whether dynamic nodes are held. */
export const PORT_CASES = [
  { name: 'a pass, whole', at: 'whole', draw: false, dynamic: false },
  { name: 'a pass and a frame, whole', at: 'whole', draw: true, dynamic: false },
  { name: 'a pass, 10 kb', at: '10 kb', draw: false, dynamic: false },
  { name: 'a pass and a frame, 10 kb', at: '10 kb', draw: true, dynamic: false },
  // ST1b-a: with dynamic nodes held, the same.
  { name: 'a pass and a frame, whole, dynamic nodes held', at: 'whole', draw: true, dynamic: true },
  { name: 'a pass and a frame, 10 kb, dynamic nodes held', at: '10 kb', draw: true, dynamic: true },
] as const;

export type PortCase = (typeof PORT_CASES)[number];

/** The genome composed with every request answered: `frame` is the step. */
export async function portFrame(at: PortCase['at'], draw = true, dynamic = false) {
  vi.resetModules();
  const { genomeFrame } = await import('./stream-fixtures.js');
  return genomeFrame(at, draw, dynamic);
}
