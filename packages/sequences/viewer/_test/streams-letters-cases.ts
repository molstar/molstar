/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The streamed letters' cases, for their heap guard and its quick twin (streams-letters.heap.ts,
// streams-letters.quick.ts): over the streamed genome with a sequence track and its regions held, a pan's pass and
// composer frame at 10 kb (the letters' cells thinner than a device px, read raw per device column from
// `StreamView.letters`) and at 100 positions (the letters drawn). Each call loads the viewer and its frame afresh
// (`vi.resetModules`).

import { vi } from 'vitest';

export const LETTER_VIEWS = ['10 kb', '100 b'] as const;

/** The genome composed with its letters streamed, at `at`, every request answered: `frame` is the step. */
export async function lettersFrame(at: (typeof LETTER_VIEWS)[number]) {
  vi.resetModules();
  const { genomeFrame } = await import('./stream-fixtures.js');
  return genomeFrame(at, true, false, true);
}
