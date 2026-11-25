/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { expect } from 'vitest';
import { turnUntil } from '../../../examples/_test/draw-checks.js';

// What the genome ring's headless tests share.

/**
 * The recorded cytobands' rows in view: their pages inflate on the platform's thread pool (no clock), so the clock is
 * turned and real I/O let through until a bands layer has rows (or a bound of turns).
 */
export async function bandsLanded(v: SeqViewer) {
  const rows = () =>
    v
      .frame!.layout.sections[0].lays.flat()
      .some((l) => l.key.endsWith('/bands/bands') && (l.lay as { rows?: { n: number } }).rows!.n > 0);
  expect(await turnUntil(rows, 200, () => v.draw())).toBe(true);
}
