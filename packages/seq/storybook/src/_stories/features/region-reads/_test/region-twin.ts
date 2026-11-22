/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DarkTokens, LightTokens, recordDrawLists } from '@molstar/seqx-core/render';
import { RegionReadsViews } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';

// What the region reads' twins share: a view of the stories opened on the fake clock as the story composes it (the
// generated source answering after `delayMs`, its coarse levels a third of it), then moved by the test; and the pass's
// primitives and the SVG export in each theme.

/** The delay of the finest level's answers: the coarse ones come after a third of it. */
export const DelayMs = 300;

export const view = (id: string) => RegionReadsViews.find((v) => v.id === id)!;

/** The view opened with the clock stopped after the load (`wait` ms), the story's opening not yet run. */
export function openRegion(id: string, wait = 0): Promise<SeqViewer> {
  return openExample(view(id), { delayMs: DelayMs, wait });
}

/** Runs the clock `ms` and draws. */
export async function advance(v: SeqViewer, ms: number): Promise<void> {
  await vi.advanceTimersByTimeAsync(ms);
  v.draw();
}

/** The view's SVG export in a theme, once its tiles are in. */
export async function exported(v: SeqViewer, dark: boolean) {
  const r = v.renderSvg({ idPrefix: dark ? 'rd-' : 'rl-', tokens: dark ? DarkTokens : LightTokens });
  await vi.advanceTimersByTimeAsync(2000);
  return r;
}

/** The data pass's primitives of one kind drawn, by their first word. */
export function drawn(v: SeqViewer, kind: string, pass: 'data' | 'labels' | 'hover' = 'data'): number {
  return recordDrawLists(v.frame!.lists, [pass])
    .split('\n')
    .filter((l) => l.startsWith(`${kind} `)).length;
}
