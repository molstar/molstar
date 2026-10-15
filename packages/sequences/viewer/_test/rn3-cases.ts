/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The cases of RN3's hot paths, for their heap guard and its quick twin (rn3.heap.ts, rn3.quick.ts): the axis view's map
// synced through a pan (the window mapped through the blocks), the navigator's axis and its window with the channels'
// marks per lane, the labels axis, the grips; the rows' tiles painted with their images made (coarse under fine, the
// overview's pending shade), the rows' own draws, and a pan's composer frames with the maps' images kept. A case is a
// name, a stable id and an `open` giving its frame; the frames are rn3-fixtures.ts's. Each `open` loads the viewer and
// its frame afresh (`vi.resetModules`), as draws.heap.ts does (RN1): a frame built here, shared by the cases, read two
// copies' double fields through one load site and boxed them (45 bytes a frame from the second case on).

import { vi } from 'vitest';

export interface Rn3Case {
  /** The case's name in the guard and the twin. */
  name: string;
  /** The case's stable id, without the file's prefix. */
  id: string;
  /** The case's frame, on a viewer loaded afresh. */
  open(): Promise<() => void>;
}

const fresh = async <T>(load: () => Promise<T>): Promise<T> => {
  vi.resetModules();
  return load();
};

export const RN3_CASES: Rn3Case[] = [];
for (const stacked of [false, true])
  for (const labels of [false, true])
    RN3_CASES.push({
      name: `the navigator, ${stacked ? 'stacked in lanes' : 'side by side'}, ${labels ? 'labels' : 'hairline'}`,
      id: `navigator, ${stacked ? 'stacked' : 'side by side'}, ${labels ? 'labels' : 'hairline'}`,
      open: () => fresh(async () => (await import('./rn3-fixtures.js')).navigatorFrame(stacked, labels)),
    });
for (const remote of [false, true])
  RN3_CASES.push({
    name: `the rows' tiles, ${remote ? 'from the region cache' : 'held rows'}`,
    id: `rows' tiles, ${remote ? 'region cache' : 'held rows'}`,
    open: () => fresh(async () => (await import('./rn3-fixtures.js')).rowsFrame(remote)),
  });
for (const kind of ['single', 'dense', 'remote'] as const)
  RN3_CASES.push({
    name: `the rows' own draws, ${kind}: cells, gaps, letters, gutters, the tree and the data columns`,
    id: `rows' own draws, ${kind}`,
    open: () => fresh(async () => (await import('./rn3-fixtures.js')).rowsDrawFrame(kind)),
  });
for (const bare of [false, true])
  RN3_CASES.push({
    name: `a pan's composer frames, ${bare ? 'no track in miniature' : "seqv:map and the navigator's map"}: images kept`,
    id: `pan's composer frames, ${bare ? 'no miniature' : 'maps'}`,
    open: () => fresh(async () => (await (await import('./rn3-fixtures.js')).mapFrame(bare)).frame),
  });

/** The positive control's case (a heap guard's), the rows' tiles of held rows. */
export const rowsHeld = () => fresh(async () => (await import('./rn3-fixtures.js')).rowsFrame(false));
