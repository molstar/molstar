/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// Where the region reads' stories open, and what their twins move to: a window as an address, run once the view is
// placed (`exampleStory`'s `opening`).

/** The Hi-C chromosome's window the story opens on (positions 100 to 250: its finest level's tiles, a few). */
export const HicWindow: readonly [number, number] = [100, 250];
/** The genotypes' window round the long deletion (200,000 to 260,000), items read at their own level. */
export const DeletionWindow: readonly [number, number] = [160_000, 300_000];
/** A stretch of the genotypes' chromosome with no link end in it (the links are in the first 200 kb). */
export const EmptyWindow: readonly [number, number] = [400_000, 700_000];
/** A stretch of the items-only genotypes small enough to read (about a hundred tiles). */
export const ReadableWindow: readonly [number, number] = [200_000, 220_000];

/** Moves the view to `[from, to)` of `system`. */
export function go(viewer: SeqViewer, system: string, [from, to]: readonly [number, number]): void {
  viewer.runCommand('address.focus', { address: `/data/${system}/(range :from ${from} :to ${to})` });
}

export const openHic = (viewer: SeqViewer) => go(viewer, 'chr1', HicWindow);
export const openDeletion = (viewer: SeqViewer) => go(viewer, 'chr', DeletionWindow);
export const openEmpty = (viewer: SeqViewer) => go(viewer, 'chr', EmptyWindow);
