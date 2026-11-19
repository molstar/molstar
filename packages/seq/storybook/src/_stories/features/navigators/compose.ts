/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// The navigators' stories' composition, DOM-free: the Strip story and its headless twin call these.

/** Where the Strip story opens: the DNA-binding domain's first half (a window smaller than the strip: it can be dragged). */
export const StripOpening = '/data/p53/(range :from 94 :to 200)';

/** Moves the viewer to where the Strip story opens. */
export function openStrip(viewer: SeqViewer): void {
  viewer.runCommand('address.focus', { address: StripOpening });
}
