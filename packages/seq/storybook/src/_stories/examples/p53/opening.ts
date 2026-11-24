/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// Gene by region opens on a region of the gene, not on all of it: the viewer's first move is a jump (the toolbar's Back
// goes to the whole gene, then Forward comes back), and the regions beyond it are read as the view moves to them.

/** The region the story opens on: the gene's first exons and their transcripts. */
export const OpeningAddress = '/data/tp53-region/(range :from 2000 :to 3500)';

/** Moves the viewer to the opening region (once the view is placed). */
export async function openOnRegion(viewer: SeqViewer): Promise<void> {
  await viewer.runCommand('address.focus', { address: OpeningAddress });
}
