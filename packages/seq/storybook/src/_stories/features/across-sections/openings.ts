/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { GrapeSection, PeachSection, ZoomedGrapeLocus, ZoomedLocus } from '@molstar/seqx-examples/across/views';
import { goToLocus } from '@molstar/seqx-viewer/locus';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// Where the across-sections stories open, run once the view is placed (`exampleStory`'s `opening`).

/** Peach zoomed into chr2 from 1.2 to 1.5 Mb, grape left at the fit (Defaults, 3). */
export const openZoomed = (viewer: SeqViewer) => goToLocus(viewer, ZoomedLocus, { section: PeachSection });

/** Grape zoomed into chr2 from 1.2 to 1.5 Mb, peach left at the fit (the dense alignment: peach on bins). */
export const openZoomedGrape = (viewer: SeqViewer) => goToLocus(viewer, ZoomedGrapeLocus, { section: GrapeSection });
