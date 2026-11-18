/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Size entry 3, the viewer in its default configuration. Bundled by `sizes.mjs`, never run.

import { createSeqViewer } from '@molstar/seqx-viewer';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';

export const create = (target: HTMLElement) => createSeqViewer(DefaultSeqViewerSpec(), target);
