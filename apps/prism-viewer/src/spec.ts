/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PrismViewerSpec } from '@molstar/prism-viewer/spec';
import { FullPrismViewerSpec } from './full-spec.js';
import { optionsFromUrl } from './url.js';

/** The spec the app ships (viewer-spec.md §5): the viewer's defaults, Mol*'s extensions and the URL's options. */
export function createAppSpec(search: string): PrismViewerSpec {
  const base = FullPrismViewerSpec();
  return {
    ...base,
    options: optionsFromUrl(search),
    // The assistant, its own chunk, loaded after the viewer is up (viewer-spec.md §2.4).
    agent: { ...(base.agent || {}), assistant: () => import('@molstar/prism-assistant/viewer') },
  };
}
