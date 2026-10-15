/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import type { PrismViewerSpec } from '@molstar/prism-viewer/spec';
import { MolstarViewerExtensions } from './extensions.js';
import { optionsFromUrl } from './url.js';

/** The spec the app ships (viewer-spec.md §5): the viewer's defaults, Mol*'s extensions and the URL's options. */
export function createAppSpec(search: string): PrismViewerSpec {
  const base = DefaultPrismViewerSpec();
  return {
    ...base,
    extensions: [...(base.extensions ?? []), ...MolstarViewerExtensions],
    options: optionsFromUrl(search),
    // The assistant, its own chunk, loaded after the viewer is up (viewer-spec.md §2.4).
    agent: { ...(base.agent || {}), assistant: () => import('@molstar/prism-assistant/viewer') },
  };
}
