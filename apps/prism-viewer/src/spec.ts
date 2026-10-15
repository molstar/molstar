/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PLUGIN_VERSION } from '@molstar/plugin/version';
import { webStore } from '@molstar/prism-viewer/persist/stores';
import type { PrismViewerSpec } from '@molstar/prism-viewer/spec';
import { FullPrismViewerSpec } from './full-spec.js';
import { optionsFromUrl } from './url.js';

/** The spec the app ships (viewer-spec.md §5): the viewer's defaults, Mol*'s extensions and the URL's options. */
export function createAppSpec(search: string): PrismViewerSpec {
  const base = FullPrismViewerSpec();
  const params = new URLSearchParams(search);
  // `?assistant`: the Assistant tab open at start, which loads it (otherwise it loads on first use: the tab, the Ask row).
  const assistant = params.has('assistant');
  // `?restore`: this tab's scene comes back after a reload (journal.md §6). Off by default: a session lives in memory.
  const restore = params.has('restore');
  return {
    ...base,
    options: optionsFromUrl(search),
    ...(assistant ? { ui: { ...base.ui, layout: { ...base.ui?.layout, start: 'assistant.tab.assistant' } } } : {}),
    // Light state in the browser (journal.md §6): preferences in localStorage; with `?restore`, this tab's light
    // restore in sessionStorage; both stamped with the app's version.
    persist: {
      key: 'prism-viewer',
      preferences: webStore('local', { prefix: 'prism-viewer', version: PLUGIN_VERSION, budget: 64 * 1024 }),
      ...(restore
        ? { restore: webStore('session', { prefix: 'prism-viewer', version: PLUGIN_VERSION, budget: 256 * 1024 }) }
        : {}),
    },
    // The assistant, its own chunk, loaded on first use (viewer-spec.md §2.4).
    agent: { ...(base.agent || {}), assistant: () => import('@molstar/prism-assistant/viewer') },
  };
}
