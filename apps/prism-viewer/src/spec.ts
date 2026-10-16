/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PLUGIN_VERSION } from '@molstar/plugin/version';
import { assistant } from '@molstar/prism-assistant/spec';
import { webStore } from '@molstar/prism-viewer/persist/stores';
import type { PrismViewerSpec } from '@molstar/prism-viewer/spec';
import { FullPrismViewerSpec } from './full-spec.js';
import { LoadUrlParameters, parseSearch } from '@molstar/prism-viewer/features/loaders';
import { optionsFromUrl } from './url.js';

/** The spec the app ships (viewer-spec.md §5): the viewer's defaults, Mol*'s extensions and the URL's options. */
export function createAppSpec(search: string): PrismViewerSpec {
  const base = FullPrismViewerSpec();
  // Names lowercased, as the loads read them (v5's matched case-insensitively).
  const params = parseSearch(search);
  // `?assistant`: the Assistant tab open at start, which loads it (otherwise it loads on first use: the tab, the Ask row).
  const openAssistant = params.has('assistant');
  // `?restore`: this tab's scene comes back after a reload (journal.md §6). Off by default: a session lives in memory.
  const restore = params.has('restore');
  // Opened with nothing to load: Home, where loading starts (a page that drives the viewer itself keeps it closed).
  const options = optionsFromUrl(search);
  const loads = LoadUrlParameters.some((p) => params.has(p));
  return {
    ...base,
    // The app is on the web: remote states (the state server's Save and Open) are on.
    options: { ...options, ui: { homeWhenEmpty: !loads, showRemoteStates: true, ...options.ui } },
    ...(openAssistant ? { ui: { ...base.ui, layout: { ...base.ui?.layout, start: 'assistant.tab.assistant' } } } : {}),
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
    agent: { ...(base.agent || {}), assistant: assistant() },
  };
}
