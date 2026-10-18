/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import './favicon.ico';
import './index.html';
import { setDebugMode, setTimingMode } from '@molstar/core/util/debug';
import { LoadUrlParameters, parseSearch, runUrlLoads } from '@molstar/prism-viewer/features/loaders';
import { createPrismViewer } from '@molstar/prism-viewer/page';
import { createAppSpec } from './spec.js';
import { modesFromUrl, unknownUrlParameters } from './url.js';

// The Prism viewer app (app.md §1): the page, the spec the app ships (the viewer's defaults, Mol*'s extensions, the
// URL's options) and the viewer, exposed on `window.molstar` for the console and for scripts on the page with Mol*'s
// library (`lib`, app.md §1, loaded beside the viewer in a chunk of its own). The URL's load parameters run once the
// viewer is created (app.md §4), as the loaders feature's operations.

async function main() {
  const target = document.getElementById('app');
  if (!target) throw new Error('The page has no #app element.');
  // Page-wide, before the plugin starts (app.md §4): v5 sets them only when asked, and so does the app.
  const modes = modesFromUrl(location.search);
  if (modes.debug) setDebugMode(true);
  if (modes.timing) setTimingMode(true);
  // The library is not part of the start bundle: its chunk loads while the viewer starts; a failed load leaves it out.
  const lib = import('./lib.js').then(
    (m) => m.lib,
    (e: unknown) => console.error('molstar.lib could not be loaded', e),
  );
  const viewer = await createPrismViewer({ target, spec: createAppSpec(location.search) });
  (window as unknown as { molstar: unknown }).molstar = { viewer, api: viewer.api, lib: await lib };

  const unknown = unknownUrlParameters(location.search, LoadUrlParameters);
  if (unknown.length) console.info(`Ignored URL parameters (not known to the viewer): ${unknown.join(', ')}`);
  // Without the loaders feature in the spec there is nothing to run them with.
  const params = parseSearch(location.search);
  const urlLoads = LoadUrlParameters.some((p) => params.has(p));
  if (urlLoads && viewer.api.describe('load.entry')) {
    runUrlLoads(viewer.api, location.search).catch((e: unknown) => console.error(e));
  } else if (!urlLoads && isReload() && viewer.restore?.available) {
    // The light restore (journal.md §6): only on a reload, and only when the URL asks for nothing to load.
    viewer.restore.run().catch((e: unknown) => console.error(e));
  }
}

function isReload(): boolean {
  const nav = performance.getEntriesByType('navigation')[0] as PerformanceNavigationTiming | undefined;
  return nav?.type === 'reload';
}

main().catch((e) => {
  console.error(e);
  const p = document.createElement('p');
  p.style.cssText = 'font: 14px system-ui; padding: 24px';
  p.textContent = `The viewer could not start: ${e instanceof Error ? e.message : String(e)}`;
  document.body.appendChild(p);
});
