/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import './favicon.ico';
import './index.html';
import { LoadUrlParameters, runUrlLoads } from '@molstar/prism-viewer/features/loaders';
import { createPrismViewer } from '@molstar/prism-viewer/page';
import { createAppSpec } from './spec.js';
import { unknownUrlParameters } from './url.js';

// The Prism viewer app (app.md §1): the page, the spec the app ships (the viewer's defaults, Mol*'s extensions, the
// URL's options) and the viewer, exposed on `window.molstar` for the console and for scripts on the page. The URL's
// load parameters run once the viewer is created (app.md §4), as the loaders feature's operations.

async function main() {
  const target = document.getElementById('app');
  if (!target) throw new Error('The page has no #app element.');
  const viewer = await createPrismViewer({ target, spec: createAppSpec(location.search) });
  (window as unknown as { molstar: unknown }).molstar = { viewer, api: viewer.api };

  const unknown = unknownUrlParameters(location.search, LoadUrlParameters);
  if (unknown.length) console.info(`Ignored URL parameters (not known to the viewer): ${unknown.join(', ')}`);
  // Without the loaders feature in the spec there is nothing to run them with.
  if (viewer.api.describe('load.entry')) {
    runUrlLoads(viewer.api, location.search).catch((e: unknown) => console.error(e));
  }
}

main().catch((e) => {
  console.error(e);
  const p = document.createElement('p');
  p.style.cssText = 'font: 14px system-ui; padding: 24px';
  p.textContent = `The viewer could not start: ${e instanceof Error ? e.message : String(e)}`;
  document.body.appendChild(p);
});
