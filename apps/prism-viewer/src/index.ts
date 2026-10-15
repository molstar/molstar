/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import './favicon.ico';
import './index.html';
import { createPrismViewer } from '@molstar/prism-viewer/page';
import { createAppSpec } from './spec.js';

// The Prism viewer app (app.md §1): the page, the spec the app ships (the viewer's defaults, Mol*'s extensions, the
// URL's options) and the viewer, exposed on `window.molstar` for the console and for scripts on the page.

async function main() {
  const target = document.getElementById('app');
  if (!target) throw new Error('The page has no #app element.');
  const viewer = await createPrismViewer({ target, spec: createAppSpec(location.search) });
  (window as unknown as { molstar: unknown }).molstar = { viewer, api: viewer.api };
}

main().catch((e) => {
  console.error(e);
  const p = document.createElement('p');
  p.style.cssText = 'font: 14px system-ui; padding: 24px';
  p.textContent = `The viewer could not start: ${e instanceof Error ? e.message : String(e)}`;
  document.body.appendChild(p);
});
