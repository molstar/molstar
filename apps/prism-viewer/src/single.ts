/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The single-file viewer (app.md §2): one ES module, every lazy part inlined, CSS included (Prism's CSS is in its JS
// modules). A page imports it and runs it with no build step:
//   import { createPrismViewer, FullPrismViewerSpec } from './prism-viewer.js';
//   const viewer = await createPrismViewer({ target, spec: FullPrismViewerSpec() });

export { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
export { createPrismViewer } from '@molstar/prism-viewer/page';
export { FullPrismViewerSpec } from './full-spec.js';
