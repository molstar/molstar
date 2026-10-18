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
// Importing it also defines `<prism-viewer>` (app.md §3), with `FullPrismViewerSpec()` while its `spec` is not set:
//   <prism-viewer pdb="1cbs" theme="dark" style="height: 480px"></prism-viewer>
// v5's `molstar.Viewer.create` script global is not carried over (parity.md §9, breaking-v6-changes.md).

import { definePrismViewerElement } from '@molstar/prism-viewer/element';
import { FullPrismViewerSpec } from './full-spec.js';

export { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
export {
  createPrismViewerElementClass,
  definePrismViewerElement,
  type PrismViewerElement,
} from '@molstar/prism-viewer/element';
export { createPrismViewer } from '@molstar/prism-viewer/page';
export { FullPrismViewerSpec };

// Once per page: a page that imports the file twice (or defined the tag itself) keeps the first definition.
definePrismViewerElement({ defaultSpec: FullPrismViewerSpec });
