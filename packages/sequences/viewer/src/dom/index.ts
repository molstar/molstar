/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `@molstar/seq-viewer/dom` (seqspec:interaction#the-dom-element): the vanilla UI, plain TypeScript over the DOM. The
// `<seq-viewer>` element (`defineSeqViewer`), `mountSeqViewer` for hosts that want no element (`createDomUI` for one that
// composes its own), and every part of seqspec:interaction#the-default-ui as a `DomPart` a host can also mount anywhere,
// with the React UI's structure, class names and skin (`@molstar/seq-viewer/skin.css`, one file for both); the default
// primitives (`Hint`, `ScrollArea`, `Menu`), the helpers a host's parts can build with (`spot`, `popover`, `hintBehavior`) and the
// browser's SVG export. Imports no React.

export type { MenuItem } from "../chrome.ts";
export * from "./core.ts";
export { downloadSvg } from "./download.ts";
export {
  createDomUI,
  defineSeqViewer,
  mountSeqViewer,
  mountSeqViewerRoot,
  type SeqHoverDetail,
  type SeqSelectionDetail,
  SeqViewerElement,
} from "./element.ts";
export { browserExportAction, canvasEmBox, compressionDeflate, withDownloadExport } from "./export.ts";
export { HoverCardPart } from "./hover-card.ts";
export * from "./parts.ts";
export { ViewSurface } from "./section.ts";
export {
  focusIndicator,
  HINT_DELAY_MS,
  HINT_GROUP_MS,
  hintBehavior,
  type Popover,
  placeHint,
  popover,
} from "./shared.ts";
