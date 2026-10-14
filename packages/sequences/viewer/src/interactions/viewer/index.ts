/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `@molstar/seq-viewer/interactions/viewer` (seqspec:viewer#extension-points): the extension's viewer half. One registry entry
// holds both halves: the format's (`InteractionsRegistry`: the provider, the layer type, the theme entry) and the
// viewer's (the `interactions:fingerprint` renderer with its draws, the `interactions:kind` theme provider). Headless:
// a registry module like the built-ins, so a canvas host and a headless one register the same entry. The viewer core
// imports nothing from here.

import type { LayerRenderer, SeqViewerRegistryEntry } from "@molstar/seq-viewer/draw";
import { InteractionsRegistry } from "../registry.js";
import { FingerprintRenderer } from "./renderer.js";
import { KindTheme } from "./themes.js";

export {
  FINGERPRINT,
  FingerprintRenderer,
  type FingerprintState,
  prepareFingerprint,
} from "./renderer.js";
export { KindTheme } from "./themes.js";

/** The extension's one registry entry: what a host adds to `spec.registry`, headless or drawing. */
export const InteractionsViewerRegistry: SeqViewerRegistryEntry = {
  svs: InteractionsRegistry,
  renderers: [FingerprintRenderer as LayerRenderer],
  themes: [KindTheme],
};
