/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer of Seq* (next): DOM-free models on @molstar/core/reactive, the Canvas2D and SVG backends, and the default
// UI on @molstar/ui. What it holds and how it is layered: _spec/README.md.

import { LanguageVersion } from '@molstar/seqx-core';

/** The SeqViewSpec language version this viewer draws. */
export const ViewerLanguageVersion = LanguageVersion;

export { type KeyState, ViewerInput, type WheelState } from './input.js';
export { createSeqViewer, SeqViewer, type ViewerCommand } from './viewer.js';
export { type NavigatorDef, NavigatorModel, type NavigatorPointer, type NavigatorValue } from './navigator.js';
export {
  type BehaviorContext,
  type BehaviorProvider,
  type SeqViewerBehavior,
  type SeqViewerBindings,
  type SeqViewerSpec,
  type SeqViewerUI,
} from './viewer-spec.js';
