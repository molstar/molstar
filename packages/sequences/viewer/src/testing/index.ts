/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export {
  type ConformanceCheck,
  type ConformanceFixtures,
  type ConformanceZoom,
  frameDiff,
  LayerFrame,
  type MirrorItem,
  mirrorItems,
  rendererConformance,
  svgDiff,
  TeeBackend,
  traceOf,
  unmatched,
} from './conformance.js';
export { type EdgeKind, EdgeKinds, edgeDocument } from './edges.js';
export { AffinityField, type Stars, StarsFormat, StarsRenderer } from './examples/stars.js';
export { hashRecording, Op, RecordedImages, RecordingBackend } from './recording.js';
export { type Streams, streamsOf } from './streams.js';
export { fixtureMeasurer, fixtureWidth } from './text.js';
