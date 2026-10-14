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
} from "./conformance.ts";
export { type EdgeKind, EdgeKinds, edgeDocument } from "./edges.ts";
export { AffinityField, type Stars, StarsFormat, StarsRenderer } from "./examples/stars.ts";
export { hashRecording, Op, RecordedImages, RecordingBackend } from "./recording.ts";
export { type Streams, streamsOf } from "./streams.ts";
export { fixtureMeasurer, fixtureWidth } from "./text.ts";
