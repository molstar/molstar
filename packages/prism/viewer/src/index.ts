/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `@molstar/prism-viewer`: what a page embedding the viewer uses (viewer-spec.md §6): create a viewer from a spec, the
// default spec, and the types of the spec, the viewer and its API. The rest is on narrower subpaths (package.json's
// exports): `/headless` (Node), `/api` (an extension's areas and results), `/spec`, `/features/<name>`, the views an
// extension draws with (`/views/structure/scene-row`).

export { type CreatePrismViewerOptions, createPrismViewer } from './page.js';
export { DefaultPrismViewerSpec } from './default-spec.js';
export type { PrismViewer } from './viewer.js';
export { type CameraSnapshot, type SetCameraOptions, syncCameras, type ViewerCamera } from './sources/camera.js';
export type {
  LazyViewerExtension,
  PrismViewerSpec,
  ViewerExtension,
  ViewerExtensionUI,
  ViewerSkill,
  ViewerUIContext,
} from './spec.js';
export {
  ApiResult,
  type ApiFailure,
  type NodeQuery,
  node,
  type OperationDescription,
  type RunOrigin,
  type ViewerApi,
} from './api.js';
export type { SceneNode } from './model/scene.js';
// The sequence's types; its values (Seq* and the views over it) are its own chunk: `@molstar/prism-viewer/sequence/index`.
export type * from './sequence/index.js';
