/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `@molstar/prism-viewer/api` (api.md §6): the API's types and helpers, and the module an extension augments with its
// areas: `declare module '@molstar/prism-viewer/api' { interface ExtensionApis { 'geo-export': GeoExportApi } }`.

// What a page, a script or an extension uses (a curated list, api.md §6): the result and operation types, an area's
// base class, the values check, node selectors, the targets' types. The runner and the API's construction stay inside.
export {
  ApiArea,
  type ApiAreaContext,
  type ApiContext,
  type ApiErrorCode,
  type ApiFailure,
  type ApiModuleKey,
  ApiResult,
  type BuiltInApis,
  type Changed,
  type CoreApi,
  type GroupOptions,
  type GroupResult,
  memberName,
  type OperationDefinition,
  type OperationDescription,
  type OperationFn,
  type RunOptions,
  type RunOrigin,
  type RunOriginKind,
  type ViewerApi,
} from './api/area.js';
export type { ViewerCommandContext, ViewerManagers } from './api/managers.js';
export {
  type NodeAttrs,
  type NodeKind,
  type NodeQuery,
  NodeKinds,
  node,
  parseNodeSelector,
  selectNodes,
} from './api/nodes.js';
export { checkValues, messageOf, type RunRecord } from './api/run.js';
export type { StructureId } from './api/scope.js';
export type { ResolvedPart, ResolvedTarget, ResolveOptions, TargetResolver } from './api/targets.js';

/** Extensions' areas, by id: an extension augments this interface here, so `viewer.api['geo-export']` is typed. */
export interface ExtensionApis {}
