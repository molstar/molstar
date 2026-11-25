/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `@molstar/ui/dom`: a JSX DOM runtime on core's reactive base (`_spec/README.md`). JSX: `"jsxImportSource":
// "@molstar/ui/dom"`. Named re-exports only (no `export *`: the cycle rule).

export {
  Dynamic,
  For,
  Match,
  nextSlice,
  Show,
  Switch,
  activity,
  sliceHasRoom,
  timeSlicing,
} from './control-flow.js';
export {
  defineElements,
  type ElementRoot,
  ViewElement,
  type ViewElementClass,
  type ViewElementProps,
} from './element.js';
export { activeElement, eventTarget, isEditable, isInside } from './events.js';
export {
  deepActiveElement,
  focusables,
  focusedByKeyboard,
  focusFirstItem,
  menuItems,
  menuKeys,
  rememberFocus,
  rovingFocus,
  rovingTabs,
  trapTab,
  uniqueId,
} from './focus.js';
export {
  cancelFrameCallback,
  FrameCounts,
  FrameState,
  flushFrame,
  pendingFrameCallbacks,
  requestFrame,
  requestFrameCallback,
} from './frame.js';
export { DomHost } from './host.js';
export type {
  Child,
  CommonProps,
  Component,
  CustomElements,
  ElementProps,
  EventProps,
  SVGProps,
} from './jsx-types.js';
export { Motion, type RevealKind } from './motion.js';
export { DomCounts, DomProbe } from './probe.js';
export {
  createElement,
  delegateEvents,
  Fragment,
  insert,
  insertRegion,
  Passthrough,
  regionWriter,
  removeBetween,
  render,
  toNode,
} from './render.js';
export {
  createRoot,
  getScope,
  type MaybeReactive,
  onActive,
  onCleanup,
  requireScope,
  runInScope,
  toStream,
  ViewScope,
} from './scope.js';
export {
  type DragHandlers,
  drag,
  frame,
  frames,
  latest,
  observeIntersection,
  observeResize,
  on,
  type Restartable,
  restartable,
  timeout,
  workQueue,
} from './scoped.js';
export { TreeNavigation, type TreeNavigationOptions } from './tree-keys.js';
export { domView, elementView, type ViewRenderer } from './view.js';
export {
  bind,
  effect,
  flush,
  follow,
  frameState,
  type WatchHandle,
  type WatchSource,
  watch,
} from './watch.js';
export {
  createLatest,
  createWorkQueue,
  type Latest,
  type WantOptions,
  type WorkLane,
  WorkQueue,
} from './work.js';
export { setAttr, setClass, setShown, setStyle, setText } from './writes.js';
