/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */
// The React bindings (react/_spec): hooks over core's reactive sources, a React view for the UI's view hosts, and a
// host for the UI's elements in a React tree. React (or preact/compat) is an optional peer dependency.

export { UiElementHost, type UiElementHostProps } from './element-host.js';
export { type ReactRoot, reactView } from './react-view.js';
export { useReactive } from './use-reactive.js';
export { useReactiveModel } from './use-reactive-model.js';
export { useReactiveProp } from './use-reactive-prop.js';
export { type ReactiveOptions, type ReactiveSource, useReactiveSelection } from './use-reactive-selection.js';
