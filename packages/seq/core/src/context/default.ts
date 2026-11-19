/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BuiltInLayerTypes, BuiltInLayouts } from '@molstar/seqx-core/render/built-ins';
import { NavigatorSvgExport } from '@molstar/seqx-core/render/navigator-export';
import { composeSeqContext, type SeqContext, type SeqContextSpec, type SeqRegistryEntry } from './context.js';

// The default headless composition (_spec/README.md "Context"): the core's layer types and layouts listed before the
// spec's own registry, and its navigators in the SVG export. A smaller composition calls `composeSeqContext` with its
// own list and bundles only that.

/** The core's layer types and section layouts as registry entries: what the default compositions list. */
// (pure: a composition that never reads it bundles none of the built-ins, whatever its bundler knows of side effects)
export const BuiltInRenderRegistry: readonly SeqRegistryEntry[] = /* @__PURE__ */ (() => [
  ...Object.values(BuiltInLayerTypes).map((type) => ({ kind: 'layer-type' as const, type })),
  ...Object.values(BuiltInLayouts).map((layout) => ({ kind: 'layout' as const, layout })),
])();

/**
 * The headless context of the default composition: `composeSeqContext` with the built-in layer types and layouts
 * registered before the spec's own entries (the same object again counts, another under a built-in's name is an error)
 * and the built-in navigators drawn in its exports.
 */
export function createSeqContext(spec: SeqContextSpec = {}): SeqContext {
  const context = composeSeqContext({ ...spec, registry: [...BuiltInRenderRegistry, ...(spec.registry ?? [])] });
  return {
    ...context,
    renderSvg: (s, options) => context.renderSvg(s, { navigatorExport: NavigatorSvgExport, ...options }),
  };
}
