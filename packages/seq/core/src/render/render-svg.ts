/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { BuiltInLayerTypes, BuiltInLayouts } from './built-ins.js';
import { exportSvg, type RenderSvgOptions, type RenderSvgResult } from './export.js';
import { NavigatorSvgExport } from './navigator-export.js';

// The headless default (_spec/README.md "Render", Export): `exportSvg` with everything the core ships, the built-in
// layer types, layouts and navigators, unless the options name their own. A composition that lists fewer calls
// `exportSvg` (a viewer, `createSeqContext`'s base) and bundles fewer.

/** A spec drawn to an SVG string once the data its view needs has settled (or the timeout passed): the built-ins. */
export function renderSvg(spec: Node, options: RenderSvgOptions): Promise<RenderSvgResult> {
  return exportSvg(spec, {
    types: BuiltInLayerTypes,
    layouts: BuiltInLayouts,
    navigatorExport: NavigatorSvgExport,
    ...options,
  });
}
