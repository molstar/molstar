/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type VisualTestCase } from '../framework/test-case.js';
import { ComputeCases } from './compute.js';
import { FallbackCases } from './fallback.js';
import { GeometryCases } from './geometry.js';
import { GpuPathCases } from './gpu-paths.js';
import { MergedCases } from './merged.js';
import { PostprocessingCases } from './postprocessing.js';
import { RendererCases } from './renderer.js';
import { ShaderVariantCases } from './shader-variants.js';
import { ThemeDataCases } from './theme-data.js';
import { TransparencyCases } from './transparency.js';

export const AllCases: VisualTestCase[] = [
  ...GeometryCases,
  ...MergedCases,
  ...ThemeDataCases,
  ...TransparencyCases,
  ...PostprocessingCases,
  ...RendererCases,
  ...FallbackCases,
  ...ComputeCases,
  ...GpuPathCases,
];

const ids = new Set<string>();
for (const c of AllCases) {
  if (ids.has(c.id)) throw new Error(`Duplicate visual test id '${c.id}'`);
  ids.add(c.id);
}

/** Cases run only when a filter or tags select cases (`--filter`, `--tags`): too slow for every run */
export const OnDemandCases: VisualTestCase[] = [...ShaderVariantCases];
