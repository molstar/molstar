/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type VisualTestCase } from '../framework/test-case.js';
import { FallbackCases } from './fallback.js';
import { GeometryCases } from './geometry.js';
import { PostprocessingCases } from './postprocessing.js';
import { RendererCases } from './renderer.js';
import { TransparencyCases } from './transparency.js';

export const AllCases: VisualTestCase[] = [
  ...GeometryCases,
  ...TransparencyCases,
  ...PostprocessingCases,
  ...RendererCases,
  ...FallbackCases,
];

const ids = new Set<string>();
for (const c of AllCases) {
  if (ids.has(c.id)) throw new Error(`Duplicate visual test id '${c.id}'`);
  ids.add(c.id);
}
