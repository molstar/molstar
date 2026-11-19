/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Size entry 2, headless rendering: reading a spec, then layout, draw lists and SVG. Bundled by `sizes.mjs`, never run.

import { parseJson } from '@molstar/core/repr/codecs/json';
import { checkSeqSpec } from '@molstar/seqx-core/lang';
import { Renderer } from '@molstar/seqx-core/render';
import { renderSvg } from '@molstar/seqx-core/render/render-svg';
import { resolveSpec } from '@molstar/seqx-core/resolve';

export function renderJson(json: string, width: number, height: number) {
  const decoded = parseJson(json);
  if (!decoded.node) return decoded.issues;
  return [checkSeqSpec(decoded.node), resolveSpec(decoded.node), renderSvg(decoded.node, { size: { width, height } })];
}

export { Renderer };
