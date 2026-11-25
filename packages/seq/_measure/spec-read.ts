/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Size entry 1, reading a spec: decode JSON, validate, resolve. Bundled by `sizes.mjs`, never run.

import { parseJson } from '@molstar/core/repr/codecs/json';
import { checkSeqSpec } from '@molstar/seqx-core/lang';
import { resolveSpec } from '@molstar/seqx-core/resolve';

export function readSpec(json: string) {
  const decoded = parseJson(json);
  if (!decoded.node) return { issues: decoded.issues };
  const issues = checkSeqSpec(decoded.node);
  return { issues, model: resolveSpec(decoded.node) };
}
