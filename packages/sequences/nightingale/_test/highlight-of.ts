/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Test support: a host's loci as one manager's `highlight` in one call, as `mount.highlight` builds it with its own
// reused extents, sink and runs (`prepareLoci`, `highlightRuns`, then the string).

import type { Closure, Projection } from "@molstar/seq-core/format";
import { Extents, type Loci } from "@molstar/seq-core/model";
import { AxisRuns, AxisSink, highlightRuns, prepareLoci } from "../src/axis.js";
import type { AxisMeta } from "../src/index.js";

export function highlightOf(p: Projection, meta: AxisMeta, loci: Loci | null, closure: Closure = "element"): string {
  const ext = Extents.create();
  if (!prepareLoci(p, ext, loci, closure)) return "";
  return AxisRuns.text(highlightRuns(p, ext, meta, new AxisSink(), AxisRuns.create()));
}
