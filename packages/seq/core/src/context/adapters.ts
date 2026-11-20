/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { FormatRuntime } from '@molstar/seqx-core/data';
import type { SeqRegistryEntry } from './context.js';

// Data adapters (_spec/README.md "Context", Registry entries): a format runtime is the default adapter of the formats
// it reads; a host replaces it for its formats. One adapter a format, so which one reads a source never depends on
// the order the registry lists them in.

/**
 * The registry with `adapter` reading its formats: each format runtime listed for one of them gives it up (a runtime
 * left with none goes, one with others keeps them under its name), then the adapter is added. The other entries are
 * kept as they are, in order.
 */
export function withAdapter(registry: readonly SeqRegistryEntry[], adapter: FormatRuntime): SeqRegistryEntry[] {
  const taken = new Set(adapter.formats);
  const out: SeqRegistryEntry[] = [];
  for (const e of registry) {
    if (e.kind !== 'format-runtime') out.push(e);
    if (e.kind !== 'format-runtime' || e.runtime === adapter) continue;
    const r = e.runtime;
    const rest = r.formats.filter((f) => !taken.has(f));
    if (rest.length === r.formats.length) out.push(e);
    else if (rest.length) out.push({ kind: 'format-runtime', runtime: narrowed(r, rest) });
  }
  out.push({ kind: 'format-runtime', runtime: adapter });
  return out;
}

/** A runtime reading only `formats` of its own, under its name. */
function narrowed(runtime: FormatRuntime, formats: readonly string[]): FormatRuntime {
  return { name: runtime.name, formats, open: (source, ctx) => runtime.open(source, ctx) };
}
