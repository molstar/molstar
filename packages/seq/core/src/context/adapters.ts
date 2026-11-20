/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { FormatRuntime } from '@molstar/seqx-core/data';
import type { SeqRegistryEntry } from './context.js';

type FormatRuntimeEntry = Extract<SeqRegistryEntry, { readonly kind: 'format-runtime' }>;

// Data adapters (_spec/README.md "Context", Registry entries): a format runtime is the default adapter of the formats
// it reads; a host replaces it for its formats. One adapter a format, so which one reads a source never depends on
// the order the registry lists them in.

/**
 * The registry with `adapter` reading its formats: each format runtime listed for one of them gives it up (a runtime
 * left with none goes, one with others keeps them under its name), then the adapter is added. The other entries are
 * kept as they are, in order. Generic over the entries a host's registry holds (a viewer's lists its own kinds).
 */
export function withAdapter<E extends { readonly kind: string }>(
  registry: readonly E[],
  adapter: FormatRuntime,
): (E | FormatRuntimeEntry)[] {
  const taken = new Set(adapter.formats);
  // one narrowing a runtime, however often it is listed: the same object, so `runtimesOf` counts it once
  const narrowings = new Map<FormatRuntime, FormatRuntime>();
  const out: (E | FormatRuntimeEntry)[] = [];
  for (const e of registry) {
    // an entry of a host's own kind, or one that names the kind without a runtime, is not a runtime: kept as it is
    if (e.kind !== 'format-runtime' || !Array.isArray((e as Partial<FormatRuntimeEntry>).runtime?.formats)) {
      out.push(e);
      continue;
    }
    const r = (e as unknown as FormatRuntimeEntry).runtime;
    if (r === adapter) continue;
    const rest = r.formats.filter((f) => !taken.has(f));
    if (rest.length === r.formats.length) out.push(e);
    else if (rest.length) {
      let n = narrowings.get(r);
      if (!n) narrowings.set(r, (n = narrowed(r, rest)));
      out.push({ kind: 'format-runtime', runtime: n });
    }
  }
  out.push({ kind: 'format-runtime', runtime: adapter });
  return out;
}

/** A runtime reading only `formats` of its own, under its name. */
function narrowed(runtime: FormatRuntime, formats: readonly string[]): FormatRuntime {
  return { name: runtime.name, formats, open: (source, ctx) => runtime.open(source, ctx) };
}
