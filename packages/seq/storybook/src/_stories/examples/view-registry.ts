/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type SeqRegistryEntry, withAdapter } from '@molstar/seqx-core/context';

// How a story composes a view's registry over the default viewer spec's (`example-root.ts` and the headless twins).

/**
 * The default registry with a view's entries: each format runtime the view lists as its formats' adapter
 * (`withAdapter`: the GFF runtime reads `gff3` in place of the records preset), its other entries added.
 */
export function withViewRegistry<E extends { readonly kind: string }>(
  base: readonly E[],
  extra: readonly SeqRegistryEntry[],
): (E | SeqRegistryEntry)[] {
  let out: (E | SeqRegistryEntry)[] = [...base];
  for (const e of extra) out = e.kind === 'format-runtime' ? withAdapter(out, e.runtime) : [...out, e];
  return out;
}
