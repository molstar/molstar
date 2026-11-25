/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { withViewKind } from '@molstar/seqx-core/resolve';

// The values layer's reference lines' door (_spec/README.md "Render", Layer types, `values`: lazy): the lines' reading,
// layout, drawing, labels and words (`values-reference.ts`) load the first time a spec holds a values layer with a
// `reference`, a model waiting for them as for the bridge chunk (`withViewKind`), so such a layer is never drawn
// without its lines. Why lazy: most specs have none; a view without one pays nothing for them.

/** The chunk's module once landed. */
export const ReferenceHooks: { chunk?: typeof import('./values-reference.js') } = {};

/** Loads the reference lines' chunk. */
export function loadReference(): Promise<void> {
  return import('./values-reference.js').then((m) => {
    ReferenceHooks.chunk = m;
  });
}

/** Whether a section's tracks (or its groups' tracks) hold a layer with `reference` lines. */
export function holdsReference(section: Node): boolean {
  const layers = (n: Node): boolean =>
    (n.children ?? []).some((c) =>
      c.kind === 'layer' ? Array.isArray(c.params?.reference) && c.params.reference.length > 0 : layers(c),
    );
  return layers(section);
}

let registered = false;

/** Registers the chunk with resolve (the values layer's module calls it): loaded when a spec holds a layer with lines. */
export function registerReference(): void {
  if (registered) return;
  registered = true;
  withViewKind('section', loadReference, holdsReference);
}
