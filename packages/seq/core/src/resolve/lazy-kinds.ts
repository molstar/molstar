/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Parts loaded when a spec holds a view-level kind (_spec/README.md "Resolve", Lazy kinds): render's bridge chunk is
// registered for `bridge` (`withViewKind`), and a model whose spec has a bridge waits for it (`StaticModel.loading`)
// as for the genome chunk, so no view with a bridge is drawn without it. Why: a reader that never draws a bridge pays
// nothing for one.

import type { Node } from '@molstar/core/repr/tree';

interface Part {
  readonly load: () => Promise<unknown>;
  landed: boolean;
  loading?: Promise<void>;
}

const Parts = new Map<string, Part[]>();

/** A part to load when a spec's view holds a node of `kind` (render's bridge chunk: `bridge`). */
export function withViewKind(kind: string, load: () => Promise<unknown>): void {
  let l = Parts.get(kind);
  if (!l) Parts.set(kind, (l = []));
  l.push({ load, landed: false });
}

/**
 * The loading of the parts a spec's views need and that have not landed, or undefined when none: it never rejects (a
 * failed load is tried again by the next model that needs it).
 */
export function viewKindsLoading(spec: Node): Promise<void> | undefined {
  if (!Parts.size) return undefined;
  let out: Promise<void>[] | undefined;
  for (const views of spec.children ?? []) {
    if (views.kind !== 'views') continue;
    for (const view of views.children ?? [])
      for (const c of view.children ?? []) {
        const parts = Parts.get(c.kind);
        if (!parts) continue;
        for (const p of parts) {
          if (p.landed) continue;
          p.loading ??= p.load().then(
            () => {
              p.landed = true;
              p.loading = undefined;
            },
            () => {
              p.loading = undefined;
            },
          );
          if (!out?.includes(p.loading)) (out ??= []).push(p.loading);
        }
      }
  }
  return out && Promise.all(out).then(() => undefined);
}
