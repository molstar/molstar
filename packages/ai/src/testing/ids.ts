/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { IdSource } from '../harness/index.js';

/** Deterministic ids: `turn-1`, `step-1`, … per kind, prefixed so two sources never collide. */
export function simulatedIds(prefix = ''): IdSource & { readonly issued: Readonly<Record<string, number>> } {
  const issued: Record<string, number> = {};
  return {
    issued,
    next(kind) {
      issued[kind] = (issued[kind] ?? 0) + 1;
      return `${prefix}${kind}-${issued[kind]}`;
    },
  };
}
