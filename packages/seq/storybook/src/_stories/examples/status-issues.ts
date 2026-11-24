/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue } from '@molstar/core/repr/tree';

/** The issues a story's status line says: the export's, and the viewer's model's it lacks (by code and address). */
export function statusIssues(exported: readonly Issue[], model: readonly Issue[]): readonly Issue[] {
  const extra = model.filter((m) => !exported.some((e) => e.code === m.code && e.address === m.address));
  return extra.length ? [...exported, ...extra] : exported;
}
