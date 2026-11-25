/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The entry by its package name (_spec/README.md "Structures"): a host imports `@molstar/seqx-core/structure`, so the
// package's exports must name it (the `./*` wildcard maps to `src/*.ts`, which holds no `structure.ts`).

import * as Structure from '@molstar/seqx-core/structure';
import { describe, expect, it } from 'vitest';

describe('@molstar/seqx-core/structure', () => {
  it('resolves by its package name and names its public surface', () => {
    expect(Object.keys(Structure).sort()).toEqual([
      'DefaultSyncLimits',
      'StructureSyncIssues',
      'colorRuns',
      'fromUnits',
      'structureFacts',
      'toUnits',
      'unitRows',
    ]);
  });
});
