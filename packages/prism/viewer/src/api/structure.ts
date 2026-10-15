/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '@molstar/prism-ui';
import type { StructureInfo } from '../sources/viewer-sources.js';
import { ApiArea, ApiResult } from './area.js';

// `structure` (operations.md §2): `get` (the scene outline) now; `preset`, `model`, `assembly` with the scene areas
// (plan phase 2.2).

export class StructureApi extends ApiArea {
  readonly id = 'structure';

  /** The loaded structures. */
  get list(): ReadonlyState<readonly StructureInfo[]> {
    return this.ctx.sources.structures;
  }

  readonly get = this.operation({
    id: 'structure.get',
    title: 'Describe the scene',
    params: {},
    read: true,
    async run(_, ctx) {
      const list = ctx.sources.structures.value;
      if (!list.length) return ApiResult.ok('Nothing is loaded.', {}, { value: list });
      const lines = list.map(
        (s) =>
          `${s.label} (${s.ref}): ${s.atoms.toLocaleString('en')} atoms, ${s.residues.toLocaleString('en')} residues, ${s.chains} chain${s.chains === 1 ? '' : 's'}${s.components.length ? `; ${s.components.join(', ')}` : ''}`,
      );
      return ApiResult.ok(lines.join('\n'), {}, { value: list });
    },
  });
}
