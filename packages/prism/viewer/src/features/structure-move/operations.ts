/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '@molstar/model/model/structure';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import { ApiArea, type ViewerApi } from '../../api/area.js';
import { labelOf } from '../../api/scope.js';
import { StructureMoveGizmos } from '../../sources/structure-move/gizmos.js';
import { watchedState } from '../../sources/watched.js';

// The `move` area (structure-move.md): a structure moved and turned by hand with handles in 3D. It adds no operations:
// a drag's release is `structure.transform` (the matrix from where it was loaded), Reset `structure.transform-reset`;
// it holds the handles and what can be moved.

export interface MovableStructure {
  readonly ref: string;
  readonly label: string;
}

export class MoveApi extends ApiArea {
  readonly id = 'move';

  /** The handles on a structure (the section's), their drags committed through the person's API */
  readonly gizmos = new StructureMoveGizmos(this.ctx.plugin);

  /** The structures loaded, in the hierarchy's order */
  readonly structures: ReadonlyState<readonly MovableStructure[]> = watchedState(
    () =>
      this.ctx.plugin.managers.structure.hierarchy.current.structures.map((s) => ({
        ref: s.cell.transform.ref,
        label: labelOf(this.ctx, s),
      })),
    () => [this.ctx.plugin.state.data.events.changed],
    { compare: (a, b) => a.length === b.length && a.every((s, i) => s.ref === b[i].ref && s.label === b[i].label) },
  );

  /** The structure the selection is in (the first such), or none */
  selectedStructure(): string | undefined {
    const selection = this.ctx.plugin.managers.structure.selection;
    for (const s of this.ctx.plugin.managers.structure.hierarchy.current.structures) {
      const data = s.cell.obj?.data;
      const loci = data && selection.getLoci(data);
      if (loci && StructureElement.Loci.is(loci) && !StructureElement.Loci.isEmpty(loci)) return s.cell.transform.ref;
    }
    return undefined;
  }

  /** Shows the handles, a drag committing through `api` (the person's: one `structure.transform`, one history step) */
  showGizmos(api: ViewerApi): boolean {
    return this.gizmos.show(async (structure, matrix) => {
      const r = await api.run('structure.transform', { structure, matrix: [...matrix] });
      return r.ok;
    });
  }
}
