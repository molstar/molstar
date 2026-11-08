/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement } from '@molstar/model/model/structure';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import { ApiArea, type ViewerApi } from '../../api/area.js';
import { labelOf, namesOf } from '../../api/scope.js';
import { type MovableKind, MoveGizmos, movableVolumes } from '../../sources/move/gizmos.js';
import { shapeEntries } from '../../sources/shapes.js';
import { watchedState } from '../../sources/watched.js';

// The `move` area (move.md): a model (a structure, a particle list, a shape, a volume) moved and turned by hand with
// handles in 3D. It adds no operations: a drag's release is its kind's transform (`structure.transform`,
// `particles.transform`, `shape.transform`, `volume.transform`: the matrix from where it was loaded), Reset their
// `transform-reset`; it holds the handles and what can be moved.

/** Each kind's transform operation and the name of its argument naming the model */
const Operations: Record<MovableKind, { readonly id: string; readonly arg: string }> = {
  structure: { id: 'structure', arg: 'structure' },
  particles: { id: 'particles', arg: 'particles' },
  shape: { id: 'shape', arg: 'node' },
  volume: { id: 'volume', arg: 'volume' },
};

export interface Movable {
  readonly ref: string;
  readonly label: string;
  readonly kind: MovableKind;
}

export class MoveApi extends ApiArea {
  readonly id = 'move';

  /** The handles on a model (the section's), their drags committed through the person's API */
  readonly gizmos = new MoveGizmos(this.ctx.plugin);

  /** The structures loaded, in the hierarchy's order, then the particle lists, the shapes, the volumes on their own */
  readonly movable: ReadonlyState<readonly Movable[]> = watchedState(
    () => {
      const { structure, particles } = this.ctx.plugin.managers;
      const lists = particles.hierarchy.current.lists;
      const names = namesOf(lists);
      const volumes = movableVolumes(this.ctx.plugin);
      return [
        ...structure.hierarchy.current.structures.map((s) => ({
          ref: s.cell.transform.ref,
          label: labelOf(this.ctx, s),
          kind: 'structure' as const,
        })),
        ...lists.map((l) => ({
          ref: l.cell.transform.ref,
          label: names.get(l.cell.transform.ref) ?? l.cell.transform.ref,
          kind: 'particles' as const,
        })),
        ...shapeEntries(this.ctx.plugin.state.data).map((e) => ({
          ref: e.cell.transform.ref,
          label: e.label,
          kind: 'shape' as const,
        })),
        ...volumes.list.map((v) => ({
          ref: v.cell.transform.ref,
          label: volumes.names.get(v.cell.transform.ref) ?? v.cell.transform.ref,
          kind: 'volume' as const,
        })),
      ];
    },
    () => [this.ctx.plugin.state.data.events.changed],
    {
      compare: (a, b) =>
        a.length === b.length &&
        a.every((s, i) => s.ref === b[i].ref && s.label === b[i].label && s.kind === b[i].kind),
    },
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

  /** Shows the handles, a drag committing through `api` (the person's: one transform operation, one history step) */
  showGizmos(api: ViewerApi): boolean {
    return this.gizmos.show(async (ref, kind, matrix) => {
      const { id, arg } = Operations[kind];
      return (await api.run(`${id}.transform`, { [arg]: ref, matrix: [...matrix] })).ok;
    });
  }

  /** Puts the model back where it was loaded */
  reset(api: ViewerApi, ref: string) {
    const { id, arg } = Operations[this.gizmos.kindOf(ref)];
    return api.run(`${id}.transform-reset`, { [arg]: ref });
  }
}
