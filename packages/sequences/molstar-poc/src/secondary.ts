/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `molstar:secondary_structure` (seqspec:libraries#the-mol-adapter): a structure's helices, strands and turns as a
// `ranges` annotation on its complex's polymer units, from the secondary structure its model carries (mmCIF's
// `struct_conf` and `struct_sheet_range`: what Mol*'s own "model" assignment reads). Coil is no range. No computation
// beyond reading, and no network.

import { type ProviderEntry, RequiredField, ref, SimpleParamsSchema } from '@molstar/seq-core/format';
import { ModelSecondaryStructure } from '@molstar/model/formats/structure/property/secondary-structure';
import { SecondaryStructureType } from '@molstar/model/model/structure/model/types';
import { chainFor, complexUnits } from './seq-units.js';
import { StructureSource, structureOf } from './source.js';

export type SecondaryKind = 'helix' | 'strand' | 'turn';

/** A residue's class from Mol*'s flags: a helix (any), a strand (beta), a turn or bend; null: coil. */
export function secondaryKind(type: SecondaryStructureType): SecondaryKind | null {
  if (
    SecondaryStructureType.is(type, SecondaryStructureType.Flag.Helix) ||
    SecondaryStructureType.is(type, SecondaryStructureType.Flag.DoubleHelix)
  )
    return 'helix';
  if (SecondaryStructureType.is(type, SecondaryStructureType.Flag.Beta)) return 'strand';
  if (
    SecondaryStructureType.is(type, SecondaryStructureType.Flag.Turn) ||
    SecondaryStructureType.is(type, SecondaryStructureType.Flag.Bend)
  )
    return 'turn';
  return null;
}

export const SecondaryStructureProvider: ProviderEntry = {
  name: 'molstar:secondary_structure',
  options: SimpleParamsSchema({ of: RequiredField(ref('complex'), 'the complex of a structure') }),
  geometry: 'ranges',
  async compute(ctx, options) {
    const of = (options as { of: string }).of;
    const parsed = ctx.source(of);
    const src = parsed?.kind === 'object' ? parsed.value : null;
    if (!StructureSource.is(src) || !ctx.set) return null;
    const structure = await structureOf(src, 0, null);
    const target: string[] = [];
    const start: number[] = [];
    const stop: number[] = [];
    const category: string[] = [];
    for (const u of complexUnits(ctx.set, of)) {
      const c = chainFor(structure, u.mvs);
      if (c?.kind !== 'polymer') continue;
      const ss = ModelSecondaryStructure.Provider.get(c.unit.model);
      if (!ss) continue;
      let open: SecondaryKind | null = null;
      let from = 0;
      let last = 0;
      const close = () => {
        if (open === null) return;
        target.push(u.ref);
        start.push(from);
        stop.push(last + 1);
        category.push(open);
        open = null;
      };
      for (let x = 0; x < c.count; x++) {
        const kind = secondaryKind(ss.type[c.residues[x]]);
        const at = c.labels[x] - u.first;
        // A run goes on while the kind does and the residues are consecutive.
        if (kind !== open || (open !== null && at !== last + 1)) {
          close();
          if (kind !== null) {
            open = kind;
            from = at;
          }
        }
        last = at;
      }
      close();
    }
    return { key: '', geometry: 'ranges', label: 'Secondary structure', columns: { target, start, stop, category } };
  },
};
