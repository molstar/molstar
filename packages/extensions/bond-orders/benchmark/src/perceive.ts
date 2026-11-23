/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { BondOrderPerceptionProvider } from '@molstar/bond-orders-extension';
import type { CifBlock } from '@molstar/io/reader/cif/data-model';
import { CCDFormat, trajectoryFromCCD } from '@molstar/model/formats/structure/mmcif';
import { Structure } from '@molstar/model/model/structure';
import type { Model } from '@molstar/model/model/structure/model';
import { bondRows, diagnostics, type OrderDiff, orderDiff } from './diagnostics.ts';
import { ccdInput } from './input.ts';
import { toMolBlock } from './sdf.ts';
import type { Tier } from './tiers.ts';

export interface Perceived {
  /** mol block with the perceived orders (tier 0: the CCD's) */
  perceived: string;
  /** mol block of the same atoms with the CCD's orders: the reference (plan §4.2); absent in tier 0 */
  reference?: string;
  /** L0: bonds whose order differs from the CCD's */
  l0: OrderDiff;
  /** plan §5.2 */
  diagnostics: string[];
}

/** A CCD component at the tier's coordinates, its orders perceived (`forceCompute`) beside the CCD's. */
export async function perceiveCcd(id: string, block: CifBlock, tier: Tier): Promise<Perceived> {
  const name = `${id} tier ${tier.name}`;
  if (!tier.perceive) {
    const structure = Structure.ofModel(await ccdModel(block, tier.coordinates));
    return { perceived: toMolBlock(structure, name), l0: { bonds: [], atoms: [] }, diagnostics: [] };
  }
  // the CCD's connectivity is kept (`chem_comp_bond`); `forceCompute` resets its orders and aromatic flags (plan §3.3)
  const model = await ccdModel(ccdInput(block, { hydrogens: !tier.hydrogens, charges: true }), tier.coordinates);
  const reference = Structure.ofModel(model);
  const perceived = Structure.ofModel(model, {
    bondProvider: BondOrderPerceptionProvider.create({ mode: 'forceCompute' }),
  });
  return {
    perceived: toMolBlock(perceived, name),
    reference: toMolBlock(reference, `${id} reference`),
    l0: orderDiff(bondRows(reference), bondRows(perceived)),
    diagnostics: diagnostics(perceived),
  };
}

async function ccdModel(block: CifBlock, coordinates: CCDFormat.CoordinateType): Promise<Model> {
  const trajectory = await trajectoryFromCCD(block).run();
  for (let i = 0; i < trajectory.frameCount; i++) {
    const model = trajectory.getFrameAtIndex(i) as Model;
    if (CCDFormat.CoordinateType.get(model) === coordinates) return model;
  }
  throw new Error(`no ${coordinates} coordinates`);
}
