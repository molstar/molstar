/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 *
 * Perception of bond orders (double / triple / aromatic) from 3D coordinates.
 * Based on the following works:
 * @link https://www.daylight.com/meetings/mug01/Sayle/m4xbondage.html
 * Roger Sayle (2001) "PDB: Cruft to Content" algorithm. The perception is based on
 * bond angles measurements to determine a first set of hybridization states,
 * followed by functional groups recognition, then a step favoring the perception
 * of aromatic rings (large tolerance from planarity, but complete Huckel compliance check)
 * and performs successive kekulization steps to assign bonds. The kekulization is not detailed.
 * It is optimized against a 17 structures dataset taken from the PDB, some authors suggest
 * that the method is fitted to the benchmark.
 *
 * @link https://doi.org/10.1021/ci049915d
 * Paul Labute (2005) "On the Perception of Molecules from 3D Atomic Coordinates"
 * In this work, no pattern recognition is performed, no preference towards aromaticity
 * and no assignment of sp2/sp hybridization done based on geometry.
 * This agnostic method proceeds by elimination (unambiguous 3-dimensional neighbours,
 * bonds too long or unambiguously non-planar) to find a set of candidate sp2/sp atoms.
 * Then a maximum weighted matching algorithm is used to assign double bonds to the candidate atoms.
 * It is benchmarked against a dataset of 177 structures taken from the PDB (166 correctly assigned),
 * containing several structures that cannot be assigned correctly using the Sayle method.
 * Some authors suggest that Labute's method is more sensitive to coordinates accuracy.
 * Some conventions are used that require post-processing for standardization and make the
 * integration into a different infrastructure more challenging (e.g. definition of SP3).
 *
 * @link https://doi.org/10.1186/1758-2946-4-26
 * Qian Zhang (2012) "A rule-based algorithm for automatic bond type perception"
 * This method is based on an extensive set of rules to assign either strict or uncertain
 * bond orders, progressively refining the assignment. Rules are based on bond lengths, angles
 * and atom types.
 * It is benchmarked against the same datasets as above and obtains similar results.
 *
 * The current implementation follows most of the Sayle method, but strives to mark ambiguous
 * cases to avoid over perception. Some of the Labute's heuristics are used to mark
 * unambiguous single bonds or sp3 atoms. Some thresholds are also taken from Zhang's work.
 * Similarly to Sayle's method, the current implementation favors aromaticity perception
 * following Huckel's rule.
 * Against Labute's dataset, current implementation has similar performance as Labute's
 * and Zhang's methods.
 * Note that the perception is inherently limited given the missing information (H suppression, missing charges, etc...).
 * As this implementation favors aromaticity perception, compounds like NDP (NADPH) or FDA (FADH2)
 * will be perceived as aromatic, i.e. like their oxidized forms. Labute notes in his article
 * that the "geometric differences of FAD and FADH2 are slight (optimized structures superpose to 0.067 Å RMSD)"
 * so we introduce here a consistent bias (we are deterministically wrong!).
 */

import { Segmentation } from '@molstar/core/data/int';
import { IndexPairBonds } from '@molstar/model/formats/structure/property/bonds/index-pair';
import { IntAdjacencyGraph } from '@molstar/core/math/graph';
import type { Unit } from '@molstar/model/model/structure';
import { BondType } from '@molstar/model/model/structure/model/types';
import type { UnitIndex } from '@molstar/model/model/structure/structure/element/util';
import type { IntraUnitBondProps, IntraUnitBonds } from '@molstar/model/model/structure/structure/unit/bonds/data';
import { computeResidueRings } from '@molstar/model/model/structure/structure/unit/rings/compute';
import { applyCachedChemCompPattern, cacheChemCompPattern, chemCompHash } from './chemcomp-cache.js';
import { State, hasPerceivableBond } from './common.js';
import type { InterUnitNeighbours } from './inter-unit-neighbours.js';
import { applyFunctionalGroups } from './functional-groups-connectivity.js';
import { perceiveAromaticRings } from './huckel.js';
import { assignHybridization } from './hybridization.js';
import { assignBondOrders } from './kekule.js';
import { KnownComponents } from './known-components.js';

/**
 * - `auto`: perceive the single bonds whose order nobody stated (`Computed`, `OrderUnknown`), in residues other than
 *   the canonical ones;
 * - `forceCompute`: perceive every covalent bond of every residue, stated orders and canonical residues included.
 */
export type BondOrdersMode = 'auto' | 'forceCompute';

/**
 * The bonds of `unit` with perceived orders: a copy of `builtIn` (same topology), written in place residue by residue.
 * A function of the model, the elements and the mode only (_spec §1.3): the inter-unit neighbours and the coordinates are
 * the model's; `getInterUnit` is called only once a residue needs perceiving (it builds the model structure).
 */
export function perceiveBondOrders(
  unit: Unit.Atomic,
  builtIn: IntraUnitBonds,
  getInterUnit: () => InterUnitNeighbours,
  mode: BondOrdersMode,
): IntraUnitBonds {
  const force = mode === 'forceCompute';
  // a bond table (MOL/SDF, MOL2, ...) is the whole built-in graph and its bonds are never unknown: nothing to
  // perceive in `auto` (_spec §1.9)
  if (!force && IndexPairBonds.Provider.get(unit.model)) return builtIn;
  const { order, flags } = builtIn.edgeProps;
  const newOrder = new Int8Array(order.length);
  const newFlags = new Uint16Array(flags.length);
  for (let t = 0, tl = order.length; t < tl; t++) {
    if (force && BondType.isCovalent(flags[t])) {
      // everything re-perceived: no stated order, no stated aromaticity
      newOrder[t] = 1;
      newFlags[t] = (flags[t] | BondType.Flag.OrderUnknown) & ~BondType.Flag.Aromatic;
    } else {
      newOrder[t] = order[t];
      newFlags[t] = flags[t];
    }
  }

  // set once perception touched a residue: perceived orders depend on the geometry, so a frame with `dynamicBonds`
  // is perceived again even where the built-in bonds would carry over (_spec §1.7)
  const props: { -readonly [K in keyof IntraUnitBondProps]: IntraUnitBondProps[K] } = {
    ...builtIn.props,
    cacheable: true,
  };
  const bonds: IntraUnitBonds = IntAdjacencyGraph.create(
    builtIn.offset,
    builtIn.a,
    builtIn.b,
    builtIn.edgeCount,
    { ...builtIn.edgeProps, order: newOrder, flags: newFlags },
    props,
  );
  if (unit.elements.length <= 1) return bonds;

  const { label_comp_id } = unit.model.atomicHierarchy.atoms;
  const prefix = `${mode}|`;
  const residuesIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements);
  while (residuesIt.hasNext) {
    const { start, end } = residuesIt.move() as unknown as { start: UnitIndex; end: UnitIndex };
    if (end - start < 2) continue;

    const compId = label_comp_id.value(unit.elements[start]);
    if (!force && KnownComponents.has(compId)) continue;
    // every order inside the residue is stated: nothing to perceive (_spec §1.9)
    if (!hasPerceivableBond(unit, bonds, start, end)) continue;

    props.canRemap = false;
    const hash = chemCompHash(unit, start, end, prefix);
    if (applyCachedChemCompPattern(unit, bonds, start, end, hash)) continue;

    // the residue's rings over the graph being perceived, never `unit.rings` (`unit.bonds` is not set yet)
    const rings = computeResidueRings(unit, start, end, bonds);
    const state = State(unit, bonds, rings, getInterUnit(), start, end);
    assignHybridization(state);
    applyFunctionalGroups(state);
    perceiveAromaticRings(state);
    assignBondOrders(state);

    cacheChemCompPattern(state, hash);
    if (state.firstAltLoc) {
      applyCachedChemCompPattern(unit, bonds, start, end, hash); // apply to handle remaining altlocs (%B, %C,...)
    }
  }
  return bonds;
}
