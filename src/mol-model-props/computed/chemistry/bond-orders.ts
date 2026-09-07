/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
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
 * and no assignment of sp2/sp3 hybridization done based on geometry.
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

import { Segmentation } from '../../../mol-data/int';
import { hasIntraBondOrderFromTable } from '../../../mol-model/structure/model/properties/atomic/bonds';
import { BondType } from '../../../mol-model/structure/model/types';
import { Unit } from '../../../mol-model/structure/structure/unit';
import { type Structure } from '../../../mol-model/structure/structure/structure';
import { type IntraUnitBonds } from '../../../mol-model/structure/structure/unit/bonds/data';
import { IntAdjacencyGraph } from '../../../mol-math/graph';
import { type UnitIndex } from '../../../mol-model/structure/structure/element/util';
import { State } from './bond-orders/common';
import { applyFunctionalGroups } from './bond-orders/functional-groups-connectivity';
import { perceiveAromaticRings } from './bond-orders/huckel';
import { applyCachedChemCompPattern, cacheChemCompPattern } from './bond-orders/chemcomp-cache';
import { assignHybridization } from './bond-orders/hybridization';
import { assignBondOrders } from './bond-orders/kekule';

/**
 * Perceive and assign bond orders in place on `bonds` for residues of `unit` whose
 * orders are not otherwise known. Mutates the `order`/`flags` edge properties.
 */
function perceiveBondOrders(structure: Structure, unit: Unit.Atomic, bonds: IntraUnitBonds, cachePrefix = '', force = false) {
    if (unit.elements.length <= 1) return;

    const model = unit.model;
    const { label_comp_id } = model.atomicHierarchy.atoms;

    const residuesIt = Segmentation.transientSegments(model.atomicHierarchy.residueAtomSegments, unit.elements);
    while (residuesIt.hasNext) {
        const { start, end } = residuesIt.move() as unknown as { start: UnitIndex, end: UnitIndex };
        if (end - start < 2) continue;

        const compId = label_comp_id.value(unit.elements[start]);
        // Canonical residues get their orders from the built-in table
        if (!force && hasIntraBondOrderFromTable(compId)) continue;

        const applied = applyCachedChemCompPattern(unit, bonds, start, end, cachePrefix);
        if (applied) continue;

        const state = State(structure, unit, bonds, start, end);

        assignHybridization(state);
        applyFunctionalGroups(state);
        perceiveAromaticRings(state);
        assignBondOrders(state);

        cacheChemCompPattern(unit, start, end, cachePrefix, state);
        if (state.firstAltLoc) {
            applyCachedChemCompPattern(unit, bonds, start, end, cachePrefix); // apply to handle remaining altlocs (%B, %C,...)
        }
    }
}

/**
 * Perceived per-unit bond-order/flag overrides.
 * Lazy initialization, cached on `unit.transientCache`.
 * */
export interface BondOrdersValue {
    getUnit(structure: Structure, unit: Unit.Atomic): IntraUnitBonds['edgeProps'] | undefined;
}

/** Bond-order perception mode (see `BondOrderProvider` params). */
export type BondOrdersMode = 'auto' | 'model' | 'force';

/**
 * Perceive bond orders for a single atomic unit, returning override arrays.
 * `unit.bonds` is not mutated. Modes:
 *  - `auto`:  perceive only `Computed` (distance/CONECT-derived) bonds, leaving authoritative orders;
 *  - `force`: reset every covalent bond to single + `Computed`, then perceive — overriding even
 *    file/table-provided orders.
 * (`model` performs no perception and is handled by `calcBondOrders`.)
 */
export function computeUnitBondOrders(structure: Structure, unit: Unit.Atomic, mode: BondOrdersMode): IntraUnitBonds['edgeProps'] {
    const src = unit.bonds;
    const newOrder: number[] = [];
    const newFlags: number[] = [];
    const { order, flags } = src.edgeProps;

    for (let i = 0, il = order.length; i < il; i++) {
        if (mode === 'force' && BondType.isCovalent(flags[i])) {
            newOrder[i] = 1;
            newFlags[i] = flags[i] | BondType.Flag.Computed;
        } else {
            newOrder[i] = order[i];
            newFlags[i] = flags[i];
        }
    }
    const edgeProps = { ...src.edgeProps, order: newOrder, flags: newFlags };
    const bonds = IntAdjacencyGraph.create(src.offset, src.a, src.b, src.edgeCount, edgeProps, src.props) as IntraUnitBonds;
    perceiveBondOrders(structure, unit, bonds, mode + '|', mode === 'force');
    return bonds.edgeProps;
}

export function calcBondOrders(mode: BondOrdersMode = 'auto'): BondOrdersValue {
    const transientKey = `bond-orders-${mode}`;
    return {
        getUnit(structure: Structure, unit: Unit.Atomic) {
            if (mode === 'model') return undefined;
            let ep = unit.transientCache.get(transientKey) as IntraUnitBonds['edgeProps'] | undefined;
            if (!ep) {
                ep = computeUnitBondOrders(structure, unit, mode);
                unit.transientCache.set(transientKey, ep);
            }
            return ep;
        }
    };
}
