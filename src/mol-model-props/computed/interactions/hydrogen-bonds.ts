/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Structure, Unit, StructureElement } from '../../../mol-model/structure';
import { AtomGeometry, AtomGeometryAngles, calcAngles, calcPlaneAngle } from '../chemistry/geometry';
import { FeatureType, FeaturesBuilder, FeatureGroup, Features } from './features';
import { MoleculeType, ProteinBackboneAtoms } from '../../../mol-model/structure/model/types';
import { typeSymbol, bondToElementCount, bondCount, formalCharge, atomId, compId, intraConnectedTo, altLoc } from '../chemistry/util';
import { Elements } from '../../../mol-model/structure/model/properties/atomic/types';
import { InteractionType, InteractionsBuilder } from './interactions';
import { ValenceModelProvider } from '../valence-model';
import { degToRad } from '../../../mol-math/misc';

export interface HydrogenBonds {

}

export const HydrogenBondsParams = {
    maxHbondDist: PD.Numeric(3.5, { min: 1, max: 5, step: 0.1 }),
    maxHbondSulfurDist: PD.Numeric(4.1, { min: 1, max: 5, step: 0.1 }),
    maxHbondAccAngleDev: PD.Numeric(45, { min: 0, max: 180, step: 1 }, { description: 'Max deviation from ideal acceptor angle' }),
    maxHbondDonAngleDev: PD.Numeric(45, { min: 0, max: 180, step: 1 }, { description: 'Max deviation from ideal donor angle' }),
    maxHbondAccOutOfPlaneAngle: PD.Numeric(90, { min: 0, max: 180, step: 1 }),
    maxHbondDonOutOfPlaneAngle: PD.Numeric(30, { min: 0, max: 180, step: 1 }),
}
export type HydrogenBondsParams = typeof HydrogenBondsParams
export type HydrogenBondsProps = PD.Values<HydrogenBondsParams>

//

// Geometric characteristics of hydrogen bonds involving sulfur atoms in proteins
// https://doi.org/10.1002/prot.22327

// Satisfying Hydrogen Bonding Potential in Proteins (HBPLUS)
// https://doi.org/10.1006/jmbi.1994.1334
// http://www.csb.yale.edu/userguides/datamanip/hbplus/hbplus_descrip.html

function getUnitValenceModel(structure: Structure, unit: Unit.Atomic) {
    const valenceModel = ValenceModelProvider.getValue(structure).value
    if (!valenceModel) throw Error('expected valence model to be available')
    const unitValenceModel = valenceModel.get(unit.id)
    if (!unitValenceModel) throw Error('expected valence model for unit to be available')
    return unitValenceModel
}

/**
 * Potential hydrogen donor
 */
export function addUnitHydrogenDonors(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
    const { totalH } = getUnitValenceModel(structure, unit)
    const { elements, conformation } = unit
    const { x, y, z } = conformation

    for (let i = 0 as StructureElement.UnitIndex, il = elements.length; i < il; ++i) {
        const element = typeSymbol(unit, i)
        if ((
                // include both nitrogen atoms in histidine due to
                // their often ambiguous protonation assignment
                isHistidineNitrogen(unit, i)
            ) || (
                totalH[i] > 0 &&
                (element === Elements.N || element === Elements.O || element === Elements.S)
            )
        ) {
            builder.addOne(FeatureType.HydrogenDonor, FeatureGroup.None, x(elements[i]), y(elements[i]), z(elements[i]), i)
        }
    }
}

/**
 * Weak hydrogen donor.
 */
export function addUnitWeakHydrogenDonors(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
    const { totalH } = getUnitValenceModel(structure, unit)
    const { elements, conformation } = unit
    const { x, y, z } = conformation

    for (let i = 0 as StructureElement.UnitIndex, il = elements.length; i < il; ++i) {
        if (
            typeSymbol(unit, i) === Elements.C &&
            totalH[i] > 0 &&
            (
                bondToElementCount(structure, unit, i, Elements.N) > 0 ||
                bondToElementCount(structure, unit, i, Elements.O) > 0 ||
                inAromaticRingWithElectronNegativeElement(structure, unit, i)
            )
        ) {
            builder.addOne(FeatureType.WeakHydrogenDonor, FeatureGroup.None, x(elements[i]), y(elements[i]), z(elements[i]), i)
        }
    }
}

function inAromaticRingWithElectronNegativeElement(structure: Structure, unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return false // TODO
    // if (!a.isAromatic()) return false

    // const ringData = a.residueType.getRings()
    // if (!ringData) return false

    // let hasElement = false
    // const rings = ringData.rings
    // rings.forEach(ring => {
    //     if (hasElement) return  // already found one
    //     if (ring.some(idx => (a.index - a.residueAtomOffset) === idx)) {  // in ring
    //         hasElement = ring.some(idx => {
    //             const atomTypeId = a.residueType.atomTypeIdList[ idx ]
    //             const number = a.atomMap.get(atomTypeId).number
    //             return number === Elements.N || number === Elements.O
    //         })
    //     }
    // })

    // return hasElement
}

/**
 * Potential hydrogen acceptor
 */
export function addUnitHydrogenAcceptors(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
    const { charge, implicitH, idealGeometry } = getUnitValenceModel(structure, unit)
    const { elements, conformation } = unit
    const { x, y, z } = conformation

    function add(i: StructureElement.UnitIndex) {
        builder.addOne(FeatureType.HydrogenAcceptor, FeatureGroup.None, x(elements[i]), y(elements[i]), z(elements[i]), i)
    }

    for (let i = 0 as StructureElement.UnitIndex, il = elements.length; i < il; ++i) {
        const element = typeSymbol(unit, i)
        if (element === Elements.O) {
            // Basically assume all oxygen atoms are acceptors!
            add(i)
        } else if (element === Elements.N) {
            if (isHistidineNitrogen(unit, i)) {
                // include both nitrogen atoms in histidine due to
                // their often ambiguous protonation assignment
                add(i)
            } else if (charge[i] < 1) {
                // Neutral nitrogen might be an acceptor
                // It must have at least one lone pair not conjugated
                const totalBonds = bondCount(structure, unit, i) + implicitH[i]
                const ig = idealGeometry[i]
                if (
                    (ig === AtomGeometry.Tetrahedral && totalBonds < 4) ||
                    (ig === AtomGeometry.Trigonal && totalBonds < 3) ||
                    (ig === AtomGeometry.Linear && totalBonds < 2)
                ) {
                    add(i)
                }
            }
        } else if (element === Elements.S) {
            const resname = compId(unit, i)
            if (resname === 'CYS' || resname === 'MET' || formalCharge(unit, i) === -1) {
                add(i)
            }
        }
    }
}


function isWater(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return unit.model.atomicHierarchy.derived.residue.moleculeType[unit.elements[index]] === MoleculeType.Water
}

function isBackbone(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return ProteinBackboneAtoms.has(atomId(unit, index))
}

function isRing(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return unit.rings.elementRingIndices.has(index)
}

function isHistidineNitrogen(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return compId(unit, index) === 'HIS' && typeSymbol(unit, index) === Elements.N && isRing(unit, index)
}

function isBackboneHydrogenBond(unitA: Unit.Atomic, indexA: StructureElement.UnitIndex, unitB: Unit.Atomic, indexB: StructureElement.UnitIndex) {
    return isBackbone(unitA, indexA) && isBackbone(unitB, indexB)
}

function isWaterHydrogenBond(unitA: Unit.Atomic, indexA: StructureElement.UnitIndex, unitB: Unit.Atomic, indexB: StructureElement.UnitIndex) {
    return isWater(unitA, indexA) && isWater(unitB, indexB)
}

function isHydrogenBond(ti: FeatureType, tj: FeatureType) {
    return (
        (ti === FeatureType.HydrogenAcceptor && tj === FeatureType.HydrogenDonor) ||
        (ti === FeatureType.HydrogenDonor && tj === FeatureType.HydrogenAcceptor)
    )
}

function isWeakHydrogenBond(ti: FeatureType, tj: FeatureType) {
    return (
        (ti === FeatureType.WeakHydrogenDonor && tj === FeatureType.HydrogenAcceptor) ||
        (ti === FeatureType.HydrogenAcceptor && tj === FeatureType.WeakHydrogenDonor)
    )
}

function getHydrogenBondType(unitA: Unit.Atomic, indexA: StructureElement.UnitIndex, unitB: Unit.Atomic, indexB: StructureElement.UnitIndex) {
    if (isWaterHydrogenBond(unitA, indexA, unitB, indexB)) {
        return InteractionType.WaterHydrogenBond
    } else if (isBackboneHydrogenBond(unitA, indexA, unitB, indexB)) {
        return InteractionType.BackboneHydrogenBond
    } else {
        return InteractionType.HydrogenBond
    }
}

/**
 * All pairs of hydrogen donor and acceptor atoms
 */
export function addHydrogenBonds(structure: Structure, unit: Unit.Atomic, features: Features, builder: InteractionsBuilder, props: HydrogenBondsProps) {
    const { maxHbondDist, maxHbondSulfurDist, maxHbondAccAngleDev, maxHbondDonAngleDev, maxHbondAccOutOfPlaneAngle, maxHbondDonOutOfPlaneAngle } = props

    const maxAccAngleDev = degToRad(maxHbondAccAngleDev)
    const maxDonAngleDev = degToRad(maxHbondDonAngleDev)
    const maxAccOutOfPlaneAngle = degToRad(maxHbondAccOutOfPlaneAngle)
    const maxDonOutOfPlaneAngle = degToRad(maxHbondDonOutOfPlaneAngle)

    const maxDist = Math.max(maxHbondDist, maxHbondSulfurDist)
    const maxHbondDistSq = maxHbondDist * maxHbondDist

    const { x, y, z, count: n, types, offsets, members, lookup3d } = features

    const valenceModel = ValenceModelProvider.getValue(structure).value
    if (!valenceModel || !valenceModel.has(unit.id)) throw new Error('valence model required')

    const { idealGeometry } = valenceModel.get(unit.id)!

    for (let i = 0; i < n; ++i) {
        const { count, indices, squaredDistances } = lookup3d.find(x[i], y[i], z[i], maxDist)
        const ti = types[i]

        for (let r = 0; r < count; ++r) {
            const j = indices[r]
            if (j <= i) continue
            const dSq = squaredDistances[r]
            const tj = types[j]

            const isWeak = isWeakHydrogenBond(ti, tj)
            if (!isWeak && !isHydrogenBond(ti, tj)) continue

            const [ l, k ] = tj === FeatureType.HydrogenAcceptor ? [ i, j ] : [ j, i ]

            const donorIdx = members[offsets[l]]
            const acceptorIdx = members[offsets[k]]

            if (acceptorIdx === donorIdx) continue // DA to self

            const altD = altLoc(unit, donorIdx)
            const altA = altLoc(unit, acceptorIdx)

            if (altD && altA && altD !== altA) continue // incompatible alternate location id
            if (unit.residueIndex[donorIdx] === unit.residueIndex[acceptorIdx]) continue // same residue
            // check if distance is ok for non-sulfur-containing hbond
            if (typeSymbol(unit, donorIdx) !== Elements.S && typeSymbol(unit, acceptorIdx) !== Elements.S && dSq > maxHbondDistSq) continue
            // no hbond if donor and acceptor are bonded
            if (intraConnectedTo(unit, donorIdx, acceptorIdx)) continue // TODO limit to covalent bonds

            const donAngles = calcAngles(structure, unit, donorIdx, unit, acceptorIdx)
            const idealDonAngle = AtomGeometryAngles.get(idealGeometry[donorIdx]) || degToRad(120)
            if (donAngles.some(donAngle => Math.abs(idealDonAngle - donAngle) > maxDonAngleDev)) continue

            if (idealGeometry[donorIdx] === AtomGeometry.Trigonal) {
                const outOfPlane = calcPlaneAngle(structure, unit, donorIdx, unit, acceptorIdx)
                if (outOfPlane !== undefined && outOfPlane > maxDonOutOfPlaneAngle) continue
            }

            const accAngles = calcAngles(structure, unit, acceptorIdx, unit, donorIdx)
            const idealAccAngle = AtomGeometryAngles.get(idealGeometry[acceptorIdx]) || degToRad(120)
            // Do not limit large acceptor angles
            if (accAngles.some(accAngle => idealAccAngle - accAngle > maxAccAngleDev)) continue

            if (idealGeometry[acceptorIdx] === AtomGeometry.Trigonal) {
                const outOfPlane = calcPlaneAngle(structure, unit, acceptorIdx, unit, donorIdx)
                if (outOfPlane !== undefined && outOfPlane > maxAccOutOfPlaneAngle) continue
            }

            const bondType = isWeak ? InteractionType.WeakHydrogenBond : getHydrogenBondType(unit, donorIdx, unit, acceptorIdx)
            builder.add(l, k, bondType)
        }
    }
}