/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Structure, Unit, StructureElement } from '../../../mol-model/structure';
import { AtomGeometry, AtomGeometryAngles, calcAngles, calcPlaneAngle } from '../chemistry/geometry';
import { FeaturesBuilder, Features } from './features';
import { MoleculeType, ProteinBackboneAtoms } from '../../../mol-model/structure/model/types';
import { typeSymbol, bondToElementCount, bondCount, formalCharge, atomId, compId, altLoc, connectedTo } from '../chemistry/util';
import { Elements } from '../../../mol-model/structure/model/properties/atomic/types';
import { ValenceModelProvider } from '../valence-model';
import { degToRad } from '../../../mol-math/misc';
import { FeatureType, FeatureGroup, InteractionType } from './common';
import { IntraLinksBuilder, InterLinksBuilder } from './builder';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';

export interface HydrogenBonds {

}

export const HydrogenBondsParams = {
    maxDist: PD.Numeric(3.5, { min: 1, max: 5, step: 0.1 }),
    maxSulfurDist: PD.Numeric(4.1, { min: 1, max: 5, step: 0.1 }),
    maxAccAngleDev: PD.Numeric(45, { min: 0, max: 180, step: 1 }, { description: 'Max deviation from ideal acceptor angle' }),
    maxDonAngleDev: PD.Numeric(45, { min: 0, max: 180, step: 1 }, { description: 'Max deviation from ideal donor angle' }),
    maxAccOutOfPlaneAngle: PD.Numeric(90, { min: 0, max: 180, step: 1 }),
    maxDonOutOfPlaneAngle: PD.Numeric(45, { min: 0, max: 180, step: 1 }),
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
    const { elements } = unit
    const { x, y, z } = unit.model.atomicConformation

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
            builder.addOne(FeatureType.HydrogenDonor, FeatureGroup.None, x[elements[i]], y[elements[i]], z[elements[i]], i)
        }
    }
}

/**
 * Weak hydrogen donor.
 */
export function addUnitWeakHydrogenDonors(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
    const { totalH } = getUnitValenceModel(structure, unit)
    const { elements } = unit
    const { x, y, z } = unit.model.atomicConformation

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
            builder.addOne(FeatureType.WeakHydrogenDonor, FeatureGroup.None, x[elements[i]], y[elements[i]], z[elements[i]], i)
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
    const { elements } = unit
    const { x, y, z } = unit.model.atomicConformation

    function add(i: StructureElement.UnitIndex) {
        builder.addOne(FeatureType.HydrogenAcceptor, FeatureGroup.None, x[elements[i]], y[elements[i]], z[elements[i]], i)
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

interface Info {
    unit: Unit.Atomic,
    types: ArrayLike<FeatureType>,
    feature: number,
    members: ArrayLike<StructureElement.UnitIndex>,
    offsets: ArrayLike<number>,
    idealGeometry: Int8Array
}
function Info(structure: Structure, unit: Unit.Atomic, features: Features) {
    const valenceModel = ValenceModelProvider.getValue(structure).value
    if (!valenceModel || !valenceModel.has(unit.id)) throw new Error('valence model required')

    return {
        unit,
        types: features.types,
        members: features.members,
        offsets: features.offsets,
        idealGeometry: valenceModel.get(unit.id)!.idealGeometry
    } as Info
}

function getOptions(props: HydrogenBondsProps) {
    return {
        maxAccAngleDev: degToRad(props.maxAccAngleDev),
        maxDonAngleDev: degToRad(props.maxDonAngleDev),
        maxAccOutOfPlaneAngle: degToRad(props.maxAccOutOfPlaneAngle),
        maxDonOutOfPlaneAngle: degToRad(props.maxDonOutOfPlaneAngle),
        maxDist: Math.max(props.maxDist, props.maxSulfurDist),
        maxHbondDistSq: props.maxDist * props.maxDist,
    }
}
type Options = ReturnType<typeof getOptions>

function testHydrogenBond(dSq: number, structure: Structure, infoA: Info, infoB: Info, opts: Options): InteractionType | undefined {

    const typeA = infoA.types[infoA.feature]
    const typeB = infoB.types[infoB.feature]

    const isWeak = isWeakHydrogenBond(typeA, typeB)
    if (!isWeak && !isHydrogenBond(typeA, typeB)) return

    const [don, acc] = typeA === FeatureType.HydrogenAcceptor ? [infoA, infoB] : [infoB, infoA]

    const donIndex = don.members[don.offsets[don.feature]]
    const accIndex = acc.members[acc.offsets[acc.feature]]

    if (accIndex === donIndex) return // DA to self

    const altD = altLoc(don.unit, donIndex)
    const altA = altLoc(acc.unit, accIndex)

    if (altD && altA && altD !== altA) return // incompatible alternate location id
    if (don.unit.residueIndex[donIndex] === acc.unit.residueIndex[accIndex]) return // same residue

    // check if distance is ok for non-sulfur-containing hbond
    if (typeSymbol(don.unit, donIndex) !== Elements.S && typeSymbol(acc.unit, accIndex) !== Elements.S && dSq > opts.maxHbondDistSq) return

    // no hbond if donor and acceptor are bonded
    if (connectedTo(structure, don.unit, donIndex, acc.unit, accIndex)) return

    const donAngles = calcAngles(structure, don.unit, donIndex, acc.unit, accIndex)
    const idealDonAngle = AtomGeometryAngles.get(don.idealGeometry[donIndex]) || degToRad(120)
    if (donAngles.some(donAngle => Math.abs(idealDonAngle - donAngle) > opts.maxDonAngleDev)) return

    if (don.idealGeometry[donIndex] === AtomGeometry.Trigonal) {
        const outOfPlane = calcPlaneAngle(structure, don.unit, donIndex, acc.unit, accIndex)
        if (outOfPlane !== undefined && outOfPlane > opts.maxDonOutOfPlaneAngle) return
    }

    const accAngles = calcAngles(structure, acc.unit, accIndex, don.unit, donIndex)
    const idealAccAngle = AtomGeometryAngles.get(acc.idealGeometry[accIndex]) || degToRad(120)

    // Do not limit large acceptor angles
    if (accAngles.some(accAngle => idealAccAngle - accAngle > opts.maxAccAngleDev)) return

    if (acc.idealGeometry[accIndex] === AtomGeometry.Trigonal) {
        const outOfPlane = calcPlaneAngle(structure, acc.unit, accIndex, don.unit, donIndex)
        if (outOfPlane !== undefined && outOfPlane > opts.maxAccOutOfPlaneAngle) return
    }

    return isWeak ? InteractionType.WeakHydrogenBond : getHydrogenBondType(don.unit, donIndex, acc.unit, accIndex)
}

/**
 * All intra-unit pairs of hydrogen donor and acceptor atoms
 */
export function addUnitHydrogenBonds(structure: Structure, unit: Unit.Atomic, features: Features, builder: IntraLinksBuilder, props: HydrogenBondsProps) {
    const opts = getOptions(props)
    const { x, y, z, count, lookup3d } = features

    const infoA = Info(structure, unit, features)
    const infoB = { ...infoA }

    for (let i = 0; i < count; ++i) {
        const { count, indices, squaredDistances } = lookup3d.find(x[i], y[i], z[i], opts.maxDist)
        if (count === 0) continue

        infoA.feature = i

        for (let r = 0; r < count; ++r) {
            const j = indices[r]
            if (j <= i) continue

            infoB.feature = j
            const bondType = testHydrogenBond(squaredDistances[r], structure, infoA, infoB, opts)
            if (bondType) builder.add(i, j, bondType)
        }
    }
}

//

const _imageTransform = Mat4()

/**
 * All inter-unit pairs of hydrogen donor and acceptor atoms
 */
export function addStructureHydrogenBonds(structure: Structure, unitA: Unit.Atomic, featuresA: Features, unitB: Unit.Atomic, featuresB: Features, builder: InterLinksBuilder, props: HydrogenBondsProps) {
    const opts = getOptions(props)

    const { count, x: xA, y: yA, z: zA } = featuresA;
    const { lookup3d } = featuresB;

    // the lookup queries need to happen in the "unitB space".
    // that means imageA = inverseOperB(operA(i))
    const imageTransform = Mat4.mul(_imageTransform, unitB.conformation.operator.inverse, unitA.conformation.operator.matrix)
    const isNotIdentity = !Mat4.isIdentity(imageTransform)
    const imageA = Vec3()

    const { center: bCenter, radius: bRadius } = lookup3d.boundary.sphere;
    const testDistanceSq = (bRadius + opts.maxDist) * (bRadius + opts.maxDist);

    const infoA = Info(structure, unitA, featuresA)
    const infoB = Info(structure, unitB, featuresB)

    builder.startUnitPair(unitA, unitB)

    for (let i = 0; i < count; ++i) {
        Vec3.set(imageA, xA[i], yA[i], zA[i])
        if (isNotIdentity) Vec3.transformMat4(imageA, imageA, imageTransform)
        if (Vec3.squaredDistance(imageA, bCenter) > testDistanceSq) continue

        const { indices, count, squaredDistances } = lookup3d.find(imageA[0], imageA[1], imageA[2], opts.maxDist)
        if (count === 0) continue

        infoA.feature = i

        for (let r = 0; r < count; ++r) {
            const j = indices[r]
            if (j <= i) continue
            infoB.feature = j
            const bondType = testHydrogenBond(squaredDistances[r], structure, infoA, infoB, opts)
            if (bondType) builder.add(i, j, bondType)
        }
    }

    builder.finishUnitPair()
}