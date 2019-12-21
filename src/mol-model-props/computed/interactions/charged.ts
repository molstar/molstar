/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Fred Ludlow <Fred.Ludlow@astx.com>
 *
 * based in part on NGL (https://github.com/arose/ngl)
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Structure, Unit, StructureElement } from '../../../mol-model/structure';
import { FeaturesBuilder, Features } from './features';
import { ProteinBackboneAtoms, PolymerNames, BaseNames, ElementSymbol } from '../../../mol-model/structure/model/types';
import { typeSymbol, atomId, altLoc, eachBondedAtom } from '../chemistry/util';
import { Elements } from '../../../mol-model/structure/model/properties/atomic/types';
import { ValenceModelProvider } from '../valence-model';
import { degToRad } from '../../../mol-math/misc';
import { FeatureType, FeatureGroup, InteractionType } from './common';
import { LinkProvider } from './links';
import { Segmentation, SortedArray } from '../../../mol-data/int';
import { isGuanidine, isAcetamidine, isPhosphate, isSulfonicAcid, isSulfate, isCarboxylate } from '../chemistry/functional-group';
import { PrincipalAxes } from '../../../mol-math/linear-algebra/matrix/principal-axes';
import { getPositions } from '../../../mol-model/structure/util';
import { Vec3 } from '../../../mol-math/linear-algebra';

export const ChargedParams = {
    piStackingDistanceMax: PD.Numeric(5.5, { min: 1, max: 8, step: 0.1 }),
    piStackingOffsetMax: PD.Numeric(2.0, { min: 0, max: 4, step: 0.1 }),
    piStackingAngleDevMax: PD.Numeric(30, { min: 0, max: 180, step: 1 }),
    cationPiDistanceMax: PD.Numeric(6.0, { min: 1, max: 8, step: 0.1 }),
    cationPiOffsetMax: PD.Numeric(2.0, { min: 0, max: 4, step: 0.1 }),
    ionicDistanceMax: PD.Numeric(5.0, { min: 0, max: 8, step: 0.1 }),
}
export type ChargedParams = typeof ChargedParams
export type ChargedProps = PD.Values<ChargedParams>

//

const PositvelyCharged = ['ARG', 'HIS', 'LYS']
const NegativelyCharged = ['GLU', 'ASP']

function getUnitValenceModel(structure: Structure, unit: Unit.Atomic) {
    const valenceModel = ValenceModelProvider.getValue(structure).value
    if (!valenceModel) throw Error('expected valence model to be available')
    const unitValenceModel = valenceModel.get(unit.id)
    if (!unitValenceModel) throw Error('expected valence model for unit to be available')
    return unitValenceModel
}

export function addUnitPositiveCharges(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
    const { charge } = getUnitValenceModel(structure, unit)
    const { elements } = unit
    const { x, y, z } = unit.model.atomicConformation

    const addedElements = new Set<StructureElement.UnitIndex>()

    const { label_comp_id } = unit.model.atomicHierarchy.residues
    const residueIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, elements)

    while (residueIt.hasNext) {
        const { index: residueIndex, start, end } = residueIt.move();
        const compId = label_comp_id.value(residueIndex)

        if (PositvelyCharged.includes(compId)) {
            builder.startState()
            for (let j = start as StructureElement.UnitIndex; j < end; ++j) {
                if (typeSymbol(unit, j) === Elements.N && !ProteinBackboneAtoms.has(atomId(unit, j))) {
                    builder.pushMember(x[elements[j]], y[elements[j]], z[elements[j]], j)
                }
            }
            builder.finishState(FeatureType.PositiveCharge, FeatureGroup.None)
        } else if (!PolymerNames.has(compId)) {
            for (let j = start as StructureElement.UnitIndex; j < end; ++j) {
                builder.startState()
                if (typeSymbol(unit, j) === Elements.N && !ProteinBackboneAtoms.has(atomId(unit, j))) {
                    builder.pushMember(x[elements[j]], y[elements[j]], z[elements[j]], j)
                }
                builder.finishState(FeatureType.PositiveCharge, FeatureGroup.None)

                let group = FeatureGroup.None
                if (isGuanidine(structure, unit, j)) {
                    group = FeatureGroup.Guanidine
                } else if (isAcetamidine(structure, unit, j)) {
                    group = FeatureGroup.Acetamidine
                }
                if (group) {
                    builder.startState()
                    eachBondedAtom(structure, unit, j, (_, k) => {
                        if (typeSymbol(unit, k) === Elements.N) {
                            addedElements.add(k)
                            builder.pushMember(x[elements[k]], y[elements[k]], z[elements[k]], k)
                        }
                    })
                    builder.finishState(FeatureType.PositiveCharge, group)
                }
            }

            for (let j = start as StructureElement.UnitIndex; j < end; ++j) {
                if (charge[j] > 0 && !addedElements.has(j)) {
                    builder.add(FeatureType.PositiveCharge, FeatureGroup.None, x[elements[j]], y[elements[j]], z[elements[j]], j)
                }
            }
        }
    }
}

export function addUnitNegativeCharges(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
    const { charge } = getUnitValenceModel(structure, unit)
    const { elements } = unit
    const { x, y, z } = unit.model.atomicConformation

    const addedElements = new Set<StructureElement.UnitIndex>()

    const { label_comp_id } = unit.model.atomicHierarchy.residues
    const residueIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, elements)

    while (residueIt.hasNext) {
        const { index: residueIndex, start, end } = residueIt.move();
        const compId = label_comp_id.value(residueIndex)

        if (NegativelyCharged.includes(compId)) {
            builder.startState()
            for (let j = start as StructureElement.UnitIndex; j < end; ++j) {
                if (typeSymbol(unit, j) === Elements.O && !ProteinBackboneAtoms.has(atomId(unit, j))) {
                    builder.pushMember(x[elements[j]], y[elements[j]], z[elements[j]], j)
                }
            }
            builder.finishState(FeatureType.NegativeCharge, FeatureGroup.None)
        } else if (BaseNames.has(compId)) {
            for (let j = start as StructureElement.UnitIndex; j < end; ++j) {
                if (isPhosphate(structure, unit, j)) {
                    builder.startState()
                    eachBondedAtom(structure, unit, j, (_, k) => {
                        if (typeSymbol(unit, k) === Elements.O) {
                            builder.pushMember(x[elements[k]], y[elements[k]], z[elements[k]], k)
                        }
                    })
                    builder.finishState(FeatureType.NegativeCharge, FeatureGroup.Phosphate)
                }
            }
        } else if (!PolymerNames.has(compId)) {
            for (let j = start as StructureElement.UnitIndex; j < end; ++j) {
                builder.startState()
                if (typeSymbol(unit, j) === Elements.N && !ProteinBackboneAtoms.has(atomId(unit, j))) {
                    builder.pushMember(x[elements[j]], y[elements[j]], z[elements[j]], j)
                }
                builder.finishState(FeatureType.PositiveCharge, FeatureGroup.None)

                let group = FeatureGroup.None
                if (isSulfonicAcid(structure, unit, j)) {
                    group = FeatureGroup.SulfonicAcid
                } else if (isPhosphate(structure, unit, j)) {
                    group = FeatureGroup.Phosphate
                } else if (isSulfate(structure, unit, j)) {
                    group = FeatureGroup.Sulfate
                } else if (isCarboxylate(structure, unit, j)) {
                    group = FeatureGroup.Carboxylate
                }
                if (group) {
                    builder.startState()
                    eachBondedAtom(structure, unit, j, (_, k) => {
                        if (typeSymbol(unit, k) === Elements.O) {
                            addedElements.add(k)
                            builder.pushMember(x[elements[k]], y[elements[k]], z[elements[k]], k)
                        }
                    })
                    builder.finishState(FeatureType.PositiveCharge, group)
                }
            }

            for (let j = start as StructureElement.UnitIndex; j < end; ++j) {
                if (charge[j] < 0 && !addedElements.has(j)) {
                    builder.add(FeatureType.NegativeCharge, FeatureGroup.None, x[elements[j]], y[elements[j]], z[elements[j]], j)
                }
            }
        }
    }
}

const AromaticRingElements = [
    Elements.B, Elements.C, Elements.N, Elements.O,
    Elements.SI, Elements.P, Elements.S,
    Elements.GE, Elements.AS,
    Elements.SN, Elements.SB,
    Elements.BI
] as ElementSymbol[]
const AromaticRingPlanarityThreshold = 0.05

function isRingAromatic(unit: Unit.Atomic, ring: SortedArray<StructureElement.UnitIndex>) {
    // TODO also check `chem_comp_bond.pdbx_aromatic_flag`
    let hasAromaticRingElement = false
    for (let i = 0, il = ring.length; i < il; ++i) {
        if (AromaticRingElements.includes(typeSymbol(unit, ring[i]))) {
            hasAromaticRingElement = true
            break
        }
    }
    if (!hasAromaticRingElement) return

    const ma = PrincipalAxes.calculateMomentsAxes(getPositions(unit, ring))
    return Vec3.magnitude(ma.dirC) < AromaticRingPlanarityThreshold
}

export function addUnitAromaticRings(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
    const { elements } = unit
    const { x, y, z } = unit.model.atomicConformation

    for (const ring of unit.rings.all) {
        if (!isRingAromatic(unit, ring)) continue

        builder.startState()
        for (let i = 0, il = ring.length; i < il; ++i) {
            const j = ring[i]
            builder.pushMember(x[elements[j]], y[elements[j]], z[elements[j]], j)
        }
        builder.finishState(FeatureType.AromaticRing, FeatureGroup.None)
    }
}

function isIonicInteraction(ti: FeatureType, tj: FeatureType) {
    return (
        (ti === FeatureType.NegativeCharge && tj === FeatureType.PositiveCharge) ||
        (ti === FeatureType.PositiveCharge && tj === FeatureType.NegativeCharge)
    )
}

function isPiStacking(ti: FeatureType, tj: FeatureType) {
    return ti === FeatureType.AromaticRing && tj === FeatureType.AromaticRing
}

function isCationPi(ti: FeatureType, tj: FeatureType) {
    return (
        (ti === FeatureType.AromaticRing && tj === FeatureType.PositiveCharge) ||
        (ti === FeatureType.PositiveCharge && tj === FeatureType.AromaticRing)
    )
}

function areFeaturesWithinDistanceSq(infoA: Features.Info, infoB: Features.Info, distanceSq: number): boolean {
    // TODO
    // const sn = atomSet1.length
    // const sm = atomSet2.length
    // for (let si = 0; si < sn; ++si) {
    //   ap1.index = atomSet1[ si ]
    //   for (let sj = 0; sj < sm; ++sj) {
    //     ap2.index = atomSet2[ sj ]
    //     if (ap1.distanceTo(ap2) <= maxDist) {
    //       return true
    //     }
    //   }
    // }
    return false
}


const tmpVecA = Vec3()
const tmpVecB = Vec3()
const tmpVecC = Vec3()
const tmpVecD = Vec3()

function getNormal(out: Vec3, info: Features.Info) {
    const { unit, feature, offsets} = info
    const { elements } = unit
    const { x, y, z } = unit.model.atomicConformation

    const i = offsets[feature]
    Vec3.set(tmpVecA, x[elements[i]], y[elements[i]], z[elements[i]])
    Vec3.set(tmpVecB, x[elements[i + 1]], y[elements[i + 1]], z[elements[i + 1]])
    Vec3.set(tmpVecC, x[elements[i + 2]], y[elements[i + 2]], z[elements[i + 2]])

    return Vec3.triangleNormal(out, tmpVecA, tmpVecB, tmpVecC)
}

const getOffset = function (infoA: Features.Info, infoB: Features.Info, normal: Vec3) {
    Vec3.set(tmpVecA, infoA.x[infoA.feature], infoA.y[infoA.feature], infoA.z[infoA.feature])
    Vec3.set(tmpVecB, infoB.x[infoB.feature], infoB.y[infoB.feature], infoB.z[infoB.feature])

    Vec3.sub(tmpVecC, tmpVecA, tmpVecB)

    Vec3.projectOnPlane(tmpVecD, tmpVecC, normal)
    Vec3.add(tmpVecD, tmpVecD, tmpVecB)
    return Vec3.distance(tmpVecD, tmpVecB)
}

function getOptions(props: ChargedProps) {
    return {
        piStackingDistanceMaxSq: props.piStackingDistanceMax * props.piStackingDistanceMax,
        piStackingOffsetMax: props.piStackingOffsetMax,
        piStackingAngleDevMax: degToRad(props.piStackingAngleDevMax),
        cationPiDistanceMaxSq: props.cationPiDistanceMax * props.cationPiDistanceMax,
        cationPiOffsetMax: props.cationPiOffsetMax,
        ionicDistanceMaxSq: props.ionicDistanceMax * props.ionicDistanceMax,
    }
}
type Options = ReturnType<typeof getOptions>

const deg180InRad = degToRad(180)
const deg90InRad = degToRad(90)

const tmpNormalA = Vec3()
const tmpNormalB = Vec3()

function testCharged(structure: Structure, infoA: Features.Info, infoB: Features.Info, distanceSq: number, opts: Options): InteractionType | undefined {
    const typeA = infoA.types[infoA.feature]
    const typeB = infoB.types[infoB.feature]

    const indexA = infoA.members[infoA.offsets[infoA.feature]]
    const indexB = infoB.members[infoB.offsets[infoB.feature]]

    if (indexA === indexB) return // to self

    const altA = altLoc(infoA.unit, indexA)
    const altB = altLoc(infoB.unit, indexB)

    if (altA && altB && altA !== altB) return // incompatible alternate location id
    if (infoA.unit.residueIndex[infoA.unit.elements[indexA]] === infoB.unit.residueIndex[infoB.unit.elements[indexB]]) return // same residue

    if (isIonicInteraction(typeA, typeB)) {
        if (areFeaturesWithinDistanceSq(infoA, infoB, opts.ionicDistanceMaxSq)) {
            return InteractionType.IonicInteraction
        }
    } else if (isPiStacking(typeA, typeB)) {
        if (distanceSq <= opts.piStackingDistanceMaxSq) {
            getNormal(tmpNormalA, infoA)
            getNormal(tmpNormalB, infoB)

            const angle = Vec3.angle(tmpNormalA, tmpNormalB)
            const offset = Math.min(getOffset(infoA, infoB, tmpNormalB), getOffset(infoB, infoA, tmpNormalA))
            if (offset <= opts.piStackingOffsetMax) {
                if (angle <= opts.piStackingAngleDevMax || angle >= deg180InRad - opts.piStackingAngleDevMax) {
                    return InteractionType.PiStacking  // parallel
                } else if (angle <= opts.piStackingAngleDevMax + deg90InRad && angle >= deg90InRad - opts.piStackingAngleDevMax) {
                    return InteractionType.PiStacking  // t-shaped
                }
            }
        }
    } else if (isCationPi(typeA, typeB)) {
        if (distanceSq <= opts.cationPiDistanceMaxSq) {
            const [infoR, infoC] = typeA === FeatureType.AromaticRing ? [infoA, infoB] : [infoB, infoA]

            getNormal(tmpNormalA, infoR)
            const offset = getOffset(infoC, infoR, tmpNormalA)
            if (offset <= opts.cationPiOffsetMax) {
                return InteractionType.CationPi
            }
        }
    }

}

export const ChargedProvider: LinkProvider<ChargedParams> = {
    name: 'charged',
    params: ChargedParams,
    createTester: (props: ChargedProps) => {
        const maxDistance = Math.max(props.ionicDistanceMax + 2, props.piStackingDistanceMax, props.cationPiDistanceMax)
        const opts = getOptions(props)
        return {
            maxDistanceSq: maxDistance * maxDistance,
            getType: (structure, infoA, infoB, distanceSq) => testCharged(structure, infoA, infoB, distanceSq, opts)
        }
    }
}