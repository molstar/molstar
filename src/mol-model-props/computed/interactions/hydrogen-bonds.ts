/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Fred Ludlow <Fred.Ludlow@astx.com>
 *
 * based in part on NGL (https://github.com/arose/ngl)
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Structure, Unit, StructureElement } from '../../../mol-model/structure';
import { AtomGeometry, AtomGeometryAngles, calcAngles, calcPlaneAngle } from '../chemistry/geometry';
import { FeaturesBuilder, Features } from './features';
import { typeSymbol, bondToElementCount, bondCount, formalCharge, compId, atomId } from '../chemistry/util';
import { Elements } from '../../../mol-model/structure/model/properties/atomic/types';
import { ValenceModelProvider } from '../valence-model';
import { degToRad } from '../../../mol-math/misc';
import { FeatureType, FeatureGroup, InteractionType } from './common';
import { ContactProvider } from './contacts';
import { MoleculeType, ProteinBackboneAtoms } from '../../../mol-model/structure/model/types';

const GeometryParams = {
    distanceMax: PD.Numeric(3.5, { min: 1, max: 5, step: 0.1 }),
    backbone: PD.Boolean(true, { description: 'Include backbone-to-backbone hydrogen bonds' }),
    accAngleDevMax: PD.Numeric(45, { min: 0, max: 180, step: 1 }, { description: 'Max deviation from ideal acceptor angle' }),
    donAngleDevMax: PD.Numeric(45, { min: 0, max: 180, step: 1 }, { description: 'Max deviation from ideal donor angle' }),
    accOutOfPlaneAngleMax: PD.Numeric(90, { min: 0, max: 180, step: 1 }),
    donOutOfPlaneAngleMax: PD.Numeric(45, { min: 0, max: 180, step: 1 }),
};
type GeometryParams = typeof GeometryParams
type GeometryProps = PD.Values<GeometryParams>

const HydrogenBondsParams = {
    ...GeometryParams,
    water: PD.Boolean(false, { description: 'Include water-to-water hydrogen bonds' }),
    sulfurDistanceMax: PD.Numeric(4.1, { min: 1, max: 5, step: 0.1 }),
};
type HydrogenBondsParams = typeof HydrogenBondsParams
type HydrogenBondsProps = PD.Values<HydrogenBondsParams>

const WeakHydrogenBondsParams = {
    ...GeometryParams,
};
type WeakHydrogenBondsParams = typeof WeakHydrogenBondsParams
type WeakHydrogenBondsProps = PD.Values<WeakHydrogenBondsParams>

//

// Geometric characteristics of hydrogen bonds involving sulfur atoms in proteins
// https://doi.org/10.1002/prot.22327

// Satisfying Hydrogen Bonding Potential in Proteins (HBPLUS)
// https://doi.org/10.1006/jmbi.1994.1334
// http://www.csb.yale.edu/userguides/datamanip/hbplus/hbplus_descrip.html

function getUnitValenceModel(structure: Structure, unit: Unit.Atomic) {
    const valenceModel = ValenceModelProvider.get(structure).value;
    if (!valenceModel) throw Error('expected valence model to be available');
    const unitValenceModel = valenceModel.get(unit.id);
    if (!unitValenceModel) throw Error('expected valence model for unit to be available');
    return unitValenceModel;
}

/**
 * Potential hydrogen donor
 */
function addUnitHydrogenDonors(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
    const { totalH } = getUnitValenceModel(structure, unit);
    const { elements } = unit;
    const { x, y, z } = unit.model.atomicConformation;
    const { elementAromaticRingIndices } = unit.rings;

    for (let i = 0 as StructureElement.UnitIndex, il = elements.length; i < il; ++i) {
        if (elementAromaticRingIndices.has(i)) continue;

        const element = typeSymbol(unit, i);
        if ((
            // include both nitrogen atoms in histidine due to
            // their often ambiguous protonation assignment
            isHistidineNitrogen(unit, i)
        ) || (
            totalH[i] > 0 &&
            (element === Elements.N || element === Elements.O || element === Elements.S)
        )) {
            builder.add(FeatureType.HydrogenDonor, FeatureGroup.None, x[elements[i]], y[elements[i]], z[elements[i]], i);
        }
    }
}

/**
 * Weak hydrogen donor.
 */
function addUnitWeakHydrogenDonors(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
    const { totalH } = getUnitValenceModel(structure, unit);
    const { elements } = unit;
    const { x, y, z } = unit.model.atomicConformation;

    for (let i = 0 as StructureElement.UnitIndex, il = elements.length; i < il; ++i) {
        if (
            typeSymbol(unit, i) === Elements.C &&
            totalH[i] > 0 &&
            (
                bondToElementCount(structure, unit, i, Elements.N) > 0 ||
                bondToElementCount(structure, unit, i, Elements.O) > 0 ||
                inAromaticRingWithElectronNegativeElement(unit, i)
            )
        ) {
            builder.add(FeatureType.WeakHydrogenDonor, FeatureGroup.None, x[elements[i]], y[elements[i]], z[elements[i]], i);
        }
    }
}

function inAromaticRingWithElectronNegativeElement(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    const { elementAromaticRingIndices, all } = unit.rings;
    const ringIndices = elementAromaticRingIndices.get(index);
    if (ringIndices === undefined) return false;

    for (let i = 0, il = ringIndices.length; i < il; ++i) {
        const ring = all[ringIndices[i]];
        for (let j = 0, jl = ring.length; j < jl; ++j) {
            const element = typeSymbol(unit, ring[j]);
            if (element === Elements.N || element === Elements.O) {
                return true;
            }
        }
    }
    return false;
}

/**
 * Potential hydrogen acceptor
 */
function addUnitHydrogenAcceptors(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
    const { charge, implicitH, idealGeometry } = getUnitValenceModel(structure, unit);
    const { elements } = unit;
    const { x, y, z } = unit.model.atomicConformation;
    const { elementAromaticRingIndices } = unit.rings;

    const add = (i: StructureElement.UnitIndex) => {
        builder.add(FeatureType.HydrogenAcceptor, FeatureGroup.None, x[elements[i]], y[elements[i]], z[elements[i]], i);
    };

    for (let i = 0 as StructureElement.UnitIndex, il = elements.length; i < il; ++i) {
        if (elementAromaticRingIndices.has(i)) continue;

        const element = typeSymbol(unit, i);
        if (element === Elements.O) {
            // Basically assume all oxygen atoms are acceptors!
            add(i);
        } else if (element === Elements.N) {
            if (isHistidineNitrogen(unit, i)) {
                // include both nitrogen atoms in histidine due to
                // their often ambiguous protonation assignment
                add(i);
            } else if (charge[i] < 1) {
                // Neutral nitrogen might be an acceptor
                // It must have at least one lone pair not conjugated
                const totalBonds = bondCount(structure, unit, i) + implicitH[i];
                const ig = idealGeometry[i];
                if (
                    (ig === AtomGeometry.Tetrahedral && totalBonds < 4) ||
                    (ig === AtomGeometry.Trigonal && totalBonds < 3) ||
                    (ig === AtomGeometry.Linear && totalBonds < 2)
                ) {
                    add(i);
                }
            }
        } else if (element === Elements.S) {
            const resname = compId(unit, i);
            if (resname === 'CYS' || resname === 'MET' || formalCharge(unit, i) === -1) {
                add(i);
            }
        }
    }
}

function isWater(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return unit.model.atomicHierarchy.derived.residue.moleculeType[unit.residueIndex[unit.elements[index]]] === MoleculeType.Water;
}

function isBackbone(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return ProteinBackboneAtoms.has(atomId(unit, index));
}

function isRing(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return unit.rings.elementRingIndices.has(index);
}

function isHistidineNitrogen(unit: Unit.Atomic, index: StructureElement.UnitIndex) {
    return compId(unit, index) === 'HIS' && typeSymbol(unit, index) === Elements.N && isRing(unit, index);
}

function isBackboneHydrogenBond(unitA: Unit.Atomic, indexA: StructureElement.UnitIndex, unitB: Unit.Atomic, indexB: StructureElement.UnitIndex) {
    return isBackbone(unitA, indexA) && isBackbone(unitB, indexB);
}

function isWaterHydrogenBond(unitA: Unit.Atomic, indexA: StructureElement.UnitIndex, unitB: Unit.Atomic, indexB: StructureElement.UnitIndex) {
    return isWater(unitA, indexA) && isWater(unitB, indexB);
}

function isHydrogenBond(ti: FeatureType, tj: FeatureType) {
    return (
        (ti === FeatureType.HydrogenAcceptor && tj === FeatureType.HydrogenDonor) ||
        (ti === FeatureType.HydrogenDonor && tj === FeatureType.HydrogenAcceptor)
    );
}

function isWeakHydrogenBond(ti: FeatureType, tj: FeatureType) {
    return (
        (ti === FeatureType.WeakHydrogenDonor && tj === FeatureType.HydrogenAcceptor) ||
        (ti === FeatureType.HydrogenAcceptor && tj === FeatureType.WeakHydrogenDonor)
    );
}

function getGeometryOptions(props: GeometryProps) {
    return {
        includeBackbone: props.backbone,
        maxAccAngleDev: degToRad(props.accAngleDevMax),
        maxDonAngleDev: degToRad(props.donAngleDevMax),
        maxAccOutOfPlaneAngle: degToRad(props.accOutOfPlaneAngleMax),
        maxDonOutOfPlaneAngle: degToRad(props.donOutOfPlaneAngleMax),
    };
}
type GeometryOptions = ReturnType<typeof getGeometryOptions>

function getHydrogenBondsOptions(props: HydrogenBondsProps) {
    return {
        ...getGeometryOptions(props),
        includeWater: props.water,
        maxSulfurDistSq: props.sulfurDistanceMax *  props.sulfurDistanceMax,
        maxDistSq: props.distanceMax *  props.distanceMax
    };
}
type HydrogenBondsOptions = ReturnType<typeof getHydrogenBondsOptions>

const deg120InRad = degToRad(120);

function checkGeometry(structure: Structure, don: Features.Info, acc: Features.Info, opts: GeometryOptions): true | undefined {

    const donIndex = don.members[don.offsets[don.feature]];
    const accIndex = acc.members[acc.offsets[acc.feature]];

    if (!opts.includeBackbone && isBackboneHydrogenBond(don.unit, donIndex, acc.unit, accIndex)) return;

    const donAngles = calcAngles(structure, don.unit, donIndex, acc.unit, accIndex);
    const idealDonAngle = AtomGeometryAngles.get(don.idealGeometry[donIndex]) || deg120InRad;
    if (donAngles.some(donAngle => Math.abs(idealDonAngle - donAngle) > opts.maxDonAngleDev)) return;

    if (don.idealGeometry[donIndex] === AtomGeometry.Trigonal) {
        const outOfPlane = calcPlaneAngle(structure, don.unit, donIndex, acc.unit, accIndex);
        if (outOfPlane !== undefined && outOfPlane > opts.maxDonOutOfPlaneAngle) return;
    }

    const accAngles = calcAngles(structure, acc.unit, accIndex, don.unit, donIndex);
    const idealAccAngle = AtomGeometryAngles.get(acc.idealGeometry[accIndex]) || deg120InRad;

    // Do not limit large acceptor angles
    if (accAngles.some(accAngle => idealAccAngle - accAngle > opts.maxAccAngleDev)) return;

    if (acc.idealGeometry[accIndex] === AtomGeometry.Trigonal) {
        const outOfPlane = calcPlaneAngle(structure, acc.unit, accIndex, don.unit, donIndex);
        if (outOfPlane !== undefined && outOfPlane > opts.maxAccOutOfPlaneAngle) return;
    }

    return true;
}

function testHydrogenBond(structure: Structure, infoA: Features.Info, infoB: Features.Info, distanceSq: number, opts: HydrogenBondsOptions): InteractionType | undefined {
    const typeA = infoA.types[infoA.feature];
    const typeB = infoB.types[infoB.feature];

    if (!isHydrogenBond(typeA, typeB)) return;

    const [don, acc] = typeB === FeatureType.HydrogenAcceptor ? [infoA, infoB] : [infoB, infoA];

    const donIndex = don.members[don.offsets[don.feature]];
    const accIndex = acc.members[acc.offsets[acc.feature]];

    // check if distance is ok depending on non-sulfur-containing hbond
    const maxDistSq = typeSymbol(don.unit, donIndex) === Elements.S || typeSymbol(acc.unit, accIndex) === Elements.S ? opts.maxSulfurDistSq : opts.maxDistSq;
    if (distanceSq > maxDistSq) return;

    if (!opts.includeWater && isWaterHydrogenBond(don.unit, donIndex, acc.unit, accIndex)) return;

    if (!checkGeometry(structure, don, acc, opts)) return;

    return InteractionType.HydrogenBond;
}

function testWeakHydrogenBond(structure: Structure, infoA: Features.Info, infoB: Features.Info, distanceSq: number, opts: GeometryOptions): InteractionType | undefined {
    const typeA = infoA.types[infoA.feature];
    const typeB = infoB.types[infoB.feature];

    if (!isWeakHydrogenBond(typeA, typeB)) return;

    const [don, acc] = typeB === FeatureType.HydrogenAcceptor ? [infoA, infoB] : [infoB, infoA];

    if (!checkGeometry(structure, don, acc, opts)) return;

    return InteractionType.WeakHydrogenBond;
}

//

export const HydrogenDonorProvider = Features.Provider([FeatureType.HydrogenDonor], addUnitHydrogenDonors);
export const WeakHydrogenDonorProvider = Features.Provider([FeatureType.WeakHydrogenDonor], addUnitWeakHydrogenDonors);
export const HydrogenAcceptorProvider = Features.Provider([FeatureType.HydrogenAcceptor], addUnitHydrogenAcceptors);

export const HydrogenBondsProvider: ContactProvider<HydrogenBondsParams> = {
    name: 'hydrogen-bonds',
    params: HydrogenBondsParams,
    createTester: (props: HydrogenBondsProps) => {
        const maxDistance = Math.max(props.distanceMax, props.sulfurDistanceMax);
        const opts = getHydrogenBondsOptions(props);
        return {
            maxDistance,
            requiredFeatures: new Set([FeatureType.HydrogenDonor, FeatureType.HydrogenAcceptor]),
            getType: (structure, infoA, infoB, distanceSq) => testHydrogenBond(structure, infoA, infoB, distanceSq, opts)
        };
    }
};

export const WeakHydrogenBondsProvider: ContactProvider<WeakHydrogenBondsParams> = {
    name: 'weak-hydrogen-bonds',
    params: WeakHydrogenBondsParams,
    createTester: (props: WeakHydrogenBondsProps) => {
        const opts = getGeometryOptions(props);
        return {
            maxDistance: props.distanceMax,
            requiredFeatures: new Set([FeatureType.WeakHydrogenDonor, FeatureType.HydrogenAcceptor]),
            getType: (structure, infoA, infoB, distanceSq) => testWeakHydrogenBond(structure, infoA, infoB, distanceSq, opts)
        };
    }
};