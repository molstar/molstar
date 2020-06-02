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
import { ProteinBackboneAtoms, PolymerNames, BaseNames } from '../../../mol-model/structure/model/types';
import { typeSymbol, atomId, eachBondedAtom } from '../chemistry/util';
import { Elements } from '../../../mol-model/structure/model/properties/atomic/types';
import { ValenceModelProvider } from '../valence-model';
import { degToRad } from '../../../mol-math/misc';
import { FeatureType, FeatureGroup, InteractionType } from './common';
import { ContactProvider } from './contacts';
import { Segmentation } from '../../../mol-data/int';
import { isGuanidine, isAcetamidine, isPhosphate, isSulfonicAcid, isSulfate, isCarboxylate } from '../chemistry/functional-group';
import { Vec3 } from '../../../mol-math/linear-algebra';

const IonicParams = {
    distanceMax: PD.Numeric(5.0, { min: 0, max: 8, step: 0.1 }),
};
type IonicParams = typeof IonicParams
type IonicProps = PD.Values<IonicParams>

const PiStackingParams = {
    distanceMax: PD.Numeric(5.5, { min: 1, max: 8, step: 0.1 }),
    offsetMax: PD.Numeric(2.0, { min: 0, max: 4, step: 0.1 }),
    angleDevMax: PD.Numeric(30, { min: 0, max: 180, step: 1 }),
};
type PiStackingParams = typeof PiStackingParams
type PiStackingProps = PD.Values<PiStackingParams>

const CationPiParams = {
    distanceMax: PD.Numeric(6.0, { min: 1, max: 8, step: 0.1 }),
    offsetMax: PD.Numeric(2.0, { min: 0, max: 4, step: 0.1 }),
};
type CationPiParams = typeof CationPiParams
type CationPiProps = PD.Values<CationPiParams>

//

const PositvelyCharged = ['ARG', 'HIS', 'LYS'];
const NegativelyCharged = ['GLU', 'ASP'];

function getUnitValenceModel(structure: Structure, unit: Unit.Atomic) {
    const valenceModel = ValenceModelProvider.get(structure).value;
    if (!valenceModel) throw Error('expected valence model to be available');
    const unitValenceModel = valenceModel.get(unit.id);
    if (!unitValenceModel) throw Error('expected valence model for unit to be available');
    return unitValenceModel;
}

function addUnitPositiveCharges(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
    const { charge } = getUnitValenceModel(structure, unit);
    const { elements } = unit;
    const { x, y, z } = unit.model.atomicConformation;

    const addedElements = new Set<StructureElement.UnitIndex>();

    const { label_comp_id } = unit.model.atomicHierarchy.atoms;
    const residueIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, elements);

    while (residueIt.hasNext) {
        const { index: residueIndex, start, end } = residueIt.move();
        const compId = label_comp_id.value(unit.model.atomicHierarchy.residueAtomSegments.offsets[residueIndex]);

        if (PositvelyCharged.includes(compId)) {
            builder.startState();
            for (let j = start as StructureElement.UnitIndex; j < end; ++j) {
                if (typeSymbol(unit, j) === Elements.N && !ProteinBackboneAtoms.has(atomId(unit, j))) {
                    builder.pushMember(x[elements[j]], y[elements[j]], z[elements[j]], j);
                }
            }
            builder.finishState(FeatureType.PositiveCharge, FeatureGroup.None);
        } else if (!PolymerNames.has(compId)) {
            addedElements.clear();

            for (let j = start as StructureElement.UnitIndex; j < end; ++j) {
                let group = FeatureGroup.None;
                if (isGuanidine(structure, unit, j)) {
                    group = FeatureGroup.Guanidine;
                } else if (isAcetamidine(structure, unit, j)) {
                    group = FeatureGroup.Acetamidine;
                }
                if (group) {
                    builder.startState();
                    eachBondedAtom(structure, unit, j, (_, k) => {
                        if (typeSymbol(unit, k) === Elements.N) {
                            addedElements.add(k);
                            builder.pushMember(x[elements[k]], y[elements[k]], z[elements[k]], k);
                        }
                    });
                    builder.finishState(FeatureType.PositiveCharge, group);
                }
            }

            for (let j = start as StructureElement.UnitIndex; j < end; ++j) {
                if (charge[j] > 0 && !addedElements.has(j)) {
                    builder.add(FeatureType.PositiveCharge, FeatureGroup.None, x[elements[j]], y[elements[j]], z[elements[j]], j);
                }
            }
        }
    }
}

function addUnitNegativeCharges(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
    const { charge } = getUnitValenceModel(structure, unit);
    const { elements } = unit;
    const { x, y, z } = unit.model.atomicConformation;

    const addedElements = new Set<StructureElement.UnitIndex>();

    const { label_comp_id } = unit.model.atomicHierarchy.atoms;
    const residueIt = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, elements);

    while (residueIt.hasNext) {
        const { index: residueIndex, start, end } = residueIt.move();
        const compId = label_comp_id.value(unit.model.atomicHierarchy.residueAtomSegments.offsets[residueIndex]);

        if (NegativelyCharged.includes(compId)) {
            builder.startState();
            for (let j = start as StructureElement.UnitIndex; j < end; ++j) {
                if (typeSymbol(unit, j) === Elements.O && !ProteinBackboneAtoms.has(atomId(unit, j))) {
                    builder.pushMember(x[elements[j]], y[elements[j]], z[elements[j]], j);
                }
            }
            builder.finishState(FeatureType.NegativeCharge, FeatureGroup.None);
        } else if (BaseNames.has(compId)) {
            for (let j = start as StructureElement.UnitIndex; j < end; ++j) {
                if (isPhosphate(structure, unit, j)) {
                    builder.startState();
                    eachBondedAtom(structure, unit, j, (_, k) => {
                        if (typeSymbol(unit, k) === Elements.O) {
                            builder.pushMember(x[elements[k]], y[elements[k]], z[elements[k]], k);
                        }
                    });
                    builder.finishState(FeatureType.NegativeCharge, FeatureGroup.Phosphate);
                }
            }
        } else if (!PolymerNames.has(compId)) {
            for (let j = start as StructureElement.UnitIndex; j < end; ++j) {
                builder.startState();
                if (typeSymbol(unit, j) === Elements.N && !ProteinBackboneAtoms.has(atomId(unit, j))) {
                    builder.pushMember(x[elements[j]], y[elements[j]], z[elements[j]], j);
                }
                builder.finishState(FeatureType.NegativeCharge, FeatureGroup.None);

                let group = FeatureGroup.None;
                if (isSulfonicAcid(structure, unit, j)) {
                    group = FeatureGroup.SulfonicAcid;
                } else if (isPhosphate(structure, unit, j)) {
                    group = FeatureGroup.Phosphate;
                } else if (isSulfate(structure, unit, j)) {
                    group = FeatureGroup.Sulfate;
                } else if (isCarboxylate(structure, unit, j)) {
                    group = FeatureGroup.Carboxylate;
                }
                if (group) {
                    builder.startState();
                    eachBondedAtom(structure, unit, j, (_, k) => {
                        if (typeSymbol(unit, k) === Elements.O) {
                            addedElements.add(k);
                            builder.pushMember(x[elements[k]], y[elements[k]], z[elements[k]], k);
                        }
                    });
                    builder.finishState(FeatureType.NegativeCharge, group);
                }
            }

            for (let j = start as StructureElement.UnitIndex; j < end; ++j) {
                if (charge[j] < 0 && !addedElements.has(j)) {
                    builder.add(FeatureType.NegativeCharge, FeatureGroup.None, x[elements[j]], y[elements[j]], z[elements[j]], j);
                }
            }
        }
    }
}

function addUnitAromaticRings(structure: Structure, unit: Unit.Atomic, builder: FeaturesBuilder) {
    const { elements } = unit;
    const { x, y, z } = unit.model.atomicConformation;

    for (const ringIndex of unit.rings.aromaticRings) {
        const ring = unit.rings.all[ringIndex];
        builder.startState();
        for (let i = 0, il = ring.length; i < il; ++i) {
            const j = ring[i];
            builder.pushMember(x[elements[j]], y[elements[j]], z[elements[j]], j);
        }
        builder.finishState(FeatureType.AromaticRing, FeatureGroup.None);
    }
}

function isIonic(ti: FeatureType, tj: FeatureType) {
    return (
        (ti === FeatureType.NegativeCharge && tj === FeatureType.PositiveCharge) ||
        (ti === FeatureType.PositiveCharge && tj === FeatureType.NegativeCharge)
    );
}

function isPiStacking(ti: FeatureType, tj: FeatureType) {
    return ti === FeatureType.AromaticRing && tj === FeatureType.AromaticRing;
}

function isCationPi(ti: FeatureType, tj: FeatureType) {
    return (
        (ti === FeatureType.AromaticRing && tj === FeatureType.PositiveCharge) ||
        (ti === FeatureType.PositiveCharge && tj === FeatureType.AromaticRing)
    );
}

const tmpPointA = Vec3();
const tmpPointB = Vec3();

function areFeaturesWithinDistanceSq(infoA: Features.Info, infoB: Features.Info, distanceSq: number): boolean {
    const { feature: featureA, offsets: offsetsA, members: membersA } = infoA;
    const { feature: featureB, offsets: offsetsB, members: membersB } = infoB;
    for (let i = offsetsA[featureA], il = offsetsA[featureA + 1]; i < il; ++i) {
        const elementA = membersA[i];
        infoA.unit.conformation.position(infoA.unit.elements[elementA], tmpPointA);
        for (let j = offsetsB[featureB], jl = offsetsB[featureB + 1]; j < jl; ++j) {
            const elementB = membersB[j];
            infoB.unit.conformation.position(infoB.unit.elements[elementB], tmpPointB);
            if (Vec3.squaredDistance(tmpPointA, tmpPointB) < distanceSq) return true;
        }
    }
    return false;
}

const tmpVecA = Vec3();
const tmpVecB = Vec3();
const tmpVecC = Vec3();
const tmpVecD = Vec3();

function getNormal(out: Vec3, info: Features.Info) {
    const { unit, feature, offsets, members } = info;
    const { elements } = unit;

    const i = offsets[feature];
    info.unit.conformation.position(elements[members[i]], tmpVecA);
    info.unit.conformation.position(elements[members[i + 1]], tmpVecB);
    info.unit.conformation.position(elements[members[i + 2]], tmpVecC);

    return Vec3.triangleNormal(out, tmpVecA, tmpVecB, tmpVecC);
}

const getOffset = function (infoA: Features.Info, infoB: Features.Info, normal: Vec3) {
    Features.position(tmpVecA, infoA);
    Features.position(tmpVecB, infoB);

    Vec3.sub(tmpVecC, tmpVecA, tmpVecB);

    Vec3.projectOnPlane(tmpVecD, tmpVecC, normal);
    Vec3.add(tmpVecD, tmpVecD, tmpVecB);
    return Vec3.distance(tmpVecD, tmpVecB);
};

function getIonicOptions(props: IonicProps) {
    return {
        distanceMaxSq: props.distanceMax * props.distanceMax,
    };
}
type IonicOptions = ReturnType<typeof getIonicOptions>

function getPiStackingOptions(props: PiStackingProps) {
    return {
        offsetMax: props.offsetMax,
        angleDevMax: degToRad(props.angleDevMax),
    };
}
type PiStackingOptions = ReturnType<typeof getPiStackingOptions>

function getCationPiOptions(props: CationPiProps) {
    return {
        offsetMax: props.offsetMax
    };
}
type CationPiOptions = ReturnType<typeof getCationPiOptions>

const deg180InRad = degToRad(180);
const deg90InRad = degToRad(90);

const tmpNormalA = Vec3();
const tmpNormalB = Vec3();

function testIonic(structure: Structure, infoA: Features.Info, infoB: Features.Info, distanceSq: number, opts: IonicOptions): InteractionType | undefined {
    const typeA = infoA.types[infoA.feature];
    const typeB = infoB.types[infoB.feature];

    if (isIonic(typeA, typeB)) {
        if (areFeaturesWithinDistanceSq(infoA, infoB, opts.distanceMaxSq)) {
            return InteractionType.Ionic;
        }
    }
}

function testPiStacking(structure: Structure, infoA: Features.Info, infoB: Features.Info, distanceSq: number, opts: PiStackingOptions): InteractionType | undefined {
    const typeA = infoA.types[infoA.feature];
    const typeB = infoB.types[infoB.feature];

    if (isPiStacking(typeA, typeB)) {
        getNormal(tmpNormalA, infoA);
        getNormal(tmpNormalB, infoB);

        const angle = Vec3.angle(tmpNormalA, tmpNormalB);
        const offset = Math.min(getOffset(infoA, infoB, tmpNormalB), getOffset(infoB, infoA, tmpNormalA));
        if (offset <= opts.offsetMax) {
            if (angle <= opts.angleDevMax || angle >= deg180InRad - opts.angleDevMax) {
                return InteractionType.PiStacking;  // parallel
            } else if (angle <= opts.angleDevMax + deg90InRad && angle >= deg90InRad - opts.angleDevMax) {
                return InteractionType.PiStacking;  // t-shaped
            }
        }
    }
}

function testCationPi(structure: Structure, infoA: Features.Info, infoB: Features.Info, distanceSq: number, opts: CationPiOptions): InteractionType | undefined {
    const typeA = infoA.types[infoA.feature];
    const typeB = infoB.types[infoB.feature];

    if (isCationPi(typeA, typeB)) {
        const [infoR, infoC] = typeA === FeatureType.AromaticRing ? [infoA, infoB] : [infoB, infoA];

        getNormal(tmpNormalA, infoR);
        const offset = getOffset(infoC, infoR, tmpNormalA);
        if (offset <= opts.offsetMax) {
            return InteractionType.CationPi;
        }
    }
}

//

export const NegativChargeProvider = Features.Provider([FeatureType.NegativeCharge], addUnitNegativeCharges);
export const PositiveChargeProvider = Features.Provider([FeatureType.PositiveCharge], addUnitPositiveCharges);
export const AromaticRingProvider = Features.Provider([FeatureType.AromaticRing], addUnitAromaticRings);

export const IonicProvider: ContactProvider<IonicParams> = {
    name: 'ionic',
    params: IonicParams,
    createTester: (props: IonicProps) => {
        const opts = getIonicOptions(props);
        return {
            maxDistance: props.distanceMax,
            requiredFeatures: new Set([FeatureType.NegativeCharge, FeatureType.PositiveCharge]),
            getType: (structure, infoA, infoB, distanceSq) => testIonic(structure, infoA, infoB, distanceSq, opts)
        };
    }
};

export const PiStackingProvider: ContactProvider<PiStackingParams> = {
    name: 'pi-stacking',
    params: PiStackingParams,
    createTester: (props: PiStackingProps) => {
        const opts = getPiStackingOptions(props);
        return {
            maxDistance: props.distanceMax,
            requiredFeatures: new Set([FeatureType.AromaticRing]),
            getType: (structure, infoA, infoB, distanceSq) => testPiStacking(structure, infoA, infoB, distanceSq, opts)
        };
    }
};

export const CationPiProvider: ContactProvider<CationPiParams> = {
    name: 'cation-pi',
    params: CationPiParams,
    createTester: (props: CationPiProps) => {
        const opts = getCationPiOptions(props);
        return {
            maxDistance: props.distanceMax,
            requiredFeatures: new Set([FeatureType.AromaticRing, FeatureType.PositiveCharge]),
            getType: (structure, infoA, infoB, distanceSq) => testCationPi(structure, infoA, infoB, distanceSq, opts)
        };
    }
};