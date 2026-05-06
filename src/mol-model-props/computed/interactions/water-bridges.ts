/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { Structure, Unit, StructureElement } from '../../../mol-model/structure';
import { IntMap } from '../../../mol-data/int';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { DataLocation } from '../../../mol-model/location';
import { DataLoci } from '../../../mol-model/loci';
import { MoleculeType } from '../../../mol-model/structure/model/types';
import { StructureLookup3DResultContext } from '../../../mol-model/structure/structure/util/lookup3d';
import { Sphere3D } from '../../../mol-math/geometry';
import { CentroidHelper } from '../../../mol-math/geometry/centroid-helper';
import { Features } from './features';
import { FeatureType } from './common';
import { GeometryOptions, checkGeometry } from './hydrogen-bonds';
import { degToRad } from '../../../mol-math/misc';
import { elementLabel } from '../../../mol-theme/label';

export type { WaterBridgeContact, WaterBridgeContacts };

interface WaterBridgeContact {
    /** non-water donor unit id */
    readonly unitA: number
    /** donor feature index in unitA */
    readonly indexA: Features.FeatureIndex
    /** non-water acceptor unit id */
    readonly unitB: number
    /** acceptor feature index in unitB */
    readonly indexB: Features.FeatureIndex
    /** bridging water unit id */
    readonly unitW: number
    /** water oxygen as HydrogenAcceptor (leg: donor → water) */
    readonly indexWA: Features.FeatureIndex
    /** water oxygen as HydrogenDonor (leg: water → acceptor) */
    readonly indexWD: Features.FeatureIndex
}

type WaterBridgeContacts = ReadonlyArray<WaterBridgeContact>;

export const WaterBridgesParams = {
    backbone: PD.Boolean(true, { description: 'Include backbone hydrogen bonds' }),
    ignoreHydrogens: PD.Boolean(true, { description: 'Ignore explicit hydrogens in geometric constraints' }),
    legDistMin: PD.Numeric(2.5, { min: 1, max: 4, step: 0.1 }, { description: 'Minimum leg distance (Å)' }),
    legDistMax: PD.Numeric(4.1, { min: 1, max: 6, step: 0.1 }, { description: 'Maximum leg distance (Å)' }),
    donAngleDevMax: PD.Numeric(80, { min: 0, max: 180, step: 1 }, { description: 'Max deviation from ideal donor angle' }),
    accAngleDevMax: PD.Numeric(50, { min: 0, max: 180, step: 1 }, { description: 'Max deviation from ideal acceptor angle' }),
    donOutOfPlaneAngleMax: PD.Numeric(45, { min: 0, max: 180, step: 1 }),
    accOutOfPlaneAngleMax: PD.Numeric(90, { min: 0, max: 180, step: 1 }),
    omegaMin: PD.Numeric(71, { min: 0, max: 180, step: 1 }, { description: 'Minimum A–W–B angle (°)' }),
    omegaMax: PD.Numeric(140, { min: 0, max: 180, step: 1 }, { description: 'Maximum A–W–B angle (°)' }),
};
export type WaterBridgesParams = typeof WaterBridgesParams;
export type WaterBridgesProps = PD.Values<WaterBridgesParams>;

function isWater(unit: Unit.Atomic, index: StructureElement.UnitIndex): boolean {
    return unit.model.atomicHierarchy.derived.residue.moleculeType[
        unit.residueIndex[unit.elements[index]]
    ] === MoleculeType.Water;
}

const _lookupCtx = StructureLookup3DResultContext();
const _vA = Vec3();
const _vB = Vec3();

function checkOmega(posA: Vec3, posW: Vec3, posB: Vec3, omegaMinRad: number, omegaMaxRad: number): boolean {
    Vec3.sub(_vA, posA, posW);
    Vec3.sub(_vB, posB, posW);
    Vec3.normalize(_vA, _vA);
    Vec3.normalize(_vB, _vB);
    const omega = Math.acos(Math.max(-1, Math.min(1, Vec3.dot(_vA, _vB))));
    return omega >= omegaMinRad && omega <= omegaMaxRad;
}

export function findWaterBridgeContacts(
    structure: Structure,
    unitsFeatures: IntMap<Features>,
    props: WaterBridgesProps
): WaterBridgeContacts {
    const legOpts: GeometryOptions = {
        ignoreHydrogens: props.ignoreHydrogens,
        includeBackbone: props.backbone,
        maxAccAngleDev: degToRad(props.accAngleDevMax),
        maxDonAngleDev: degToRad(props.donAngleDevMax),
        maxAccOutOfPlaneAngle: degToRad(props.accOutOfPlaneAngleMax),
        maxDonOutOfPlaneAngle: degToRad(props.donOutOfPlaneAngleMax),
    };
    const legDistMinSq = props.legDistMin * props.legDistMin;
    const legDistMaxSq = props.legDistMax * props.legDistMax;
    const omegaMinRad = degToRad(props.omegaMin);
    const omegaMaxRad = degToRad(props.omegaMax);

    type Candidate = { unit: Unit.Atomic; featureIdx: Features.FeatureIndex; pos: Vec3; distSq: number };

    // Best bridge per unique (donor, acceptor) feature pair across all water molecules.
    const best = new Map<string, { contact: WaterBridgeContact; combinedDistSq: number }>();
    const wPos = Vec3();

    for (const unitW of structure.units) {
        if (!Unit.isAtomic(unitW)) continue;
        const featW = unitsFeatures.get(unitW.id);
        if (!featW || featW.count === 0) continue;

        // Map each water-oxygen local index to its acceptor and donor feature indices.
        const waterMap = new Map<StructureElement.UnitIndex, {
            acc: Features.FeatureIndex | undefined,
            don: Features.FeatureIndex | undefined
        }>();

        for (let fi = 0 as Features.FeatureIndex; fi < featW.count; fi++) {
            const mi = featW.members[featW.offsets[fi]] as StructureElement.UnitIndex;
            if (!isWater(unitW, mi)) continue;
            const t = featW.types[fi];
            if (t !== FeatureType.HydrogenAcceptor && t !== FeatureType.HydrogenDonor) continue;
            let e = waterMap.get(mi);
            if (!e) waterMap.set(mi, (e = { acc: undefined, don: undefined }));
            if (t === FeatureType.HydrogenAcceptor) e.acc = fi;
            else e.don = fi;
        }

        for (const [waterAtomIdx, { acc: accFW, don: donFW }] of waterMap) {
            if (accFW === undefined || donFW === undefined) continue;

            unitW.conformation.position(unitW.elements[waterAtomIdx], wPos);

            const infoWAcc = Features.Info(structure, unitW, featW);
            infoWAcc.feature = accFW;
            const infoWDon = Features.Info(structure, unitW, featW);
            infoWDon.feature = donFW;

            const { count, indices, units: hitUnits, squaredDistances } =
                structure.lookup3d.find(wPos[0], wPos[1], wPos[2], props.legDistMax, _lookupCtx);

            const donors: Candidate[] = [];
            const acceptors: Candidate[] = [];

            for (let r = 0; r < count; r++) {
                const hitUnit = hitUnits[r];
                if (!Unit.isAtomic(hitUnit)) continue;
                if (hitUnit === unitW) continue;

                const hitLocalIdx = indices[r] as StructureElement.UnitIndex;
                if (isWater(hitUnit as Unit.Atomic, hitLocalIdx)) continue;

                const dSq = squaredDistances[r];
                if (dSq < legDistMinSq || dSq > legDistMaxSq) continue;

                const hitFeat = unitsFeatures.get(hitUnit.id);
                if (!hitFeat || hitFeat.count === 0) continue;

                const { indices: fIdxs, offsets: fOff } = hitFeat.elementsIndex;
                for (let k = fOff[hitLocalIdx], kl = fOff[hitLocalIdx + 1]; k < kl; k++) {
                    const fi = fIdxs[k] as Features.FeatureIndex;
                    const fType = hitFeat.types[fi];
                    if (fType !== FeatureType.HydrogenDonor && fType !== FeatureType.HydrogenAcceptor) continue;

                    const atomicUnit = hitUnit as Unit.Atomic;
                    const memberIdx = hitFeat.members[hitFeat.offsets[fi]] as StructureElement.UnitIndex;
                    const candidatePos = Vec3();
                    atomicUnit.conformation.position(atomicUnit.elements[memberIdx], candidatePos);

                    if (fType === FeatureType.HydrogenDonor) {
                        const infoDon = Features.Info(structure, atomicUnit, hitFeat);
                        infoDon.feature = fi;
                        if (checkGeometry(structure, infoDon, infoWAcc, legOpts)) {
                            donors.push({ unit: atomicUnit, featureIdx: fi, pos: candidatePos, distSq: dSq });
                        }
                    } else {
                        const infoAcc = Features.Info(structure, atomicUnit, hitFeat);
                        infoAcc.feature = fi;
                        if (checkGeometry(structure, infoWDon, infoAcc, legOpts)) {
                            acceptors.push({ unit: atomicUnit, featureIdx: fi, pos: candidatePos, distSq: dSq });
                        }
                    }
                }
            }

            for (const don of donors) {
                for (const acc of acceptors) {
                    if (don.unit === acc.unit && don.featureIdx === acc.featureIdx) continue;
                    if (!checkOmega(don.pos, wPos, acc.pos, omegaMinRad, omegaMaxRad)) continue;

                    const combinedDistSq = don.distSq + acc.distSq;
                    const pairKey = `${don.unit.id}|${don.featureIdx}|${acc.unit.id}|${acc.featureIdx}`;

                    const existing = best.get(pairKey);
                    if (!existing || combinedDistSq < existing.combinedDistSq) {
                        best.set(pairKey, {
                            contact: {
                                unitA: don.unit.id, indexA: don.featureIdx,
                                unitB: acc.unit.id, indexB: acc.featureIdx,
                                unitW: unitW.id, indexWA: accFW, indexWD: donFW,
                            },
                            combinedDistSq,
                        });
                    }
                }
            }
        }
    }

    return Array.from(best.values()).map(e => e.contact);
}

// ---------------------------------------------------------------------------
// Location / Loci for use by the renderer and color theme.
// ---------------------------------------------------------------------------

export { WaterBridges };

namespace WaterBridges {
    export interface Data {
        readonly structure: Structure
        readonly waterBridges: WaterBridgeContacts
        readonly unitsFeatures: IntMap<Features>
    }

    export interface Element { bridgeIndex: number }

    export interface Location extends DataLocation<Data, Element> {}

    export function Location(data: Data, bridgeIndex = 0): Location {
        return DataLocation('water-bridges', data, { bridgeIndex });
    }

    export function isLocation(x: any): x is Location {
        return !!x && x.kind === 'data-location' && x.tag === 'water-bridges';
    }

    export interface Loci extends DataLoci<Data, Element> {}

    export function Loci(data: Data, elements: ReadonlyArray<Element>): Loci {
        return DataLoci('water-bridges', data, elements,
            bs => getBoundingSphere(data, elements, bs),
            () => getLabel(data, elements));
    }

    function getLabel(data: Data, elements: ReadonlyArray<Element>): string {
        const e = elements[0];
        if (e === undefined) return '';
        const { structure, waterBridges, unitsFeatures } = data;
        const wb = waterBridges[e.bridgeIndex];

        const uA = structure.unitMap.get(wb.unitA) as Unit.Atomic;
        const fA = unitsFeatures.get(wb.unitA);
        const uW = structure.unitMap.get(wb.unitW) as Unit.Atomic;
        const fW = unitsFeatures.get(wb.unitW);
        const uB = structure.unitMap.get(wb.unitB) as Unit.Atomic;
        const fB = unitsFeatures.get(wb.unitB);

        const locA = StructureElement.Location.create(structure, uA, uA.elements[fA.members[fA.offsets[wb.indexA]]]);
        const locW = StructureElement.Location.create(structure, uW, uW.elements[fW.members[fW.offsets[wb.indexWA]]]);
        const locB = StructureElement.Location.create(structure, uB, uB.elements[fB.members[fB.offsets[wb.indexB]]]);

        return [
            'Water Bridge',
            `${elementLabel(locA, { granularity: 'element' })} → ${elementLabel(locW, { granularity: 'residue' })} → ${elementLabel(locB, { granularity: 'element' })}`,
        ].join('</br>');
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'data-loci' && x.tag === 'water-bridges';
    }

    function getBoundingSphere(data: Data, elements: ReadonlyArray<Element>, boundingSphere: Sphere3D) {
        return CentroidHelper.fromPairProvider(elements.length, (i, pA, pB) => {
            const wb = data.waterBridges[elements[i].bridgeIndex];
            const uA = data.structure.unitMap.get(wb.unitA) as Unit.Atomic;
            const fA = data.unitsFeatures.get(wb.unitA);
            uA.conformation.position(uA.elements[fA.members[fA.offsets[wb.indexA]]], pA);
            const uB = data.structure.unitMap.get(wb.unitB) as Unit.Atomic;
            const fB = data.unitsFeatures.get(wb.unitB);
            uB.conformation.position(uB.elements[fB.members[fB.offsets[wb.indexB]]], pB);
        }, boundingSphere);
    }
}
