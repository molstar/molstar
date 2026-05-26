/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { Structure, Unit, StructureElement } from '../../../mol-model/structure';
import { IntMap, OrderedSet } from '../../../mol-data/int';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { DataLocation } from '../../../mol-model/location';
import { DataLoci } from '../../../mol-model/loci';
import { MoleculeType, NucleicBackboneAtoms, ProteinBackboneAtoms } from '../../../mol-model/structure/model/types';
import { StructureLookup3DResultContext } from '../../../mol-model/structure/structure/util/lookup3d';
import { Sphere3D } from '../../../mol-math/geometry';
import { CentroidHelper } from '../../../mol-math/geometry/centroid-helper';
import { Features } from './features';
import { FeatureType, InteractionType, InteractionFlag } from './common';
import { GeometryOptions, checkGeometry } from './hydrogen-bonds';
import { degToRad } from '../../../mol-math/misc';
import { bundleLabel, LabelGranularity } from '../../../mol-theme/label';
import { cantorPairing } from '../../../mol-data/util/hash-functions';

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
    props: { type: InteractionType.WaterBridge, flag: InteractionFlag }
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

export const WaterBridgesProvider = {
    requiredFeatures: new Set([FeatureType.HydrogenDonor, FeatureType.HydrogenAcceptor]),
    params: WaterBridgesParams,
    find: findWaterBridgeContacts,
};

function isWater(unit: Unit.Atomic, index: StructureElement.UnitIndex): boolean {
    return unit.model.atomicHierarchy.derived.residue.moleculeType[
        unit.residueIndex[unit.elements[index]]
        ] === MoleculeType.Water;
}

function isBackboneAtom(unit: Unit.Atomic, index: StructureElement.UnitIndex): boolean {
    const element = unit.elements[index];
    const moleculeType = unit.model.atomicHierarchy.derived.residue.moleculeType[unit.residueIndex[element]];
    if (moleculeType !== MoleculeType.Protein && moleculeType !== MoleculeType.RNA && moleculeType !== MoleculeType.DNA) {
        return false;
    }

    const atomId = unit.model.atomicHierarchy.atoms.label_atom_id.value(element);
    if (moleculeType === MoleculeType.Protein) {
        return ProteinBackboneAtoms.has(atomId);
    }

    return NucleicBackboneAtoms.has(atomId);
}

const _lookupCtx = StructureLookup3DResultContext();

type Candidate = {
    unit: Unit.Atomic
    featureIdx: Features.FeatureIndex
    memberIdx: StructureElement.UnitIndex
    x: number
    y: number
    z: number
    distSq: number
};

type FeatureKey = number;

function featureKey(unitId: number, featureIndex: Features.FeatureIndex): FeatureKey {
    return cantorPairing(unitId, featureIndex);
}

type BestBridge = { contact: WaterBridgeContact; combinedDistSq: number };
type BestBridgeMap = Map<FeatureKey, Map<FeatureKey, BestBridge>>;

function getBestBridge(best: BestBridgeMap, donorKey: FeatureKey, acceptorKey: FeatureKey): BestBridge | undefined {
    return best.get(donorKey)?.get(acceptorKey);
}

function setBestBridge(best: BestBridgeMap, donorKey: FeatureKey, acceptorKey: FeatureKey, value: BestBridge) {
    let acceptors = best.get(donorKey);
    if (acceptors === undefined) {
        acceptors = new Map();
        best.set(donorKey, acceptors);
    }
    acceptors.set(acceptorKey, value);
}

function bestBridgeValues(best: BestBridgeMap): BestBridge[] {
    const values: BestBridge[] = [];
    for (const acceptors of best.values()) {
        for (const value of acceptors.values()) values.push(value);
    }
    return values;
}

function checkOmega(don: Candidate, posW: Vec3, acc: Candidate, cosOmegaMin: number, cosOmegaMax: number): boolean {
    const ax = don.x - posW[0];
    const ay = don.y - posW[1];
    const az = don.z - posW[2];

    const bx = acc.x - posW[0];
    const by = acc.y - posW[1];
    const bz = acc.z - posW[2];

    const aLenSq = ax * ax + ay * ay + az * az;
    const bLenSq = bx * bx + by * by + bz * bz;

    if (aLenSq === 0 || bLenSq === 0) return false;

    const cosOmega = (ax * bx + ay * by + az * bz) / Math.sqrt(aLenSq * bLenSq);

    // cos decreases monotonically on [0, pi], so:
    // omega >= omegaMin && omega <= omegaMax
    // is equivalent to:
    // cos(omega) <= cos(omegaMin) && cos(omega) >= cos(omegaMax)
    return cosOmega <= cosOmegaMin && cosOmega >= cosOmegaMax;
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

    if (omegaMinRad > omegaMaxRad) return [];

    const cosOmegaMin = Math.cos(omegaMinRad);
    const cosOmegaMax = Math.cos(omegaMaxRad);

    // Best bridge per unique donor/acceptor feature pair across all water molecules.
    const best: BestBridgeMap = new Map();

    const wPos = Vec3();
    const candidatePos = Vec3();

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

        if (waterMap.size === 0) continue;

        const infoWAcc = Features.Info(structure, unitW, featW);
        const infoWDon = Features.Info(structure, unitW, featW);

        for (const [waterAtomIdx, { acc: accFW, don: donFW }] of waterMap) {
            if (accFW === undefined || donFW === undefined) continue;

            unitW.conformation.position(unitW.elements[waterAtomIdx], wPos);

            infoWAcc.feature = accFW;
            infoWDon.feature = donFW;

            const { count, indices, units: hitUnits } =
                structure.lookup3d.find(wPos[0], wPos[1], wPos[2], props.legDistMax, _lookupCtx);

            const donors: Candidate[] = [];
            const acceptors: Candidate[] = [];

            const donorKeys = new Set<FeatureKey>();
            const acceptorKeys = new Set<FeatureKey>();

            for (let r = 0; r < count; r++) {
                const hitUnit = hitUnits[r];
                if (!Unit.isAtomic(hitUnit)) continue;

                const atomicUnit = hitUnit as Unit.Atomic;
                const hitLocalIdx = indices[r] as StructureElement.UnitIndex;

                // Only skip the water atom itself. Other atoms in the same unit can still be valid.
                if (atomicUnit === unitW && hitLocalIdx === waterAtomIdx) continue;
                if (isWater(atomicUnit, hitLocalIdx)) continue;

                const hitFeat = unitsFeatures.get(atomicUnit.id);
                if (!hitFeat || hitFeat.count === 0) continue;

                const infoHit = Features.Info(structure, atomicUnit, hitFeat);

                const { indices: fIdxs, offsets: fOff } = hitFeat.elementsIndex;
                for (let k = fOff[hitLocalIdx], kl = fOff[hitLocalIdx + 1]; k < kl; k++) {
                    const fi = fIdxs[k] as Features.FeatureIndex;
                    const fType = hitFeat.types[fi];

                    if (fType !== FeatureType.HydrogenDonor && fType !== FeatureType.HydrogenAcceptor) continue;

                    const memberIdx = hitFeat.members[hitFeat.offsets[fi]] as StructureElement.UnitIndex;

                    if (!props.backbone && isBackboneAtom(atomicUnit, memberIdx)) continue;

                    atomicUnit.conformation.position(atomicUnit.elements[memberIdx], candidatePos);

                    const distSq = Vec3.squaredDistance(candidatePos, wPos);
                    if (distSq < legDistMinSq || distSq > legDistMaxSq) continue;

                    infoHit.feature = fi;

                    if (fType === FeatureType.HydrogenDonor) {
                        const key = featureKey(atomicUnit.id, fi);
                        if (donorKeys.has(key)) continue;

                        if (checkGeometry(structure, infoHit, infoWAcc, legOpts)) {
                            donorKeys.add(key);
                            donors.push({
                                unit: atomicUnit,
                                featureIdx: fi,
                                memberIdx,
                                x: candidatePos[0],
                                y: candidatePos[1],
                                z: candidatePos[2],
                                distSq,
                            });
                        }
                    } else {
                        const key = featureKey(atomicUnit.id, fi);
                        if (acceptorKeys.has(key)) continue;

                        if (checkGeometry(structure, infoWDon, infoHit, legOpts)) {
                            acceptorKeys.add(key);
                            acceptors.push({
                                unit: atomicUnit,
                                featureIdx: fi,
                                memberIdx,
                                x: candidatePos[0],
                                y: candidatePos[1],
                                z: candidatePos[2],
                                distSq,
                            });
                        }
                    }
                }
            }

            for (const don of donors) {
                for (const acc of acceptors) {
                    // Reject bridges where donor and acceptor are the same physical atom
                    // represented by different feature indices.
                    if (don.unit === acc.unit && don.memberIdx === acc.memberIdx) continue;

                    if (!checkOmega(don, wPos, acc, cosOmegaMin, cosOmegaMax)) continue;

                    const combinedDistSq = don.distSq + acc.distSq;
                    const donorKey = featureKey(don.unit.id, don.featureIdx);
                    const acceptorKey = featureKey(acc.unit.id, acc.featureIdx);

                    const existing = getBestBridge(best, donorKey, acceptorKey);
                    if (!existing || combinedDistSq < existing.combinedDistSq) {
                        setBestBridge(best, donorKey, acceptorKey, {
                            contact: {
                                unitA: don.unit.id,
                                indexA: don.featureIdx,
                                unitB: acc.unit.id,
                                indexB: acc.featureIdx,
                                unitW: unitW.id,
                                indexWA: accFW,
                                indexWD: donFW,
                                props: { type: InteractionType.WaterBridge, flag: InteractionFlag.None },
                            },
                            combinedDistSq,
                        });
                    }
                }
            }
        }
    }

    return bestBridgeValues(best).map(e => e.contact);
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

        const options = { granularity: 'element' as LabelGranularity };
        if (fA.offsets[wb.indexA + 1] - fA.offsets[wb.indexA] > 1 ||
                fB.offsets[wb.indexB + 1] - fB.offsets[wb.indexB] > 1) {
            options.granularity = 'residue';
        }

        return [
            'Water Bridge',
            bundleLabel({ loci: [
                StructureElement.Loci(structure, [{ unit: uA, indices: OrderedSet.ofSingleton(fA.members[fA.offsets[wb.indexA]] as StructureElement.UnitIndex) }]),
                StructureElement.Loci(structure, [{ unit: uW, indices: OrderedSet.ofSingleton(fW.members[fW.offsets[wb.indexWA]] as StructureElement.UnitIndex) }]),
                StructureElement.Loci(structure, [{ unit: uB, indices: OrderedSet.ofSingleton(fB.members[fB.offsets[wb.indexB]] as StructureElement.UnitIndex) }]),
            ] }, options),
        ].join('</br>');
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'data-loci' && x.tag === 'water-bridges';
    }

    function getBoundingSphere(data: Data, elements: ReadonlyArray<Element>, boundingSphere: Sphere3D) {
        return CentroidHelper.fromPairProvider(elements.length * 2, (i, pA, pB) => {
            const wb = data.waterBridges[elements[i >> 1].bridgeIndex];

            const uA = data.structure.unitMap.get(wb.unitA) as Unit.Atomic;
            const fA = data.unitsFeatures.get(wb.unitA);

            const uW = data.structure.unitMap.get(wb.unitW) as Unit.Atomic;
            const fW = data.unitsFeatures.get(wb.unitW);

            const uB = data.structure.unitMap.get(wb.unitB) as Unit.Atomic;
            const fB = data.unitsFeatures.get(wb.unitB);

            const aIdx = fA.members[fA.offsets[wb.indexA]] as StructureElement.UnitIndex;
            const wIdx = fW.members[fW.offsets[wb.indexWA]] as StructureElement.UnitIndex;
            const bIdx = fB.members[fB.offsets[wb.indexB]] as StructureElement.UnitIndex;

            if ((i & 1) === 0) {
                uA.conformation.position(uA.elements[aIdx], pA);
                uW.conformation.position(uW.elements[wIdx], pB);
            } else {
                uW.conformation.position(uW.elements[wIdx], pA);
                uB.conformation.position(uB.elements[bIdx], pB);
            }
        }, boundingSphere);
    }
}