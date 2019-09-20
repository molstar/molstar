/**
 * Copyright (c) 2017-2019 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { LinkType } from '../../../model/types';
import Structure from '../../structure';
import Unit from '../../unit';
import { getElementIdx, getElementPairThreshold, getElementThreshold, isHydrogen, LinkComputationProps, MetalsSet, DefaultLinkComputationProps } from './common';
import { InterUnitBonds } from './data';
import { UniqueArray } from '../../../../../mol-data/generic';
import { SortedArray } from '../../../../../mol-data/int';
import { Vec3, Mat4 } from '../../../../../mol-math/linear-algebra';
import StructureElement from '../../element';
import { StructConn } from '../../../../../mol-model-formats/structure/mmcif/bonds';

const MAX_RADIUS = 4;

function addMapEntry<A, B>(map: Map<A, B[]>, a: A, b: B) {
    if (map.has(a)) map.get(a)!.push(b);
    else map.set(a, [b]);
}

interface PairState {
    mapAB: Map<number, InterUnitBonds.BondInfo[]>,
    mapBA: Map<number, InterUnitBonds.BondInfo[]>,
    bondedA: UniqueArray<StructureElement.UnitIndex, StructureElement.UnitIndex>,
    bondedB: UniqueArray<StructureElement.UnitIndex, StructureElement.UnitIndex>
}

function addLink(indexA: number, indexB: number, order: number, flag: LinkType.Flag, state: PairState) {
    addMapEntry(state.mapAB, indexA, { indexB, order, flag });
    addMapEntry(state.mapBA, indexB, { indexB: indexA, order, flag });
    UniqueArray.add(state.bondedA, indexA, indexA);
    UniqueArray.add(state.bondedB, indexB, indexB);
}

const _imageTransform = Mat4.zero();

function findPairLinks(unitA: Unit.Atomic, unitB: Unit.Atomic, props: LinkComputationProps, map: Map<number, InterUnitBonds.UnitPairBonds[]>) {
    const state: PairState = { mapAB: new Map(), mapBA: new Map(), bondedA: UniqueArray.create(), bondedB: UniqueArray.create() };
    let bondCount = 0;

    const { elements: atomsA } = unitA;
    const { x: xA, y: yA, z: zA } = unitA.model.atomicConformation;
    const { elements: atomsB } = unitB;
    const atomCount = unitA.elements.length;

    const { type_symbol: type_symbolA, label_alt_id: label_alt_idA } = unitA.model.atomicHierarchy.atoms;
    const { type_symbol: type_symbolB, label_alt_id: label_alt_idB } = unitB.model.atomicHierarchy.atoms;
    const { lookup3d } = unitB;
    const structConn = unitA.model === unitB.model && unitA.model.sourceData.kind === 'mmCIF' ? StructConn.get(unitA.model) : void 0;

    // the lookup queries need to happen in the "unitB space".
    // that means imageA = inverseOperB(operA(aI))
    const imageTransform = Mat4.mul(_imageTransform, unitB.conformation.operator.inverse, unitA.conformation.operator.matrix);
    const isNotIdentity = !Mat4.isIdentity(imageTransform);
    const imageA = Vec3.zero();

    const { center: bCenter, radius: bRadius } = lookup3d.boundary.sphere;
    const testDistanceSq = (bRadius + MAX_RADIUS) * (bRadius + MAX_RADIUS);

    for (let _aI = 0; _aI < atomCount; _aI++) {
        const aI = atomsA[_aI];
        Vec3.set(imageA, xA[aI], yA[aI], zA[aI]);
        if (isNotIdentity) Vec3.transformMat4(imageA, imageA, imageTransform);
        if (Vec3.squaredDistance(imageA, bCenter) > testDistanceSq) continue;

        const structConnEntries = props.forceCompute ? void 0 : structConn && structConn.getAtomEntries(aI);
        if (structConnEntries && structConnEntries.length) {
            for (const se of structConnEntries) {
                if (se.distance < MAX_RADIUS) continue;

                for (const p of se.partners) {
                    const _bI = SortedArray.indexOf(unitB.elements, p.atomIndex);
                    if (_bI < 0) continue;
                    addLink(_aI, _bI, se.order, se.flags, state);
                    bondCount++;
                }
            }
        }

        const { indices, count, squaredDistances } = lookup3d.find(imageA[0], imageA[1], imageA[2], MAX_RADIUS);
        if (count === 0) continue;

        const aeI = getElementIdx(type_symbolA.value(aI));
        const isHa = isHydrogen(aeI);
        const thresholdA = getElementThreshold(aeI);
        const altA = label_alt_idA.value(aI);
        const metalA = MetalsSet.has(aeI);

        for (let ni = 0; ni < count; ni++) {
            const _bI = indices[ni];
            const bI = atomsB[_bI];

            const altB = label_alt_idB.value(bI);
            // TODO: check if they have the same model?
            if (altA && altB && altA !== altB) continue;

            const beI = getElementIdx(type_symbolB.value(bI)!);
            const isMetal = metalA || MetalsSet.has(beI);

            const isHb = isHydrogen(beI);
            if (isHa && isHb) continue;

            const dist = Math.sqrt(squaredDistances[ni]);
            if (dist === 0) continue;

            // handle "struct conn" bonds.
            if (structConnEntries && structConnEntries.length) {
                let added = false;
                for (const se of structConnEntries) {
                    for (const p of se.partners) {
                        if (p.atomIndex === bI) {
                            addLink(_aI, _bI, se.order, se.flags, state);
                            bondCount++;
                            added = true;
                            break;
                        }
                    }
                    if (added) break;
                }
                if (added) continue;
            }

            if (isHa || isHb) {
                if (dist < props.maxCovalentHydrogenBondingLength) {
                    addLink(_aI, _bI, 1, LinkType.Flag.Covalent | LinkType.Flag.Computed, state); // TODO: check if correct
                    bondCount++;
                }
                continue;
            }

            const thresholdAB = getElementPairThreshold(aeI, beI);
            const pairingThreshold = thresholdAB > 0
                ? thresholdAB
                : beI < 0 ? thresholdA : Math.max(thresholdA, getElementThreshold(beI));

            if (dist <= pairingThreshold) {
                addLink(_aI, _bI, 1, (isMetal ? LinkType.Flag.MetallicCoordination : LinkType.Flag.Covalent) | LinkType.Flag.Computed, state);
                bondCount++;
            }
        }
    }

    if (bondCount) {
        addMapEntry(map, unitA.id, new InterUnitBonds.UnitPairBonds(unitA, unitB, bondCount, state.bondedA.array, state.mapAB));
        addMapEntry(map, unitB.id, new InterUnitBonds.UnitPairBonds(unitB, unitA, bondCount, state.bondedB.array, state.mapBA));
    }

    return bondCount;
}

export interface InterLinkComputationProps extends LinkComputationProps {
    validUnitPair: (unitA: Unit, unitB: Unit) => boolean
}

function findLinks(structure: Structure, props: InterLinkComputationProps) {
    const map = new Map<number, InterUnitBonds.UnitPairBonds[]>();
    if (!structure.units.some(u => Unit.isAtomic(u))) return new InterUnitBonds(map);

    const { validUnitPair } = props;
    const lookup = structure.lookup3d;
    const imageCenter = Vec3.zero();

    for (const unit of structure.units) {
        if (!Unit.isAtomic(unit)) continue;

        const bs = unit.lookup3d.boundary.sphere;
        Vec3.transformMat4(imageCenter, bs.center, unit.conformation.operator.matrix);
        const closeUnits = lookup.findUnitIndices(imageCenter[0], imageCenter[1], imageCenter[2], bs.radius + MAX_RADIUS);
        for (let i = 0; i < closeUnits.count; i++) {
            const other = structure.units[closeUnits.indices[i]];
            if (!Unit.isAtomic(other) || unit.id >= other.id || !validUnitPair(unit, other)) continue;

            if (other.elements.length >= unit.elements.length) findPairLinks(unit, other, props, map);
            else findPairLinks(other, unit, props, map);
        }
    }

    return new InterUnitBonds(map);
}

function ValidUnitPair(structure: Structure) {
    const { masterModel } = structure
    if (masterModel) {
        return (a: Unit, b: Unit) => a.model === b.model || a.model === masterModel || b.model === masterModel
    } else {
        return (a: Unit, b: Unit) => a.model === b.model
    }
}

function computeInterUnitBonds(structure: Structure, props?: Partial<InterLinkComputationProps>): InterUnitBonds {
    return findLinks(structure, {
        ...DefaultLinkComputationProps,
        validUnitPair: (props && props.validUnitPair) || ValidUnitPair(structure),
    });
}

export { computeInterUnitBonds };
