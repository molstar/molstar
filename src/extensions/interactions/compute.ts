/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FeatureType } from '../../mol-model-props/computed/interactions/common';
import { computeInteractions as _compute, Interactions } from '../../mol-model-props/computed/interactions/interactions';
import { Structure, StructureElement } from '../../mol-model/structure';
import { RuntimeContext } from '../../mol-task';
import { AssetManager } from '../../mol-util/assets';
import { InteractionInfo, InteractionTypeToKind, StructureInteractions } from './model';

export interface ComputeInteractionsOptions {
    // includeInterStructure?: boolean,
    // kinds
}

export async function computeInteractions(
    ctx: RuntimeContext,
    loci: [ref: string, StructureElement.Loci][],
    options?: ComputeInteractionsOptions
): Promise<StructureInteractions> {
    const unitIdToRef = new Map<number, string>();

    const builder = Structure.Builder({ masterModel: loci[0][1].structure.models[0] });
    for (const [ref, l] of loci) {
        const s = StructureElement.Loci.toStructure(l);
        const unitIds: number[] = [];
        for (const unit of s.units) {
            const newUnit = builder.addUnit(unit.kind, unit.model, unit.conformation.operator, unit.elements, unit.traits);
            unitIds.push(newUnit.id);
            unitIdToRef.set(newUnit.id, ref);
        }
    }

    const structure = builder.getStructure();
    const interactions = await _compute({ runtime: ctx, assetManager: new AssetManager() }, structure, { });

    const { edges } = interactions.contacts;
    const result: StructureInteractions = { kind: 'structure-interactions', elements: [] };
    for (const e of edges) {
        if (e.unitA > e.unitB) continue;

        const [a, aType] = processFeature(structure, interactions, e.unitA, e.indexA);
        const [b] = processFeature(structure, interactions, e.unitB, e.indexB);

        const kind = InteractionTypeToKind[e.props.type] ?? 'unknown';
        const info: InteractionInfo = { kind };

        if (kind === 'hydrogen-bond' || kind === 'weak-hydrogen-bond') {
            const isADonor = aType === FeatureType.HydrogenDonor || aType === FeatureType.WeakHydrogenDonor;

            result.elements.push({
                info,
                aStructureRef: isADonor ? unitIdToRef.get(e.unitA)! : unitIdToRef.get(e.unitB)!,
                bStructureRef: isADonor ? unitIdToRef.get(e.unitB)! : unitIdToRef.get(e.unitA)!,
                a: isADonor ? a : b,
                b: isADonor ? b : a,
            });
        } else {
            result.elements.push({
                info,
                aStructureRef: unitIdToRef.get(e.unitA)!,
                bStructureRef: unitIdToRef.get(e.unitB)!,
                a,
                b,
            });
        }
    }

    return result;
}

const _loc = StructureElement.Location.create();
function processFeature(structure: Structure, interactions: Interactions, unitId: number, featureIndex: number) {
    _loc.structure = structure;
    _loc.unit = structure.unitMap.get(unitId);
    const xs = interactions.unitsFeatures.get(unitId)!;

    let type: FeatureType = FeatureType.None;
    const builder = structure.subsetBuilder(false);
    builder.beginUnit(_loc.unit.id);
    for (let o = xs.offsets[featureIndex], uIEnd = xs.offsets[featureIndex + 1]; o < uIEnd; o++) {
        const unitIndex = xs.members[o];
        _loc.element = _loc.unit.elements[unitIndex];
        builder.addElement(_loc.element);
        type = xs.types[o];
    }
    builder.commitUnit();

    return [Structure.toStructureElementLoci(builder.getStructure()), type] as const;
}
