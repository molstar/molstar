/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { computeInteractions as _compute, Interactions } from '../../mol-model-props/computed/interactions/interactions';
import { Structure, StructureElement } from '../../mol-model/structure';
import { RuntimeContext } from '../../mol-task';
import { AssetManager } from '../../mol-util/assets';
import { InteractionTypeToKind, StructureInteractions } from './model';

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
    const result: StructureInteractions = { elements: [] };
    for (const e of edges) {
        if (e.unitA > e.unitB) continue;

        result.elements.push({
            info: { kind: InteractionTypeToKind[e.props.type] },
            aStructureRef: unitIdToRef.get(e.unitA)!,
            bStructureRef: unitIdToRef.get(e.unitB)!,
            a: toLoci(structure, interactions, e.unitA, e.indexA),
            b: toLoci(structure, interactions, e.unitB, e.indexB),
        });
    }

    return result;
}

const _loc = StructureElement.Location.create();
function toLoci(structure: Structure, interactions: Interactions, unitId: number, featureIndex: number) {
    _loc.structure = structure;
    _loc.unit = structure.unitMap.get(unitId);
    const xs = interactions.unitsFeatures.get(unitId)!;

    const builder = structure.subsetBuilder(false);
    builder.beginUnit(_loc.unit.id);
    for (let o = xs.offsets[featureIndex], uIEnd = xs.offsets[featureIndex + 1]; o < uIEnd; o++) {
        const unitIndex = xs.members[o];
        _loc.element = _loc.unit.elements[unitIndex];
        builder.addElement(_loc.element);
    }
    builder.commitUnit();

    return Structure.toStructureElementLoci(builder.getStructure());
}
