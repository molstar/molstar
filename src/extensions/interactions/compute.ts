/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { computeInteractions as _compute, Interactions } from '../../mol-model-props/computed/interactions/interactions';
import { Structure, StructureElement } from '../../mol-model/structure';
import { structureElementLocationToSchemaItem, StructureElementSchema, StructureElementSchemaItem } from '../../mol-model/structure/query/schema';
import { RuntimeContext } from '../../mol-task';
import { AssetManager } from '../../mol-util/assets';
import { InteractionTypeToKind, StructureInteractions } from './model';

export interface ComputeInteractionsOptions {
    includeInterStructure?: boolean,
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

        const [lociA, a] = toSchema(structure, interactions, e.unitA, e.indexA);
        const [lociB, b] = toSchema(structure, interactions, e.unitB, e.indexB);

        result.elements.push({
            schema: {
                kind: InteractionTypeToKind[e.props.type],
                aStructureRef: unitIdToRef.get(e.unitA),
                a,
                bStructureRef: unitIdToRef.get(e.unitB),
                b,
            },
            a: lociA,
            b: lociB,
        });
    }

    return result;
}

const _loc = StructureElement.Location.create();
function toSchema(structure: Structure, interactions: Interactions, unitId: number, featureIndex: number): [StructureElement.Loci, StructureElementSchema] {
    _loc.structure = structure;
    _loc.unit = structure.unitMap.get(unitId);
    const fx = interactions.unitsFeatures.get(unitId)!;

    const builder = structure.subsetBuilder(false);
    builder.beginUnit(_loc.unit.id);
    const items: StructureElementSchemaItem[] = [];
    for (let o = fx.offsets[featureIndex], uIEnd = fx.offsets[featureIndex + 1]; o < uIEnd; o++) {
        const unitIndex = fx.members[o];
        _loc.element = _loc.unit.elements[unitIndex];
        items.push(structureElementLocationToSchemaItem(_loc));
        builder.addElement(_loc.element);
    }
    builder.commitUnit();

    const schema = items.length === 1 ? items[0] : items;
    const loci = Structure.toStructureElementLoci(builder.getStructure());

    return [loci, schema];
}
