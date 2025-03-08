/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { computeInteractions as _compute } from '../../mol-model-props/computed/interactions/interactions';
import { Structure, StructureElement, StructureProperties } from '../../mol-model/structure';
import { StructureElementSchema } from '../../mol-model/structure/query/schema';
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

        result.elements.push({
            kind: InteractionTypeToKind[e.props.type],
            aStructureRef: unitIdToRef.get(e.unitA),
            a: toSchema(structure, e.unitA, e.indexA),
            bStructureRef: unitIdToRef.get(e.unitB),
            b: toSchema(structure, e.unitB, e.indexB),
        });
    }

    return result;
}

const _loc = StructureElement.Location.create();
function toSchema(structure: Structure, unitIndex: number, elementIndex: number): StructureElementSchema {
    _loc.structure = structure;
    _loc.unit = structure.units[unitIndex];
    _loc.element = _loc.unit.elements[elementIndex];
    return {
        label_entity_id: StructureProperties.chain.label_entity_id(_loc),
        label_asym_id: StructureProperties.chain.label_asym_id(_loc),
        label_seq_id: StructureProperties.residue.label_seq_id(_loc),
        // TODO: use atom name if unique within the residue
        atom_index: StructureProperties.atom.sourceIndex(_loc),
    };
}
