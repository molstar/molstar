/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureSelection } from '../../mol-model/structure';
import { StructureElementSchema } from '../../mol-model/structure/query/schema';
import { StructureQueryHelper } from '../../mol-plugin-state/helpers/structure-query';
import { InteractionSchema, StructureInteractionElement, StructureInteractions } from './model';

export function getCustomInteractionData(interactions: InteractionSchema[], structures: { [ref: string]: Structure }): StructureInteractions {
    const elements: StructureInteractionElement[] = [];

    for (const schema of interactions) {
        elements.push({
            sourceSchema: schema,
            info: { kind: schema.kind },
            aStructureRef: schema.aStructureRef,
            a: resolveLoci(structures[schema.aStructureRef!], schema.a),
            bStructureRef: schema.bStructureRef,
            b: resolveLoci(structures[schema.bStructureRef!], schema.b)
        });
    }

    return { elements };
}

function resolveLoci(structure: Structure, schema: StructureElementSchema) {
    const expr = StructureElementSchema.toExpression(schema);
    const selection = StructureQueryHelper.createAndRun(structure, expr).selection;
    return StructureSelection.toLociWithSourceUnits(selection);
}