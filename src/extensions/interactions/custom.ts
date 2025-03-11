/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureSelection } from '../../mol-model/structure';
import { StructureElementSchema } from '../../mol-model/structure/query/schema';
import { StructureQueryHelper } from '../../mol-plugin-state/helpers/structure-query';
import { InteractionInfo, InteractionElementSchema, StructureInteractionElement, StructureInteractions } from './model';

export function getCustomInteractionData(interactions: InteractionElementSchema[], structures: { [ref: string]: Structure }): StructureInteractions {
    const elements: StructureInteractionElement[] = [];

    for (const schema of interactions) {
        let info: InteractionInfo;
        if (schema.kind === 'hydrogen-bond' || schema.kind === 'weak-hydrogen-bond') {
            info = {
                kind: schema.kind,
                hydrogenStructureRef: schema.hydrogenStructureRef,
                hydrogen: schema.hydrogen ? resolveLoci(structures[schema.hydrogenStructureRef!], schema.hydrogen) : undefined,
            };
        } else if (schema.kind === 'covalent') {
            info = { kind: schema.kind, degree: schema.degree, aromatic: schema.aromatic };
        } else {
            info = { kind: schema.kind };
        }
        elements.push({
            sourceSchema: schema,
            info,
            aStructureRef: schema.aStructureRef,
            a: resolveLoci(structures[schema.aStructureRef!], schema.a),
            bStructureRef: schema.bStructureRef,
            b: resolveLoci(structures[schema.bStructureRef!], schema.b),
        });
    }

    return { kind: 'structure-interactions', elements };
}

function resolveLoci(structure: Structure, schema: StructureElementSchema) {
    const expr = StructureElementSchema.toExpression(schema);
    const selection = StructureQueryHelper.createAndRun(structure, expr).selection;
    return StructureSelection.toLociWithSourceUnits(selection);
}