/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureElement } from '../../mol-model/structure';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { computeInteractions, ComputeInteractionsOptions } from './compute';

const Factory = StateTransformer.builderFactory('interactions-extension');

export class InteractionData extends PluginStateObject.Create<{ interactions: any }>({ name: 'Interactions', typeClass: 'Data' }) { }

export interface ComputeInteractionSource {
    structureRef: string,
    bundle?: StructureElement.Bundle,
    // schema?: StructureElementSchema,
}

export const ComputeInteractions = Factory({
    name: 'compute-interactions',
    display: 'Compute Interactions',
    from: PluginStateObject.Root,
    to: InteractionData,
    params: {
        sources: PD.Value<ComputeInteractionSource[]>([], { isHidden: true }),
        options: PD.Value<ComputeInteractionsOptions>({}, { isHidden: true }),
    },
})({
    apply({ params, dependencies }) {
        return Task.create('Compute Interactions', async ctx => {
            const loci: [string, StructureElement.Loci][] = [];
            for (const src of params.sources) {
                const structure = dependencies?.[src.structureRef].data as Structure;
                if (src.bundle) {
                    loci.push([src.structureRef, StructureElement.Bundle.toLoci(src.bundle, structure)]);
                } else {
                    loci.push([src.structureRef, Structure.toStructureElementLoci(structure)]);
                }
            }

            const interactions = await computeInteractions(ctx, loci, params.options);
            return new InteractionData({ interactions }, { label: 'Interactions' });
        });
    }
});