/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { Structure, StructureElement } from '../../mol-model/structure';
import { PluginStateObject as SO } from '../../mol-plugin-state/objects';
import { StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { computeInteractions, ComputeInteractionsOptions } from './compute';
import { getCustomInteractionData } from './custom';
import { InteractionSchema, StructureInteractions } from './model';
import { buildInteractionsShape } from './visuals';

const Factory = StateTransformer.builderFactory('interactions-extension');

export class InteractionData extends SO.Create<{ interactions: StructureInteractions }>({ name: 'Interactions', typeClass: 'Data' }) { }

export interface ComputeInteractionSource {
    structureRef: string,
    bundle?: StructureElement.Bundle,
    // schema?: StructureElementSchema,
}

export const ComputeInteractions = Factory({
    name: 'compute-interactions',
    display: 'Compute Interactions',
    from: SO.Root,
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

export const CustomInteractions = Factory({
    name: 'custom-interactions',
    display: 'Custom Interactions',
    from: SO.Root,
    to: InteractionData,
    params: {
        interactions: PD.Value<InteractionSchema[]>([], { isHidden: true }),
    },
})({
    apply({ params, dependencies }) {
        return Task.create('Custom Interactions', async ctx => {
            const structures: { [ref: string]: Structure } = {};
            for (const [k, v] of Object.entries(dependencies ?? {})) {
                structures[k] = v.data as Structure;
            }
            const interactions = getCustomInteractionData(params.interactions, structures);
            return new InteractionData({ interactions }, { label: 'Custom Interactions' });
        });
    }
});

export const InteractionsShape = Factory({
    name: 'interactions-shape',
    display: { name: 'Interactions Shape' },
    from: InteractionData,
    to: SO.Shape.Provider,
})({
    apply({ a }) {
        return new SO.Shape.Provider({
            label: 'Interactions Shape Provider',
            data: a.data.interactions,
            params: PD.withDefaults(Mesh.Params, { }),
            getShape: (_, data: StructureInteractions, __, prev: any) => buildInteractionsShape(data, prev?.geometry),
            geometryUtils: Mesh.Utils,
        }, { label: 'Interactions Shape Provider' });
    }
});
