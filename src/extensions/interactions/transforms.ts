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
import { computeContacts, ComputeInteractionsOptions } from './compute';
import { getCustomInteractionData } from './custom';
import { InteractionElementSchema, StructureInteractions } from './model';
import { buildInteractionsShape, InteractionVisualParams } from './visuals';

const Factory = StateTransformer.builderFactory('interactions-extension');

export class InteractionData extends SO.Create<{ interactions: StructureInteractions }>({ name: 'Interactions', typeClass: 'Data' }) { }

export interface ComputeInteractionSource {
    structureRef: string,
    bundle?: StructureElement.Bundle,
    // schema?: StructureElementSchema,
}

export const ComputeContacts = Factory({
    name: 'compute-contacts',
    display: 'Compute Contacts',
    from: SO.Root,
    to: InteractionData,
    params: {
        sources: PD.Value<ComputeInteractionSource[]>([], { isHidden: true }),
        options: PD.Value<ComputeInteractionsOptions>({}, { isHidden: true }),
    },
})({
    apply({ params, dependencies }) {
        return Task.create('Compute Contacts', async ctx => {
            const loci: [string, StructureElement.Loci][] = [];
            for (const src of params.sources) {
                const structure = dependencies?.[src.structureRef].data as Structure;
                if (src.bundle) {
                    loci.push([src.structureRef, StructureElement.Bundle.toLoci(src.bundle, structure)]);
                } else {
                    loci.push([src.structureRef, Structure.toStructureElementLoci(structure)]);
                }
            }

            const interactions = await computeContacts(ctx, loci, params.options);
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
        interactions: PD.Value<InteractionElementSchema[]>([], { isHidden: true }),
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
    params: InteractionVisualParams
})({
    canAutoUpdate: () => true,
    apply({ a, params }) {
        return new SO.Shape.Provider({
            label: 'Interactions Shape Provider',
            data: { interactions: a.data.interactions, params },
            params: PD.withDefaults(Mesh.Params, { }),
            getShape: (_, data: { interactions: StructureInteractions, params: InteractionVisualParams }, __, prev: any) => buildInteractionsShape(data.interactions, data.params, prev?.geometry),
            geometryUtils: Mesh.Utils,
        }, { label: 'Interactions Shape Provider' });
    }
});
