/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { addFixedCountDashedCylinder } from '../../mol-geo/geometry/mesh/builder/cylinder';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { Vec3 } from '../../mol-math/linear-algebra';
import { Shape } from '../../mol-model/shape';
import { Structure, StructureElement } from '../../mol-model/structure';
import { PluginStateObject as SO } from '../../mol-plugin-state/objects';
import { StateTransformer } from '../../mol-state';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { computeInteractions, ComputeInteractionsOptions } from './compute';
import { StructureInteractions } from './model';

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

export const InteractionsShape = Factory({
    name: 'interactions-shape',
    display: { name: 'Interactions Shape' },
    from: InteractionData,
    to: SO.Shape.Provider,
})({
    apply({ a }) {
        return new SO.Shape.Provider({
            label: 'Interactions Shape Provider',
            data: a.data,
            params: PD.withDefaults(Mesh.Params, { }),
            getShape: (_, data, __, prev: any) => buildPrimitiveMesh(data, prev?.geometry),
            geometryUtils: Mesh.Utils,
        }, { label: 'Interactions Shape Provider' });
    }
});

function buildPrimitiveMesh(data: InteractionData['data'], prev?: Mesh): Shape<Mesh> {
    const mesh = MeshBuilder.createState(1024, 1024, prev);

    mesh.currentGroup = -1;
    const tooltips = new Map<number, string>();

    for (const interaction of data.interactions.elements) {
        mesh.currentGroup++;
        tooltips.set(mesh.currentGroup, interaction.schema.kind);

        const a = StructureElement.Loci.getBoundary(interaction.a);
        const b = StructureElement.Loci.getBoundary(interaction.b);

        const radius = 0.1;
        const dist = Vec3.distance(a.sphere.center, b.sphere.center);
        const count = Math.ceil(dist / (2 * radius));
        addFixedCountDashedCylinder(mesh, a.sphere.center, b.sphere.center, 1.0, count, true, {
            radiusBottom: radius,
            radiusTop: radius,
            topCap: true,
            bottomCap: true,
        });
    }

    return Shape.create(
        'Interactions',
        data,
        MeshBuilder.getMesh(mesh),
        (g) => 0x0 as any,
        (g) => 1,
        (g) => tooltips.get(g) ?? ''
    );
}
