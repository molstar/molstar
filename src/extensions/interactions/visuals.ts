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
import { StructureElement } from '../../mol-model/structure';
import { StructureInteractions } from './model';

export function buildInteractionsShape(interactions: StructureInteractions, prev?: Mesh): Shape<Mesh> {
    const mesh = MeshBuilder.createState(1024, 1024, prev);

    mesh.currentGroup = -1;
    const tooltips = new Map<number, string>();

    for (const interaction of interactions.elements) {
        mesh.currentGroup++;
        let tooltip = interaction.info.kind;
        if (interaction.sourceSchema?.description) {
            tooltip += ` (${interaction.sourceSchema.description})`;
        }
        tooltips.set(mesh.currentGroup, tooltip);

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
        interactions,
        MeshBuilder.getMesh(mesh),
        (g) => 0x0 as any,
        (g) => 1,
        (g) => tooltips.get(g) ?? ''
    );
}
