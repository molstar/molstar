/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from 'mol-math/linear-algebra';
import { Unit, Element } from 'mol-model/structure';
import { SizeTheme } from '../../../../theme';
import { RuntimeContext } from 'mol-task';
import { icosahedronVertexCount } from '../../../../primitive/icosahedron';
import { Mesh } from '../../../../shape/mesh';
import { MeshBuilder } from '../../../../shape/mesh-builder';
import { ValueCell, defaults } from 'mol-util';
import { TextureImage } from 'mol-gl/renderable/util';
import { Loci, isEveryLoci } from 'mol-model/loci';
import { MarkerAction, applyMarkerAction } from '../../../../util/marker-data';
import { Interval } from 'mol-data/int';
import { getPhysicalRadius } from '../../../../theme/structure/size/physical';

export function getElementRadius(unit: Unit, props: SizeTheme): Element.Property<number> {
    switch (props.name) {
        case 'uniform':
            return () => props.value
        case 'physical':
            const radius = getPhysicalRadius(unit)
            const factor = defaults(props.factor, 1)
            return (l) => radius(l) * factor
    }
}

export async function createElementSphereMesh(ctx: RuntimeContext, unit: Unit, radius: Element.Property<number>, detail: number, mesh?: Mesh) {
    const { elements } = unit;
    const elementCount = elements.length;
    const vertexCount = elementCount * icosahedronVertexCount(detail)
    const meshBuilder = MeshBuilder.create(vertexCount, vertexCount / 2, mesh)

    const v = Vec3.zero()
    const pos = unit.conformation.invariantPosition
    const l = Element.Location()
    l.unit = unit

    for (let i = 0; i < elementCount; i++) {
        l.element = elements[i]
        pos(elements[i], v)

        meshBuilder.setId(i)
        meshBuilder.addIcosahedron(v, radius(l), detail)

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Sphere mesh', current: i, max: elementCount });
        }
    }

    return meshBuilder.getMesh()
}


export function markElement(tMarker: ValueCell<TextureImage>, group: Unit.SymmetryGroup, loci: Loci, action: MarkerAction) {
    let changed = false
    const elementCount = group.elements.length
    const instanceCount = group.units.length
    const array = tMarker.ref.value.array
    if (isEveryLoci(loci)) {
        applyMarkerAction(array, 0, elementCount * instanceCount, action)
        changed = true
    } else if (Element.isLoci(loci)) {
        for (const e of loci.elements) {
            const unitIdx = Unit.findUnitById(e.unit.id, group.units)
            if (unitIdx !== -1) {
                if (Interval.is(e.indices)) {
                    const idxStart = unitIdx * elementCount + Interval.start(e.indices);
                    const idxEnd = unitIdx * elementCount + Interval.end(e.indices);
                    if (applyMarkerAction(array, idxStart, idxEnd, action) && !changed) {
                        changed = true
                    }
                } else {
                    for (let i = 0, _i = e.indices.length; i < _i; i++) {
                        const idx = unitIdx * elementCount + e.indices[i];
                        if (applyMarkerAction(array, idx, idx + 1, action) && !changed) {
                            changed = true
                        }
                    }
                }
            }
        }
    } else {
        return
    }
    if (changed) {
        ValueCell.update(tMarker, tMarker.ref.value)
    }
}
