/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Unit, Element } from 'mol-model/structure';
import { DefaultStructureProps, UnitsVisual } from '../index';
import { RuntimeContext } from 'mol-task'
import { createTransforms, createColors } from './util/common';
import { markElement } from './util/element';
import { deepEqual } from 'mol-util';
import { MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../util/mesh-data';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { OrderedSet } from 'mol-data/int';
import { createMarkers, MarkerAction } from '../../../util/marker-data';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { SizeTheme } from '../../../theme';
import { createMeshValues, updateMeshValues, updateRenderableState, createRenderableState, DefaultMeshProps } from '../../util';
import { MeshBuilder } from '../../../shape/mesh-builder';
import { getPolymerElementCount, PolymerTraceIterator } from './util/polymer';
import { Vec3 } from 'mol-math/linear-algebra';

// export function spline(target: THREE.Vector3, p1: THREE.Vector3, p2: THREE.Vector3, p3: THREE.Vector3, t: number) {
//     let a = Math.pow(1 - t, 2) / 2;
//     let c = Math.pow(t, 2) / 2;
//     let b = 1 - a - c;

//     let x = a * p1.x + b * p2.x + c * p3.x;
//     let y = a * p1.y + b * p2.y + c * p3.y;
//     let z = a * p1.z + b * p2.z + c * p3.z;

//     target.set(x, y, z);
// }

async function createPolymerTraceMesh(ctx: RuntimeContext, unit: Unit, mesh?: Mesh) {
    const polymerElementCount = getPolymerElementCount(unit)
    console.log('polymerElementCount', polymerElementCount)
    if (!polymerElementCount) return Mesh.createEmpty(mesh)

    // TODO better vertex count estimates
    const builder = MeshBuilder.create(polymerElementCount * 30, polymerElementCount * 30 / 2, mesh)
    const linearSegmentCount = 10

    const v0 = Vec3.zero()
    const v1 = Vec3.zero()

    let i = 0
    const polymerTraceIt = PolymerTraceIterator(unit)
    while (polymerTraceIt.hasNext) {
        const v = polymerTraceIt.move()

        Vec3.spline(v1, v.c0, v.c1, v.c2, v.c3, 0.5, 0.5)

        builder.setId(v.index)
        for (let j = 1; j <= linearSegmentCount; ++j) {
            let t = j * 1.0 / linearSegmentCount;
            Vec3.copy(v0, v1)
            // if ((v.last && t > 0.5) || (v.first && t < 0.5)) break
            if (t < 0.5) {
                Vec3.spline(v1, v.c0, v.c1, v.c2, v.c3, t + 0.5, 0.5)
            } else {
                Vec3.spline(v1, v.c1, v.c2, v.c3, v.c4, t - 0.5, 0.5)
            }
            // TODO size theme
            builder.addCylinder(v0, v1, 1.0, { radiusTop: 0.1, radiusBottom: 0.1 })
        }

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Backbone mesh', current: i, max: polymerElementCount });
        }
        ++i
    }

    return builder.getMesh()
}

export const DefaultPolymerTraceProps = {
    ...DefaultMeshProps,
    ...DefaultStructureProps,
    sizeTheme: { name: 'physical', factor: 1 } as SizeTheme,
    detail: 0,
    unitKinds: [ Unit.Kind.Atomic, Unit.Kind.Spheres ] as Unit.Kind[]
}
export type PolymerTraceProps = Partial<typeof DefaultPolymerTraceProps>

export function PolymerTraceVisual(): UnitsVisual<PolymerTraceProps> {
    let renderObject: MeshRenderObject
    let currentProps: typeof DefaultPolymerTraceProps
    let mesh: Mesh
    let currentGroup: Unit.SymmetryGroup

    return {
        get renderObject () { return renderObject },
        async create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: PolymerTraceProps = {}) {
            currentProps = Object.assign({}, DefaultPolymerTraceProps, props)
            currentGroup = group

            const { colorTheme, unitKinds } = { ...DefaultPolymerTraceProps, ...props }
            const instanceCount = group.units.length
            const elementCount = group.elements.length
            const unit = group.units[0]

            mesh = unitKinds.includes(unit.kind)
                ? await createPolymerTraceMesh(ctx, unit, mesh)
                : Mesh.createEmpty(mesh)

            const transforms = createTransforms(group)
            const color = createColors(group, elementCount, colorTheme)
            const marker = createMarkers(instanceCount * elementCount)

            const counts = { drawCount: mesh.triangleCount * 3, elementCount, instanceCount }

            const values: MeshValues = {
                ...getMeshData(mesh),
                ...color,
                ...marker,
                aTransform: transforms,
                elements: mesh.indexBuffer,
                ...createMeshValues(currentProps, counts),
            }
            const state = createRenderableState(currentProps)

            renderObject = createMeshRenderObject(values, state)
        },
        async update(ctx: RuntimeContext, props: PolymerTraceProps) {
            const newProps = Object.assign({}, currentProps, props)

            if (!renderObject) return false

            let updateColor = false

            if (newProps.detail !== currentProps.detail) {
                const unit = currentGroup.units[0]
                mesh = await createPolymerTraceMesh(ctx, unit, mesh)
                ValueCell.update(renderObject.values.drawCount, mesh.triangleCount * 3)
                updateColor = true
            }

            if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) {
                updateColor = true
            }

            if (updateColor) {
                const elementCount = currentGroup.elements.length
                if (ctx.shouldUpdate) await ctx.update('Computing trace colors');
                createColors(currentGroup, elementCount, newProps.colorTheme, renderObject.values)
            }

            updateMeshValues(renderObject.values, newProps)
            updateRenderableState(renderObject.state, newProps)

            currentProps = newProps
            return true
        },
        getLoci(pickingId: PickingId) {
            const { objectId, instanceId, elementId } = pickingId
            if (renderObject.id === objectId) {
                const unit = currentGroup.units[instanceId]
                const indices = OrderedSet.ofSingleton(elementId as Element.Index);
                return Element.Loci([{ unit, indices }])
            }
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            markElement(renderObject.values.tMarker, currentGroup, loci, action)
        },
        destroy() {
            // TODO
        }
    }
}
