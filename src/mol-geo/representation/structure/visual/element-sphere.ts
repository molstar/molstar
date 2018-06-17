/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ValueCell } from 'mol-util/value-cell'

import { RenderObject, createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Unit, Element } from 'mol-model/structure';
import { DefaultStructureProps, UnitsVisual } from '../index';
import { RuntimeContext } from 'mol-task'
import { createTransforms, createColors, createElementSphereMesh, markElement, getElementRadius } from '../utils';
import VertexMap from '../../../shape/vertex-map';
import { deepEqual, defaults } from 'mol-util';
import { fillSerial } from 'mol-gl/renderable/util';
import { RenderableState, MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../../util/mesh-data';
import { Mesh } from '../../../shape/mesh';
import { PickingId } from '../../../util/picking';
import { OrderedSet } from 'mol-data/int';
import { createMarkers, MarkerAction } from '../../../util/marker-data';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { SizeTheme } from '../../../theme';

export const DefaultElementSphereProps = {
    ...DefaultStructureProps,
    sizeTheme: { name: 'physical', factor: 1 } as SizeTheme,
    flipSided: false,
    flatShaded: false,
    detail: 0,
}
export type ElementSphereProps = Partial<typeof DefaultElementSphereProps>

export function ElementSphereVisual(): UnitsVisual<ElementSphereProps> {
    const renderObjects: RenderObject[] = []
    let spheres: MeshRenderObject
    let currentProps: typeof DefaultElementSphereProps
    let mesh: Mesh
    let currentGroup: Unit.SymmetryGroup
    let vertexMap: VertexMap

    return {
        renderObjects,
        async create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: ElementSphereProps = {}) {
            currentProps = Object.assign({}, DefaultElementSphereProps, props)

            renderObjects.length = 0 // clear
            currentGroup = group

            const { detail, colorTheme, sizeTheme } = { ...DefaultElementSphereProps, ...props }
            const instanceCount = group.units.length
            const elementCount = group.elements.length
            const unit = group.units[0]

            const radius = getElementRadius(unit, sizeTheme)
            mesh = await createElementSphereMesh(ctx, unit, radius, detail, mesh)
            // console.log(mesh)
            vertexMap = VertexMap.fromMesh(mesh)

            if (ctx.shouldUpdate) await ctx.update('Computing spacefill transforms');
            const transforms = createTransforms(group)

            if (ctx.shouldUpdate) await ctx.update('Computing spacefill colors');
            const color = createColors(group, vertexMap, colorTheme)

            if (ctx.shouldUpdate) await ctx.update('Computing spacefill marks');
            const marker = createMarkers(instanceCount * elementCount)

            const values: MeshValues = {
                ...getMeshData(mesh),
                aTransform: transforms,
                aInstanceId: ValueCell.create(fillSerial(new Float32Array(instanceCount))),
                ...color,
                ...marker,

                uAlpha: ValueCell.create(defaults(props.alpha, 1.0)),
                uInstanceCount: ValueCell.create(instanceCount),
                uElementCount: ValueCell.create(elementCount),

                elements: mesh.indexBuffer,

                drawCount: ValueCell.create(mesh.triangleCount * 3),
                instanceCount: ValueCell.create(instanceCount),

                dDoubleSided: ValueCell.create(defaults(props.doubleSided, true)),
                dFlatShaded: ValueCell.create(defaults(props.flatShaded, false)),
                dFlipSided: ValueCell.create(defaults(props.flipSided, false)),
                dUseFog: ValueCell.create(defaults(props.useFog, true)),
            }
            const state: RenderableState = {
                depthMask: defaults(props.depthMask, true),
                visible: defaults(props.visible, true)
            }

            spheres = createMeshRenderObject(values, state)
            renderObjects.push(spheres)
        },
        async update(ctx: RuntimeContext, props: ElementSphereProps) {
            const newProps = Object.assign({}, currentProps, props)

            if (!spheres) return false

            let updateColor = false

            if (newProps.detail !== currentProps.detail) {
                const unit = currentGroup.units[0]
                const radius = getElementRadius(unit, newProps.sizeTheme)
                mesh = await createElementSphereMesh(ctx, unit, radius, newProps.detail, mesh)
                ValueCell.update(spheres.values.drawCount, mesh.triangleCount * 3)
                // TODO update in-place
                vertexMap = VertexMap.fromMesh(mesh)
                updateColor = true
            }

            if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) {
                updateColor = true
            }

            if (updateColor) {
                if (ctx.shouldUpdate) await ctx.update('Computing spacefill colors');
                createColors(currentGroup, vertexMap, newProps.colorTheme, spheres.values)
            }

            ValueCell.updateIfChanged(spheres.values.uAlpha, newProps.alpha)
            ValueCell.updateIfChanged(spheres.values.dDoubleSided, newProps.doubleSided)
            ValueCell.updateIfChanged(spheres.values.dFlipSided, newProps.flipSided)
            ValueCell.updateIfChanged(spheres.values.dFlatShaded, newProps.flatShaded)

            spheres.state.visible = newProps.visible
            spheres.state.depthMask = newProps.depthMask

            currentProps = newProps
            return true
        },
        getLoci(pickingId: PickingId) {
            const { objectId, instanceId, elementId } = pickingId
            if (spheres.id === objectId) {
                const unit = currentGroup.units[instanceId]
                const indices = OrderedSet.ofSingleton(elementId as Element.Index);
                return Element.Loci([{ unit, indices }])
            }
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            markElement(spheres.values.tMarker, currentGroup, loci, action)
        },
        destroy() {
            // TODO
        }
    }
}
