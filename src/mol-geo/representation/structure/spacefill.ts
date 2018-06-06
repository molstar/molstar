/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ValueCell } from 'mol-util/value-cell'

import { RenderObject, createMeshRenderObject, MeshRenderObject } from 'mol-gl/render-object'
import { Unit, Element, Queries } from 'mol-model/structure';
import { UnitsRepresentation, DefaultStructureProps } from './index';
import { Task } from 'mol-task'
import { createTransforms, createColors, createFlags, createEmptyFlags, createSphereMesh } from './utils';
import VertexMap from '../../shape/vertex-map';
import { deepEqual, defaults } from 'mol-util';
import { fillSerial } from 'mol-gl/renderable/util';
import { RenderableState, MeshValues } from 'mol-gl/renderable';
import { getMeshData } from '../../util/mesh-data';
import { Mesh } from '../../shape/mesh';
import { PickingId } from '../../util/picking';
import { SortedArray } from 'mol-data/int';

function createSpacefillMesh(unit: Unit, detail: number, mesh?: Mesh) {
    let radius: Element.Property<number>
    if (Unit.isAtomic(unit)) {
        radius = Queries.props.atom.vdw_radius
    } else if (Unit.isSpheres(unit)) {
        radius = Queries.props.coarse.sphere_radius
    } else {
        console.warn('Unsupported unit type')
        return Task.constant('Empty mesh', Mesh.createEmpty(mesh))
    }
    return createSphereMesh(unit, (l) => radius(l) * 0.3, detail, mesh)
}

export const DefaultSpacefillProps = {
    ...DefaultStructureProps,
    flipSided: false,
    flatShaded: false,
    detail: 0,
}
export type SpacefillProps = Partial<typeof DefaultSpacefillProps>

export default function Spacefill(): UnitsRepresentation<SpacefillProps> {
    const renderObjects: RenderObject[] = []
    let spheres: MeshRenderObject
    let currentProps: typeof DefaultSpacefillProps
    let mesh: Mesh
    let currentGroup: Unit.SymmetryGroup
    let vertexMap: VertexMap

    return {
        renderObjects,
        create(group: Unit.SymmetryGroup, props: SpacefillProps = {}) {
            currentProps = Object.assign({}, DefaultSpacefillProps, props)

            return Task.create('Spacefill.create', async ctx => {
                renderObjects.length = 0 // clear
                currentGroup = group

                const { detail, colorTheme, hoverSelection } = { ...DefaultSpacefillProps, ...props }

                mesh = await createSpacefillMesh(group.units[0], detail).runAsChild(ctx, 'Computing spacefill mesh')
                // console.log(mesh)
                vertexMap = VertexMap.fromMesh(mesh)

                await ctx.update('Computing spacefill transforms');
                const transforms = createTransforms(group)

                await ctx.update('Computing spacefill colors');
                const color = createColors(group, vertexMap, colorTheme)

                await ctx.update('Computing spacefill flags');
                const flag = createFlags(group, hoverSelection.instanceId, hoverSelection.elementId)

                const instanceCount = group.units.length

                const values: MeshValues = {
                    ...getMeshData(mesh),
                    aTransform: transforms,
                    aInstanceId: ValueCell.create(fillSerial(new Float32Array(instanceCount))),
                    ...color,
                    ...flag,

                    uAlpha: ValueCell.create(defaults(props.alpha, 1.0)),
                    uInstanceCount: ValueCell.create(instanceCount),
                    uElementCount: ValueCell.create(group.elements.length),

                    elements: mesh.indexBuffer,

                    drawCount: ValueCell.create(mesh.triangleCount * 3),
                    instanceCount: ValueCell.create(instanceCount),

                    dDoubleSided: ValueCell.create(defaults(props.doubleSided, true)),
                    dFlatShaded: ValueCell.create(defaults(props.flatShaded, false)),
                    dFlipSided: ValueCell.create(defaults(props.flipSided, false)),
                }
                const state: RenderableState = {
                    depthMask: defaults(props.depthMask, true),
                    visible: defaults(props.visible, true)
                }

                spheres = createMeshRenderObject(values, state)
                renderObjects.push(spheres)
            })
        },
        update(props: SpacefillProps) {
            const newProps = Object.assign({}, currentProps, props)

            return Task.create('Spacefill.update', async ctx => {
                if (!spheres) return false

                let updateColor = false

                if (newProps.detail !== currentProps.detail) {
                    mesh = await createSpacefillMesh(currentGroup.units[0], newProps.detail, mesh).runAsChild(ctx, 'Computing spacefill mesh')
                    ValueCell.update(spheres.values.drawCount, mesh.triangleCount * 3)
                    // TODO update in-place
                    vertexMap = VertexMap.fromMesh(mesh)
                    updateColor = true
                }

                if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) {
                    updateColor = true
                }

                if (updateColor) {
                    await ctx.update('Computing spacefill colors');
                    createColors(currentGroup, vertexMap, newProps.colorTheme, spheres.values)
                }

                if (newProps.hoverSelection !== currentProps.hoverSelection) {
                    await ctx.update('Computing spacefill flags');
                    if (newProps.hoverSelection.objectId === spheres.id) {
                        createFlags(currentGroup, newProps.hoverSelection.instanceId, newProps.hoverSelection.elementId, spheres.values)
                    } else {
                        createEmptyFlags(spheres.values)
                    }
                }

                ValueCell.updateIfChanged(spheres.values.uAlpha, newProps.alpha)
                ValueCell.updateIfChanged(spheres.values.dDoubleSided, newProps.doubleSided)
                ValueCell.updateIfChanged(spheres.values.dFlipSided, newProps.flipSided)
                ValueCell.updateIfChanged(spheres.values.dFlatShaded, newProps.flatShaded)

                spheres.state.visible = newProps.visible
                spheres.state.depthMask = newProps.depthMask

                currentProps = newProps
                return true
            })
        },
        getLoci(pickingId: PickingId) {
            const { objectId, instanceId, elementId } = pickingId
            if (spheres.id === objectId) {
                const unit = currentGroup.units[instanceId]
                const elements = SortedArray.ofSingleton(currentGroup.elements[elementId])
                return Element.Loci([{ unit, elements }])
            }
            return null
        }
    }
}
