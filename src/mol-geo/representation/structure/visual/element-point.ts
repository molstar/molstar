/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ValueCell } from 'mol-util/value-cell'
import { createPointRenderObject, RenderObject, PointRenderObject } from 'mol-gl/render-object'
import { Unit, Element } from 'mol-model/structure';
import { RuntimeContext } from 'mol-task'
import { fillSerial } from 'mol-gl/renderable/util';

import { UnitsVisual, DefaultStructureProps } from '../index';
import VertexMap from '../../../shape/vertex-map';
import { SizeTheme } from '../../../theme';
import { createTransforms, createColors, createSizes, markElement } from '../utils';
import { deepEqual, defaults } from 'mol-util';
import { SortedArray, OrderedSet } from 'mol-data/int';
import { RenderableState, PointValues } from 'mol-gl/renderable';
import { PickingId } from '../../../util/picking';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { MarkerAction, createMarkers } from '../../../util/marker-data';
import { Vec3 } from 'mol-math/linear-algebra';

export const DefaultPointProps = {
    ...DefaultStructureProps,
    sizeTheme: { name: 'physical' } as SizeTheme
}
export type PointProps = Partial<typeof DefaultPointProps>

export function createPointVertices(unit: Unit) {
    const elements = unit.elements
    const elementCount = elements.length
    const vertices = new Float32Array(elementCount * 3)

    const pos = unit.conformation.invariantPosition

    const p = Vec3.zero()
    for (let i = 0; i < elementCount; i++) {
        const i3 = i * 3
        pos(elements[i], p)
        vertices[i3] = p[0]
        vertices[i3 + 1] = p[1]
        vertices[i3 + 2] = p[2]
    }
    return vertices
}

export default function PointVisual(): UnitsVisual<PointProps> {
    const renderObjects: RenderObject[] = []
    let points: PointRenderObject
    let currentProps = DefaultPointProps
    let currentGroup: Unit.SymmetryGroup

    let _units: ReadonlyArray<Unit>
    let _elements: SortedArray

    return {
        renderObjects,
        async create(ctx: RuntimeContext, group: Unit.SymmetryGroup, props: PointProps = {}) {
            currentProps = Object.assign({}, DefaultPointProps, props)

            renderObjects.length = 0 // clear
            currentGroup = group

            _units = group.units
            _elements = group.elements;

            const { colorTheme, sizeTheme } = currentProps
            const elementCount = _elements.length
            const instanceCount = group.units.length

            const vertexMap = VertexMap.create(
                elementCount,
                elementCount + 1,
                fillSerial(new Uint32Array(elementCount)),
                fillSerial(new Uint32Array(elementCount + 1))
            )

            if (ctx.shouldUpdate) await ctx.update('Computing point vertices');
            const vertices = createPointVertices(_units[0])

            if (ctx.shouldUpdate) await ctx.update('Computing point transforms');
            const transforms = createTransforms(group)

            if (ctx.shouldUpdate) await ctx.update('Computing point colors');
            const color = createColors(group, elementCount, colorTheme)

            if (ctx.shouldUpdate) await ctx.update('Computing point sizes');
            const size = createSizes(group, vertexMap, sizeTheme)

            if (ctx.shouldUpdate) await ctx.update('Computing spacefill marks');
            const marker = createMarkers(instanceCount * elementCount)

            const values: PointValues = {
                aPosition: ValueCell.create(vertices),
                aElementId: ValueCell.create(fillSerial(new Float32Array(elementCount))),
                aTransform: transforms,
                aInstanceId: ValueCell.create(fillSerial(new Float32Array(instanceCount))),
                ...color,
                ...marker,
                ...size,

                uAlpha: ValueCell.create(defaults(props.alpha, 1.0)),
                uInstanceCount: ValueCell.create(instanceCount),
                uElementCount: ValueCell.create(group.elements.length),

                drawCount: ValueCell.create(vertices.length / 3),
                instanceCount: ValueCell.create(instanceCount),

                dPointSizeAttenuation: ValueCell.create(true),
                dUseFog: ValueCell.create(defaults(props.useFog, true)),
            }
            const state: RenderableState = {
                depthMask: defaults(props.depthMask, true),
                visible: defaults(props.visible, true)
            }

            points = createPointRenderObject(values, state)
            renderObjects.push(points)
        },
        async update(ctx: RuntimeContext, props: PointProps) {
            if (!points || !_units || !_elements) return false

            const newProps = { ...currentProps, ...props }
            if (deepEqual(currentProps, newProps)) {
                console.log('props identical, nothing to change')
                return true
            }

            // const elementCount = OrderedSet.size(_elementGroup.elements)
            // const unitCount = _units.length

            // const vertexMap = VertexMap.create(
            //     elementCount,
            //     elementCount + 1,
            //     fillSerial(new Uint32Array(elementCount)),
            //     fillSerial(new Uint32Array(elementCount + 1))
            // )

            if (!deepEqual(currentProps.colorTheme, newProps.colorTheme)) {
                console.log('colorTheme changed', currentProps.colorTheme, newProps.colorTheme)
                // if (ctx.shouldUpdate) await ctx.update('Computing point colors');
                // const color = createColors(_units, _elementGroup, vertexMap, newProps.colorTheme)
                // ValueCell.update(points.props.color, color)
            }

            if (!deepEqual(currentProps.sizeTheme, newProps.sizeTheme)) {
                console.log('sizeTheme changed', currentProps.sizeTheme, newProps.sizeTheme)
            }

            currentProps = newProps
            return false
        },
        getLoci(pickingId: PickingId) {
            const { objectId, instanceId, elementId } = pickingId
            if (points.id === objectId) {
                const unit = currentGroup.units[instanceId]
                const indices = OrderedSet.ofSingleton(elementId as Element.Index)
                return Element.Loci([{ unit, indices }])
            }
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            markElement(points.values.tMarker, currentGroup, loci, action)
        },
        destroy() {
            // TODO
        }
    }
}
