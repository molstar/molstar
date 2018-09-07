/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ValueCell } from 'mol-util/value-cell'
import { PointRenderObject } from 'mol-gl/render-object'
import { Unit, Structure } from 'mol-model/structure';
import { RuntimeContext } from 'mol-task'
import { UnitsVisual, DefaultStructureProps } from '..';
import { getElementLoci, StructureElementIterator, markElement } from './util/element';
import { createColors, createSizes, createUnitsPointRenderObject } from './util/common';
import { deepEqual, UUID } from 'mol-util';
import { Interval } from 'mol-data/int';
import { PickingId } from '../../../util/picking';
import { Loci, EmptyLoci, isEveryLoci } from 'mol-model/loci';
import { MarkerAction, createMarkers, applyMarkerAction } from '../../../util/marker-data';
import { Vec3 } from 'mol-math/linear-algebra';
import { fillSerial } from 'mol-util/array';
import { SizeThemeProps } from 'mol-view/theme/size';
import { LocationIterator } from '../../../util/location-iterator';
import { createTransforms } from '../../../util/transform-data';
import { StructureGroup } from '../units-visual';
import { updateRenderableState } from '../../util';

export const DefaultElementPointProps = {
    ...DefaultStructureProps,

    sizeTheme: { name: 'uniform', value: 0.2 } as SizeThemeProps,
    pointSizeAttenuation: true,
}
export type ElementPointProps = Partial<typeof DefaultElementPointProps>

export async function createElementPointVertices(ctx: RuntimeContext, unit: Unit, vertices?: ValueCell<Float32Array>) {
    const elements = unit.elements
    const n = elements.length * 3
    const array = vertices && vertices.ref.value.length >= n ? vertices.ref.value : new Float32Array(n)

    const pos = unit.conformation.invariantPosition

    const p = Vec3.zero()
    for (let i = 0; i < n; i += 3) {
        pos(elements[i / 3], p)
        array[i] = p[0]
        array[i + 1] = p[1]
        array[i + 2] = p[2]

        if (i % 10000 === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Creating points', current: i / 3, max: elements.length });
        }
        ++i
    }
    return vertices ? ValueCell.update(vertices, array) : ValueCell.create(array)
}

export function ElementPointVisual(): UnitsVisual<ElementPointProps> {
    let renderObject: PointRenderObject | undefined
    let currentProps = DefaultElementPointProps
    let currentGroup: Unit.SymmetryGroup
    let currentStructure: Structure
    let locationIt: LocationIterator
    let vertices: ValueCell<Float32Array>
    let currentConformationId: UUID

    return {
        get renderObject () { return renderObject },
        async createOrUpdate(ctx: RuntimeContext, props: ElementPointProps = {}, structureGroup?: StructureGroup) {
            if (structureGroup) currentStructure = structureGroup.structure
            const group = structureGroup ? structureGroup.group : undefined
            if (!group && !currentGroup) {
                throw new Error('missing group')
            } else if (group && !currentGroup) {
                currentProps = Object.assign({}, DefaultElementPointProps, props)
                currentProps.colorTheme.structure = currentStructure
                currentGroup = group
                locationIt = StructureElementIterator.fromGroup(group)

                const unit = group.units[0]
                currentConformationId = Unit.conformationId(unit)
                vertices = await createElementPointVertices(ctx, unit, vertices)

                renderObject = await createUnitsPointRenderObject(ctx, group, vertices, locationIt, currentProps)
            } else if (renderObject) {
                if (group) currentGroup = group

                const newProps = { ...currentProps, ...props }
                const unit = currentGroup.units[0]

                let updateTransform = false
                let createVertices = false
                let updateColor = false
                let updateSize = false

                const newConformationId = Unit.conformationId(unit)
                if (newConformationId !== currentConformationId) {
                    currentConformationId = newConformationId
                    createVertices = true
                }

                if (currentGroup.units.length !== locationIt.instanceCount) updateTransform = true

                if (!deepEqual(newProps.sizeTheme, currentProps.sizeTheme)) createVertices = true
                if (!deepEqual(newProps.colorTheme, currentProps.colorTheme)) updateColor = true
                if (!deepEqual(newProps.sizeTheme, currentProps.sizeTheme)) updateSize = true

                if (updateTransform) {
                    locationIt = StructureElementIterator.fromGroup(currentGroup)
                    const { instanceCount, groupCount } = locationIt
                    createTransforms(currentGroup, renderObject.values)
                    createMarkers(instanceCount * groupCount, renderObject.values)
                    updateColor = true
                    updateSize = true
                }

                if (createVertices) {
                    await createElementPointVertices(ctx, unit, vertices)
                    ValueCell.update(renderObject.values.aGroup, fillSerial(new Float32Array(locationIt.groupCount))) // TODO reuse array
                    ValueCell.update(renderObject.values.drawCount, locationIt.groupCount)
                    updateColor = true
                    updateSize = true
                }

                if (updateColor) {
                    await createColors(ctx, locationIt, newProps.colorTheme, renderObject.values)
                }

                if (updateSize) {
                    await createSizes(ctx, locationIt, newProps.sizeTheme, renderObject.values)
                }

                updateRenderableState(renderObject.state, newProps)

                currentProps = newProps
            }
        },
        getLoci(pickingId: PickingId) {
            return renderObject ? getElementLoci(pickingId, currentGroup, renderObject.id) : EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            if (!renderObject) return false
            const { tMarker } = renderObject.values
            const { groupCount, instanceCount } = locationIt

            function apply(interval: Interval) {
                const start = Interval.start(interval)
                const end = Interval.end(interval)
                return applyMarkerAction(tMarker.ref.value.array, start, end, action)
            }

            let changed = false
            if (isEveryLoci(loci)) {
                apply(Interval.ofBounds(0, groupCount * instanceCount))
                changed = true
            } else {
                changed = markElement(loci, currentGroup, apply)
            }
            if (changed) {
                ValueCell.update(tMarker, tMarker.ref.value)
            }
            return changed
        },
        destroy() {
            // TODO
            renderObject = undefined
        }
    }
}
