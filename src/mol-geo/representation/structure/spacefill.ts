/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createRenderObject, RenderObject } from 'mol-gl/renderer'
import { createColorTexture } from 'mol-gl/util';
import { Vec3, Mat4 } from 'mol-math/linear-algebra'
import { OrderedSet } from 'mol-data/int'
import { Element, Unit, ElementSet } from 'mol-model/structure';
import P from 'mol-model/structure/query/properties';
import { RepresentationProps, UnitRepresentation } from './index';
import { Task } from 'mol-task'
import { MeshBuilder } from '../../shape/mesh-builder';

export default function Spacefill(): UnitRepresentation {
    const renderObjects: RenderObject[] = []

    // unit: Unit, atomGroup: AtomGroup

    return {
        create: (units: ReadonlyArray<Unit>, elements: ElementSet, props: Partial<RepresentationProps> = {}) => Task.create('Spacefill', async ctx => {
            const l = Element.Location();
            const meshBuilder = MeshBuilder.create()

            const unitIds = ElementSet.unitIds(elements);
            for (let i = 0, _i = unitIds.length; i < _i; i++) {
                const unitId = unitIds[i];
                const unit = units[unitId];
                const elementGroup = ElementSet.unitGetByIndex(elements, i);

                const elementCount = OrderedSet.size(elementGroup.elements)

                l.unit = unit;

                const v = Vec3.zero()
                const m = Mat4.identity()

                for (let i = 0; i < elementCount; i++) {
                    l.element = OrderedSet.getAt(elementGroup.elements, i)

                    v[0] = P.atom.x(l)
                    v[1] = P.atom.y(l)
                    v[2] = P.atom.z(l)
                    Mat4.setTranslation(m, v)

                    meshBuilder.addIcosahedron(m, { radius: P.atom.vdw(l), detail: 1 })
                }

                if (i % 10 === 0 && ctx.shouldUpdate) {
                    await ctx.update({ message: 'Spacefill', current: i, max: _i });
                }
            }

            const transformArray = new Float32Array(32)
            const m4 = Mat4.identity()
            Mat4.toArray(m4, transformArray, 0)

            const color = ValueCell.create(createColorTexture(1))
            color.ref.value.set([ 0, 0, 255 ])

            const mesh = meshBuilder.getMesh()

            // console.log(mesh)

            const spheres = createRenderObject('mesh', {
                position: mesh.vertexBuffer,
                normal: mesh.normalBuffer,
                color,
                transform: ValueCell.create(transformArray),
                elements: mesh.indexBuffer,

                instanceCount: transformArray.length / 16,
                elementCount: mesh.triangleCount,
                positionCount: mesh.vertexCount
            }, {})
            renderObjects.push(spheres)

            return renderObjects
        }),
        update: (props: RepresentationProps) => false
    }
}
