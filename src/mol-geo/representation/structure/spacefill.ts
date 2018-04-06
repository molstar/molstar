/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { createRenderObject, RenderObject } from 'mol-gl/renderer'
import { createColorTexture } from 'mol-gl/util';
import Icosahedron from 'mol-geo/primitive/icosahedron'
import { Vec3, Mat4 } from 'mol-math/linear-algebra'
import { OrderedSet } from 'mol-data/int'
import { Atom, AtomGroup, Unit } from 'mol-model/structure';
import P from 'mol-model/structure/query/properties';
import { RepresentationProps, UnitRepresentation } from './index';
import { Task } from 'mol-task'

export default function Spacefill(): UnitRepresentation {
    let vertices: Float32Array
    let normals: Float32Array

    const renderObjects: RenderObject[] = []

    // unit: Unit, atomGroup: AtomGroup

    return {
        create: (unit: Unit, atomGroup: AtomGroup, props: Partial<RepresentationProps> = {}) => Task.create('Spacefill', async ctx => {
            const atomCount = OrderedSet.size(atomGroup.atoms)

            const l = Atom.Location();
            l.unit = unit;

            const sphere = Icosahedron({ radius: 1, detail: 0 })
            const vertexCount = sphere.vertices.length / 3

            vertices = new Float32Array(atomCount * vertexCount * 3)
            normals = new Float32Array(atomCount * vertexCount * 3)

            const v = Vec3.zero()
            const m = Mat4.identity()

            for (let i = 0; i < atomCount; i++) {
                l.atom = OrderedSet.getAt(atomGroup.atoms, i)

                v[0] = P.atom.x(l)
                v[1] = P.atom.y(l)
                v[2] = P.atom.z(l)
                Mat4.setTranslation(m, v)

                for (let j = 0; j < vertexCount; ++j) {
                    Vec3.fromArray(v, sphere.vertices, j * 3)
                    Vec3.transformMat4(v, v, m)
                    Vec3.toArray(v, vertices, i * vertexCount * 3 + j * 3)
                }

                normals.set(sphere.normals, i * vertexCount * 3);

                if (i % 100 === 0 && ctx.shouldUpdate) {
                    await ctx.update({ message: 'Spacefill', current: i, max: atomCount });
                }
            }

            const transformArray = new Float32Array(16)
            const m4 = Mat4.identity()
            Mat4.toArray(m4, transformArray, 0)

            const color = ValueCell.create(createColorTexture(1))
            color.ref.value.set([ 0, 0, 255 ])

            const spheres = createRenderObject(
                'mesh',
                {
                    position: ValueCell.create(new Float32Array(vertices)),
                    normal: ValueCell.create(new Float32Array(normals)),
                    color,
                    transform: ValueCell.create(transformArray)
                }
            )

            // console.log({ vertices, normals, vertexCount, atomCount })

            renderObjects.push(spheres)

            return renderObjects
        }),
        update: (props: RepresentationProps) => false
    }
}
