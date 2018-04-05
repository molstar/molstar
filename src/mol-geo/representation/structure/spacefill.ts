/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { ValueBox } from 'mol-util/value-cell'

import { MeshRenderable, Renderable } from 'mol-gl/renderable'
import { calculateTextureInfo } from 'mol-gl/util';
import Icosahedron from 'mol-geo/primitive/icosahedron'
import { Vec3, Mat4 } from 'mol-math/linear-algebra'
import { OrderedSet } from 'mol-data/int'
import { Atom, AtomGroup, Unit } from 'mol-model/structure';
import P from 'mol-model/structure/query/properties';
import { RepresentationProps, UnitRepresentation } from './index';

export default function Spacefill(regl: REGL.Regl): UnitRepresentation {
    let vertices: Float32Array
    let normals: Float32Array

    const renderables: Renderable<any>[] = []

    return {
        create: (unit: Unit, atomGroup: AtomGroup, props: RepresentationProps) => {
            const atomCount = OrderedSet.size(atomGroup.atoms)

            const l = Atom.Location();
            l.unit = unit;

            const sphere = Icosahedron(1, 0)
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

                normals.set(sphere.normals, i * vertexCount * 3)
            }

            const transformArray = new Float32Array(16)
            const m4 = Mat4.identity()
            Mat4.toArray(m4, transformArray, 0)

            const colorTexInfo = calculateTextureInfo(3, 3)
            const color = new Uint8Array(colorTexInfo.length)
            color.set([
                0, 0, 255,
                0, 255, 0,
                255, 0, 0
            ])
            // console.log(color, colorTexInfo)
            const colorTex = regl.texture({
                width: colorTexInfo.width,
                height: colorTexInfo.height,
                format: 'rgb',
                type: 'uint8',
                wrapS: 'clamp',
                wrapT: 'clamp',
                data: color
            })

            const spheres = MeshRenderable.create(regl,
                {
                    position: ValueBox(new Float32Array(vertices)),
                    normal: ValueBox(new Float32Array(normals)),
                    transform: ValueBox(transformArray)
                },
                {
                    colorTex,
                    colorTexSize: [ colorTexInfo.width, colorTexInfo.height ],
                    'light.position': Vec3.create(0, 0, -100),
                    'light.color': Vec3.create(1.0, 1.0, 1.0),
                    'light.ambient': Vec3.create(0.5, 0.5, 0.5),
                    'light.falloff': 0,
                    'light.radius': 500
                }
            )

            // console.log({ vertices, normals, vertexCount, atomCount })

            renderables.push(spheres)

            return true
        },
        update: (props: RepresentationProps) => false,
        draw: () => renderables.forEach(r => r.draw())
    }
}
