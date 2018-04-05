/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { ValueBox } from 'mol-util/value-cell'
import * as glContext from 'mol-gl/context'
import { Camera } from 'mol-gl/camera'
import { Vec3, Mat4 } from 'mol-math/linear-algebra'
import { PointRenderable, MeshRenderable } from 'mol-gl/renderable'
import Model from 'mol-gl/model';
import { calculateTextureInfo } from 'mol-gl/util';
import Icosahedron from 'mol-geo/primitive/icosahedron'
import Box from 'mol-geo/primitive/box'
import Spacefill from 'mol-geo/representation/structure/spacefill'

import CIF from 'mol-io/reader/cif'
import { Run } from 'mol-task'
// import Computation from 'mol-util/computation'
import { AtomSet, Structure } from 'mol-model/structure'
import { UnitRepresentation } from 'mol-geo/representation/structure';

async function parseCif(data: string|Uint8Array) {
    const comp = CIF.parse(data)
    const parsed = await Run(comp);
    if (parsed.isError) throw parsed;
    return parsed
}

async function getPdb(pdb: string) {
    const data = await fetch(`https://files.rcsb.org/download/${pdb}.cif`)
    const parsed = await parseCif(await data.text())
    const structure = Structure.ofData({ kind: 'mmCIF', data: CIF.schema.mmCIF(parsed.result.blocks[0]) })
    return structure
}

import mcubes from './mcubes'

export default class State {
    regl: REGL.Regl

    async initRegl (container: HTMLDivElement) {
        const regl = glContext.create({
            container,
            extensions: [
                'OES_texture_float',
                'OES_texture_float_linear',
                'OES_element_index_uint',
                // 'EXT_disjoint_timer_query',
                'EXT_blend_minmax',
                'ANGLE_instanced_arrays'
            ],
            // profile: true
        })

        const camera = Camera.create(regl, container, {
            center: Vec3.create(0, 0, 0),
            near: 0.01,
            far: 1000
        })

        const p1 = Vec3.create(0, 4, 0)
        const p2 = Vec3.create(-3, 0, 0)

        const model1 = Model(regl)
        const model2 = Model(regl, { position: p1 })
        const model3 = Model(regl, { position: p2 })

        const position = ValueBox(new Float32Array([0, -1, 0, -1, 0, 0, 1, 1, 0]))
        const normal = ValueBox(new Float32Array([0, 0, 0, 0, 0, 0, 0, 0, 0]))

        const transformArray1 = ValueBox(new Float32Array(16))
        const transformArray2 = ValueBox(new Float32Array(16 * 3))
        const m4 = Mat4.identity()
        Mat4.toArray(m4, transformArray1.value, 0)
        Mat4.toArray(m4, transformArray2.value, 0)
        Mat4.setTranslation(m4, p1)
        Mat4.toArray(m4, transformArray2.value, 16)
        Mat4.setTranslation(m4, p2)
        Mat4.toArray(m4, transformArray2.value, 32)

        const colorTexInfo = calculateTextureInfo(3, 3)
        const color = new Uint8Array(colorTexInfo.length)
        color.set([
            0, 0, 255,
            0, 255, 0,
            255, 0, 0
        ])
        console.log(color, colorTexInfo)
        const colorTex = regl.texture({
            width: colorTexInfo.width,
            height: colorTexInfo.height,
            format: 'rgb',
            type: 'uint8',
            wrapS: 'clamp',
            wrapT: 'clamp',
            data: color
        })

        // position.update((array: Float32Array) => {
        //     positionFromModel({}, array, 0)
        // })

        const points = PointRenderable.create(regl, {
            position,
            transform: transformArray1
        })
        const mesh = MeshRenderable.create(regl,
            {
                position,
                normal,
                transform: transformArray2
            },
            {
                colorTex,
                colorTexSize: [ colorTexInfo.width, colorTexInfo.height ]
            }
        )

        const sphere = Icosahedron(1, 1)
        console.log(sphere)

        const box = Box(1, 1, 1, 1, 1, 1)
        console.log(box)

        const points2 = PointRenderable.create(regl, {
            position: ValueBox(new Float32Array(box.vertices)),
            transform: transformArray1
        })

        let rr = 1;
        function cubesF(x: number, y: number, z: number) {
            return x * x + y * y + z * z - rr * rr;
            // const a = ca;
            // const t = (x + y + z + a);
            // return x * x * x + y * y * y + z * z * z + a * a * a - t * t * t;
        }

        let cubes = await mcubes(cubesF);

        const makeCubesMesh = () => MeshRenderable.create(regl,
            {
                position: cubes.surface.vertexBuffer,
                normal: cubes.surface.normalBuffer as any,
                transform: transformArray1
            },
            {
                colorTex,
                colorTexSize: [ colorTexInfo.width, colorTexInfo.height ],
                'light.position': Vec3.create(0, 0, -100),
                'light.color': Vec3.create(1.0, 1.0, 1.0),
                'light.ambient': Vec3.create(0.5, 0.5, 0.5),
                'light.falloff': 0,
                'light.radius': 500
            },
            cubes.surface.indexBuffer.value
        );

        let mesh2 = makeCubesMesh();

        // const makeCubes = async () => {
        //     rr = Math.random();
        //     cubes = await mcubes(cubesF, cubes);
        //     mesh2 = makeCubesMesh();
        //     setTimeout(makeCubes, 1000 / 15);
        // };
        // makeCubes();

        // const mesh2 = MeshRenderable.create(regl,
        //     {
        //         position: Attribute.create(regl, new Float32Array(box.vertices), { size: 3 }),
        //         normal: Attribute.create(regl, new Float32Array(box.normals), { size: 3 }),
        //         ...createTransformAttributes(regl, transformArray1)
        //     },
        //     {
        //         colorTex,
        //         colorTexSize: [ colorTexInfo.width, colorTexInfo.height ],
        //         'light.position': Vec3.create(0, 0, -20),
        //         'light.color': Vec3.create(1.0, 1.0, 1.0),
        //         'light.ambient': Vec3.create(0.5, 0.5, 0.5),
        //         'light.falloff': 0,
        //         'light.radius': 500
        //     },
        //     box.indices
        // )

        function createSpacefills (structure: Structure) {
            const spacefills: UnitRepresentation[] = []
            const { atoms, units } = structure;
            const unitIds = AtomSet.unitIds(atoms);
            for (let i = 0, _i = unitIds.length; i < _i; i++) {
                const unitId = unitIds[i];
                const unit = units[unitId];
                const atomGroup = AtomSet.unitGetByIndex(atoms, i);

                const spacefill = Spacefill(regl)
                spacefill.create(unit, atomGroup, {})
                console.log('spacefill', spacefill)
                spacefills.push(spacefill)
            }
            return spacefills
        }

        const structures = await getPdb('1crn')
        const spacefills = createSpacefills(structures[0])

        structures[0]

        const baseContext = regl({
            context: {
                model: Mat4.identity(),
                transform: Mat4.setTranslation(Mat4.identity(), Vec3.create(6, 0, 0))
            },
            uniforms: {
                model: regl.context('model' as any),
                transform: regl.context('transform' as any),
            }
        })

        regl.frame((ctx) => {
            camera.update((state: any) => {
                if (!camera.isDirty()) return
                baseContext(() => {
                    // console.log(ctx)
                    regl.clear({color: [0, 0, 0, 1]})
                    spacefills.forEach(r => r.draw())
                    // position.update(array => { array[0] = Math.random() })
                    // points.update(a => { a.position[0] = Math.random() })
                    // mesh.draw()
                    // points.draw()
                    mesh2.draw()
                    // points2.draw()
                    // model1({}, ({ transform }) => {
                    //     points.draw()
                    // })
                    // model2({}, ({ transform }) => {
                    //     points.draw()
                    //     model3({ transform }, () => {
                    //         points.draw()
                    //     })
                    // })
                })
            }, undefined)
        })

        this.regl = regl
    }

    constructor() {

    }
}
