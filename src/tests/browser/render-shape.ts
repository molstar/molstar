/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html'
import { Canvas3D } from 'mol-canvas3d/canvas3d';
import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
import { Sphere } from 'mol-geo/primitive/sphere';
import { Mat4, Vec3 } from 'mol-math/linear-algebra';
import { Shape } from 'mol-model/shape';
import { ShapeRepresentation } from 'mol-repr/shape/representation';
import { ColorNames } from 'mol-util/color/tables';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { labelFirst } from 'mol-theme/label';

const parent = document.getElementById('app')!
parent.style.width = '100%'
parent.style.height = '100%'

const canvas = document.createElement('canvas')
canvas.style.width = '100%'
canvas.style.height = '100%'
parent.appendChild(canvas)

const info = document.createElement('div')
info.style.position = 'absolute'
info.style.fontFamily = 'sans-serif'
info.style.fontSize = '24pt'
info.style.bottom = '20px'
info.style.right = '20px'
info.style.color = 'white'
parent.appendChild(info)

const canvas3d = Canvas3D.create(canvas, parent)
canvas3d.animate()
canvas3d.input.move.subscribe(async ({x, y}) => {
    const pickingId = await canvas3d.identify(x, y)
    let label = ''
    if (pickingId) {
        const { loci } = canvas3d.getLoci(pickingId)
        label = labelFirst(loci)
    }
    info.innerText = label
})

const builderState = MeshBuilder.createState()
const t = Mat4.identity()
const sphere = Sphere(2)
builderState.currentGroup = 0
MeshBuilder.addPrimitive(builderState, t, sphere)
const mesh = MeshBuilder.getMesh(builderState)

const myData = {
    mesh,
    groupCount: 1,
    colors: [ColorNames.tomato, ColorNames.springgreen],
    labels: ['FooBaz0', 'FooBaz1'],
    transforms: [Mat4.identity(), Mat4.fromTranslation(Mat4.zero(), Vec3.create(3, 0, 0))]
}
type MyData = typeof myData
function getShape(data: MyData, props: {}, shape?: Shape<Mesh>) {
    const { mesh, colors, labels, transforms, groupCount } = data
    return shape || Shape.create(
        'test', mesh,
        (groupId: number, instanceId: number) => colors[instanceId * groupCount + groupId],
        (groupId: number, instanceId: number) => labels[instanceId * groupCount + groupId],
        transforms
    )
}

const repr = ShapeRepresentation(getShape, Mesh.Utils)

async function add() {
    await repr.createOrUpdate({}, myData).run()
    console.log(repr)
    canvas3d.add(repr)
    canvas3d.resetCamera()
}
add()

setTimeout(async () => {
    myData.colors[0] = ColorNames.darkmagenta
    await repr.createOrUpdate({}, myData).run()
}, 1000)