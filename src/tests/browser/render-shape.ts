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
import { RuntimeContext, Progress } from 'mol-task';
import { Representation } from 'mol-repr/representation';
import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { EveryLoci } from 'mol-model/loci';




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

let prevReprLoci = Representation.Loci.Empty
const canvas3d = Canvas3D.create(canvas, parent)
canvas3d.animate()
canvas3d.input.move.subscribe(async ({x, y}) => {
    const pickingId = await canvas3d.identify(x, y)
    let label = ''
    if (pickingId) {
        const reprLoci = canvas3d.getLoci(pickingId)
        label = labelFirst(reprLoci.loci)
        if (!Representation.Loci.areEqual(prevReprLoci, reprLoci)) {
            canvas3d.mark(prevReprLoci, MarkerAction.RemoveHighlight)
            canvas3d.mark(reprLoci, MarkerAction.Highlight)
            prevReprLoci = reprLoci
        }
    } else {
        canvas3d.mark({ loci: EveryLoci }, MarkerAction.RemoveHighlight)
        prevReprLoci = Representation.Loci.Empty
    }
    info.innerText = label
})

/**
 * Create a mesh of spheres at given centers
 * - asynchronous (using async/await)
 * - progress tracking (via `ctx.update`)
 * - re-use storage from an existing mesh if given
 */
async function getSphereMesh(ctx: RuntimeContext, centers: number[], mesh?: Mesh) {
    const builderState = MeshBuilder.createState(centers.length * 128, centers.length * 128 / 2, mesh)
    const t = Mat4.identity()
    const v = Vec3.zero()
    const sphere = Sphere(4)
    builderState.currentGroup = 0
    for (let i = 0, il = centers.length / 3; i < il; ++i) {
        // for production, calls to update should be guarded by `if (ctx.shouldUpdate)`
        await ctx.update({ current: i, max: il, message: `adding sphere ${i}` })
        builderState.currentGroup = i
        Mat4.setTranslation(t, Vec3.fromArray(v, centers, i * 3))
        MeshBuilder.addPrimitive(builderState, t, sphere)
    }
    return MeshBuilder.getMesh(builderState)
}

const myData = {
    centers: [0, 0, 0, 0, 3, 0, 1, 0 , 4],
    colors: [ColorNames.tomato, ColorNames.springgreen,ColorNames.springgreen],
    labels: ['Sphere 0, Instance A', 'Sphere 1, Instance A', 'Sphere 0, Instance B', 'Sphere 1, Instance B'],
    transforms: [Mat4.identity(), Mat4.fromTranslation(Mat4.zero(), Vec3.create(3, 0, 0))]
}
type MyData = typeof myData

/**
 * Get shape from `MyData` object
 */
async function getShape(ctx: RuntimeContext, data: MyData, props: {}, shape?: Shape<Mesh>) {
    await ctx.update('async creation of shape from  myData')
    const { centers, colors, labels, transforms } = data
    const mesh = await getSphereMesh(ctx, centers, shape && shape.geometry)
    const groupCount = centers.length / 3
    return shape || Shape.create(
        'test', mesh,
        (groupId: number) => colors[groupId], // color: per group, same for instances
        () => 1, // size: constant
        (groupId: number, instanceId: number) => labels[instanceId * groupCount + groupId], // label: per group and instance
        transforms
    )
}

// Init ShapeRepresentation container
const repr = ShapeRepresentation(getShape, Mesh.Utils)

export async function init() {
    // Create shape from myData and add to canvas3d
    await repr.createOrUpdate({}, myData).run((p: Progress) => console.log(Progress.format(p)))
    canvas3d.add(repr)
    canvas3d.resetCamera()

    // Change color after 1s
    setTimeout(async () => {
        myData.colors[0] = ColorNames.darkmagenta
        // Calling `createOrUpdate` with `data` will trigger color and transform update
        await repr.createOrUpdate({}, myData).run()
    }, 1000)
}
export default init();