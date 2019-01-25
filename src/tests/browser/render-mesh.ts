/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html'
import { Canvas3D } from 'mol-canvas3d/canvas3d';
import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
import { Sphere } from 'mol-geo/primitive/sphere';
import { Mat4 } from 'mol-math/linear-algebra';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { Representation } from 'mol-repr/representation';
import { Color } from 'mol-util/color';
import { createRenderObject } from 'mol-gl/render-object';

const parent = document.getElementById('app')!
parent.style.width = '100%'
parent.style.height = '100%'

const canvas = document.createElement('canvas')
canvas.style.width = '100%'
canvas.style.height = '100%'
parent.appendChild(canvas)

const canvas3d = Canvas3D.create(canvas, parent)
canvas3d.animate()

function meshRepr() {
    const builderState = MeshBuilder.createState()
    const t = Mat4.identity()
    const sphere = Sphere(2)
    MeshBuilder.addPrimitive(builderState, t, sphere)
    const mesh = MeshBuilder.getMesh(builderState)

    const values = Mesh.Utils.createValuesSimple(mesh, {}, Color(0xFF0000), 1)
    const state = Mesh.Utils.createRenderableState({})
    const renderObject = createRenderObject('mesh', values, state)
    const repr = Representation.fromRenderObject('sphere-mesh', renderObject)
    return repr
}

canvas3d.add(meshRepr())
canvas3d.resetCamera()