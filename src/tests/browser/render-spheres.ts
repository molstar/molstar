/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html'
import { Canvas3D } from 'mol-canvas3d/canvas3d';
import { SpheresBuilder } from 'mol-geo/geometry/spheres/spheres-builder';
import { Geometry } from 'mol-geo/geometry/geometry';
import { createSpheresRenderObject } from 'mol-gl/render-object';
import { Representation } from 'mol-repr/representation';
import { Spheres } from 'mol-geo/geometry/spheres/spheres';
import { Color } from 'mol-util/color';

const parent = document.getElementById('app')!
parent.style.width = '100%'
parent.style.height = '100%'

const canvas = document.createElement('canvas')
canvas.style.width = '100%'
canvas.style.height = '100%'
parent.appendChild(canvas)

const canvas3d = Canvas3D.create(canvas, parent)
canvas3d.animate()

function spheresRepr() {
    const spheresBuilder = SpheresBuilder.create(3, 1)
    spheresBuilder.add(0, 0, 0, 0)
    spheresBuilder.add(5, 0, 0, 0)
    spheresBuilder.add(-4, 1, 0, 0)
    const spheres = spheresBuilder.getSpheres()

    const values = Spheres.createValuesSimple(spheres, {}, Color(0xFF0000), 1)
    const state = Geometry.createRenderableState()
    const renderObject = createSpheresRenderObject(values, state)
    console.log(renderObject)
    const repr = Representation.fromRenderObject('spheres', renderObject)
    return repr
}

canvas3d.add(spheresRepr())
canvas3d.resetCamera()