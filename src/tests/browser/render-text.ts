/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html'
import { Canvas3D } from 'mol-canvas3d/canvas3d';
import { Geometry } from 'mol-geo/geometry/geometry';
import { TextBuilder } from 'mol-geo/geometry/text/text-builder';
import { Text } from 'mol-geo/geometry/text/text';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Color } from 'mol-util/color';
import { createTextRenderObject, createSpheresRenderObject } from 'mol-gl/render-object';
import { Representation } from 'mol-repr/representation';
import { SpheresBuilder } from 'mol-geo/geometry/spheres/spheres-builder';
import { Spheres } from 'mol-geo/geometry/spheres/spheres';

const parent = document.getElementById('app')!
parent.style.width = '100%'
parent.style.height = '100%'

const canvas = document.createElement('canvas')
canvas.style.width = '100%'
canvas.style.height = '100%'
parent.appendChild(canvas)

const canvas3d = Canvas3D.create(canvas, parent)
canvas3d.animate()

function textRepr() {
    const props: PD.Values<Text.Params> = {
        ...PD.getDefaultValues(Text.Params),
        attachment: 'middle-center',
        fontSize: 96,
        fontWeight: 'bold',
    }

    const textBuilder = TextBuilder.create(props, 1, 1)
    textBuilder.add('Hello world', 0, 0, 0, 0)
    textBuilder.add('Добрый день', 0, 1, 0, 0)
    textBuilder.add('美好的一天', 0, 2, 0, 0)
    textBuilder.add('¿Cómo estás?', 0, -1, 0, 0)
    textBuilder.add('αβγ Å', 0, -2, 0, 0)
    const text = textBuilder.getText()

    const values = Text.createValuesSimple(text, props, Color(0xFFDD00), 1)
    const state = Text.createRenderableState(props)
    const renderObject = createTextRenderObject(values, state)
    console.log('text', renderObject)
    const repr = Representation.fromRenderObject('text', renderObject)
    return repr
}

function spheresRepr() {
    const spheresBuilder = SpheresBuilder.create(2, 1)
    spheresBuilder.add(5, 0, 0, 0)
    spheresBuilder.add(-4, 1, 0, 0)
    const spheres = spheresBuilder.getSpheres()

    const values = Spheres.createValuesSimple(spheres, {}, Color(0xFF0000), 1)
    const state = Geometry.createRenderableState()
    const renderObject = createSpheresRenderObject(values, state)
    console.log('spheres', renderObject)
    const repr = Representation.fromRenderObject('spheres', renderObject)
    return repr
}

canvas3d.add(textRepr())
canvas3d.add(spheresRepr())
canvas3d.resetCamera()