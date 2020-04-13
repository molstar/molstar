/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html';
import { Canvas3D } from '../../mol-canvas3d/canvas3d';
import { TextBuilder } from '../../mol-geo/geometry/text/text-builder';
import { Text } from '../../mol-geo/geometry/text/text';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Color } from '../../mol-util/color';
import { Representation } from '../../mol-repr/representation';
import { SpheresBuilder } from '../../mol-geo/geometry/spheres/spheres-builder';
import { createRenderObject } from '../../mol-gl/render-object';
import { Spheres } from '../../mol-geo/geometry/spheres/spheres';
import { resizeCanvas } from '../../mol-canvas3d/util';

const parent = document.getElementById('app')!;
parent.style.width = '100%';
parent.style.height = '100%';

const canvas = document.createElement('canvas');
parent.appendChild(canvas);
resizeCanvas(canvas, parent);

const canvas3d = Canvas3D.fromCanvas(canvas);
canvas3d.animate();

function textRepr() {
    const props: PD.Values<Text.Params> = {
        ...PD.getDefaultValues(Text.Params),
        attachment: 'top-right',
        fontQuality: 3,
        fontWeight: 'normal',
        borderWidth: 0.3,
        background: true,
        backgroundOpacity: 0.5,
        tether: true,
        tetherLength: 1.5,
        tetherBaseWidth: 0.5,
    };

    const textBuilder = TextBuilder.create(props, 1, 1);
    textBuilder.add('Hello world', 0, 0, 0, 1, 1, 0);
    // textBuilder.add('Добрый день', 0, 1, 0, 0, 0)
    // textBuilder.add('美好的一天', 0, 2, 0, 0, 0)
    // textBuilder.add('¿Cómo estás?', 0, -1, 0, 0, 0)
    // textBuilder.add('αβγ Å', 0, -2, 0, 0, 0)
    const text = textBuilder.getText();

    const values = Text.Utils.createValuesSimple(text, props, Color(0xFFDD00), 1);
    const state = Text.Utils.createRenderableState(props);
    const renderObject = createRenderObject('text', values, state, -1);
    console.log('text', renderObject, props);
    const repr = Representation.fromRenderObject('text', renderObject);
    return repr;
}

function spheresRepr() {
    const spheresBuilder = SpheresBuilder.create(1, 1);
    spheresBuilder.add(0, 0, 0, 0);
    spheresBuilder.add(5, 0, 0, 0);
    spheresBuilder.add(-4, 1, 0, 0);
    const spheres = spheresBuilder.getSpheres();

    const values = Spheres.Utils.createValuesSimple(spheres, {}, Color(0xFF0000), 0.2);
    const state = Spheres.Utils.createRenderableState({});
    const renderObject = createRenderObject('spheres', values, state, -1);
    console.log('spheres', renderObject);
    const repr = Representation.fromRenderObject('spheres', renderObject);
    return repr;
}

canvas3d.add(textRepr());
canvas3d.add(spheresRepr());
canvas3d.requestCameraReset();