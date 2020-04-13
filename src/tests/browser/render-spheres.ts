/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html';
import { resizeCanvas } from '../../mol-canvas3d/util';
import { Canvas3D } from '../../mol-canvas3d/canvas3d';
import { SpheresBuilder } from '../../mol-geo/geometry/spheres/spheres-builder';
import { Spheres } from '../../mol-geo/geometry/spheres/spheres';
import { Color } from '../../mol-util/color';
import { createRenderObject } from '../../mol-gl/render-object';
import { Representation } from '../../mol-repr/representation';

const parent = document.getElementById('app')!;
parent.style.width = '100%';
parent.style.height = '100%';

const canvas = document.createElement('canvas');
parent.appendChild(canvas);
resizeCanvas(canvas, parent);

const canvas3d = Canvas3D.fromCanvas(canvas);
canvas3d.animate();

function spheresRepr() {
    const spheresBuilder = SpheresBuilder.create(3, 1);
    spheresBuilder.add(0, 0, 0, 0);
    spheresBuilder.add(5, 0, 0, 0);
    spheresBuilder.add(-4, 1, 0, 0);
    const spheres = spheresBuilder.getSpheres();

    const values = Spheres.Utils.createValuesSimple(spheres, {}, Color(0xFF0000), 1);
    const state = Spheres.Utils.createRenderableState({});
    const renderObject = createRenderObject('spheres', values, state, -1);
    console.log(renderObject);
    const repr = Representation.fromRenderObject('spheres', renderObject);
    return repr;
}

canvas3d.add(spheresRepr());
canvas3d.requestCameraReset();