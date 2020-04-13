/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html';
import { resizeCanvas } from '../../mol-canvas3d/util';
import { Canvas3D } from '../../mol-canvas3d/canvas3d';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { Mat4 } from '../../mol-math/linear-algebra';
import { HexagonalPrismCage } from '../../mol-geo/primitive/prism';
import { SpikedBall } from '../../mol-geo/primitive/spiked-ball';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
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

function meshRepr() {
    const builderState = MeshBuilder.createState();

    const t = Mat4.identity();
    MeshBuilder.addCage(builderState, t, HexagonalPrismCage(), 0.005, 2, 20);

    const t2 = Mat4.identity();
    Mat4.scaleUniformly(t2, t2, 0.1);
    MeshBuilder.addPrimitive(builderState, t2, SpikedBall(3));

    const mesh = MeshBuilder.getMesh(builderState);

    const values = Mesh.Utils.createValuesSimple(mesh, {}, Color(0xFF0000), 1);
    const state = Mesh.Utils.createRenderableState({});
    const renderObject = createRenderObject('mesh', values, state, -1);
    const repr = Representation.fromRenderObject('mesh', renderObject);
    return repr;
}

canvas3d.add(meshRepr());
canvas3d.requestCameraReset();