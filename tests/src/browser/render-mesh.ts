/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html';
import { resizeCanvas } from '@molstar/graphics/canvas3d/util';
import { Canvas3D, Canvas3DContext } from '@molstar/graphics/canvas3d/canvas3d';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import { HexagonalPrismCage } from '@molstar/graphics/geo/primitive/prism';
import { SpikedBall } from '@molstar/graphics/geo/primitive/spiked-ball';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { Color } from '@molstar/core/util/color';
import { createRenderObject } from '@molstar/graphics/gl/render-object';
import { Representation } from '@molstar/graphics/repr/representation';
import { Torus } from '@molstar/graphics/geo/primitive/torus';
import { ParamDefinition } from '@molstar/core/util/param-definition';
import { AssetManager } from '@molstar/core/util/assets';

const parent = document.getElementById('app')!;
parent.style.width = '100%';
parent.style.height = '100%';

const canvas = document.createElement('canvas');
parent.appendChild(canvas);

const assetManager = new AssetManager();

const canvas3dContext = Canvas3DContext.fromCanvas(canvas, assetManager);
const canvas3d = Canvas3D.create(canvas3dContext);
resizeCanvas(canvas, parent, canvas3dContext.pixelScale);
canvas3dContext.syncPixelScale();
canvas3d.requestResize();
canvas3d.animate();

canvas3d.input.resize.subscribe(() => {
    resizeCanvas(canvas, parent, canvas3dContext.pixelScale);
    canvas3dContext.syncPixelScale();
    canvas3d.requestResize();
});

function meshRepr() {
    const builderState = MeshBuilder.createState();

    const t = Mat4.identity();
    Mat4.scaleUniformly(t, t, 10);
    MeshBuilder.addCage(builderState, t, HexagonalPrismCage(), 0.05, 2, 20);

    const t2 = Mat4.identity();
    Mat4.scaleUniformly(t2, t2, 1);
    MeshBuilder.addPrimitive(builderState, t2, SpikedBall(3));

    const t3 = Mat4.identity();
    Mat4.scaleUniformly(t3, t3, 8);
    MeshBuilder.addPrimitive(builderState, t3, Torus({ tubularSegments: 64, radialSegments: 32, tube: 0.1 }));

    const mesh = MeshBuilder.getMesh(builderState);

    const props = ParamDefinition.getDefaultValues(Mesh.Utils.Params);
    const values = Mesh.Utils.createValuesSimple(mesh, props, Color(0xFF4433), 1);
    const state = Mesh.Utils.createRenderableState(props);
    const renderObject = createRenderObject('mesh', values, state, -1);
    console.log('mesh', renderObject);
    const repr = Representation.fromRenderObject('mesh', renderObject);
    return repr;
}

canvas3d.add(meshRepr());
canvas3d.requestCameraReset();