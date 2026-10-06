/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html';
import { resizeCanvas } from '@molstar/graphics/canvas3d/util';
import { Canvas3D, Canvas3DContext } from '@molstar/graphics/canvas3d/canvas3d';
import { SpheresBuilder } from '@molstar/graphics/geo/geometry/spheres/spheres-builder';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { Color } from '@molstar/core/util/color';
import { createRenderObject } from '@molstar/graphics/gl/render-object';
import { Representation } from '@molstar/graphics/repr/representation';
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

function spheresRepr() {
  const spheresBuilder = SpheresBuilder.create(3, 1);
  spheresBuilder.add(0, 0, 0, 0);
  spheresBuilder.add(5, 0, 0, 0);
  spheresBuilder.add(-4, 1, 0, 0);
  const spheres = spheresBuilder.getSpheres();

  const props = ParamDefinition.getDefaultValues(Spheres.Utils.Params);
  const values = Spheres.Utils.createValuesSimple(spheres, {}, Color(0xff0000), 1);
  const state = Spheres.Utils.createRenderableState(props);
  const renderObject = createRenderObject('spheres', values, state, -1);
  console.log(renderObject);
  const repr = Representation.fromRenderObject('spheres', renderObject);
  return repr;
}

canvas3d.add(spheresRepr());
canvas3d.requestCameraReset();
