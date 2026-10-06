/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import './index.html';
import { Canvas3D, Canvas3DContext } from '@molstar/graphics/canvas3d/canvas3d';
import { TextBuilder } from '@molstar/graphics/geo/geometry/text/text-builder';
import { Text } from '@molstar/graphics/geo/geometry/text/text';
import { ParamDefinition, ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Color } from '@molstar/core/util/color';
import { Representation } from '@molstar/graphics/repr/representation';
import { SpheresBuilder } from '@molstar/graphics/geo/geometry/spheres/spheres-builder';
import { createRenderObject } from '@molstar/graphics/gl/render-object';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { resizeCanvas } from '@molstar/graphics/canvas3d/util';
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

  const values = Text.Utils.createValuesSimple(text, props, Color(0xffdd00), 1);
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

  const props = ParamDefinition.getDefaultValues(Spheres.Utils.Params);
  const values = Spheres.Utils.createValuesSimple(spheres, props, Color(0xff0000), 0.2);
  const state = Spheres.Utils.createRenderableState(props);
  const renderObject = createRenderObject('spheres', values, state, -1);
  console.log('spheres', renderObject);
  const repr = Representation.fromRenderObject('spheres', renderObject);
  return repr;
}

canvas3d.add(textRepr());
canvas3d.add(spheresRepr());
canvas3d.requestCameraReset();
