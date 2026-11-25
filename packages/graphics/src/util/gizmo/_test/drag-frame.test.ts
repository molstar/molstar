/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { ButtonsType, InputObserver, ModifiersKeys } from '@molstar/core/util/input/input-observer';
import { describe, expect, it } from 'vitest';
import { Camera } from '../../../canvas3d/camera.js';
import type { Canvas3D } from '../../../canvas3d/canvas3d.js';
import { TrackballControls } from '../../../canvas3d/controls/trackball.js';
import type { Scene } from '../../../scene/scene.js';
import { trackballTurnRate } from '../drag-frame.js';
import { type GizmoPose, turn } from '../gestures.js';

// a turn handle and the camera's trackball alike: a drag turns what is seen by as much, and the same way
describe('a turn at the trackball rate', () => {
  it('turns an item as the trackball turns the scene, on a canvas of any shape, cropped or not', () => {
    for (const [width, height, cropped] of [
      [800, 600, false],
      [300, 900, false],
      // a tile of an image twice the size: the field of view spans the whole image
      [800, 600, true],
    ] as const) {
      const camera = new Camera({
        position: Vec3.create(0, 0, 50),
        target: Vec3(),
        up: Vec3.unitY,
      });
      camera.setState({ radius: 10, radiusMax: 100 }, 0);
      camera.viewport = { x: 0, y: 0, width, height };
      if (cropped) {
        Camera.copyCrop(camera.crop, {
          enabled: true,
          fullWidth: 2 * width,
          fullHeight: 2 * height,
          offsetX: width / 2,
          offsetY: height / 2,
        });
      }
      const input = InputObserver.create();
      Object.defineProperties(input, {
        width: { value: width },
        height: { value: height },
      });
      const controls = TrackballControls.create(input, camera, { boundingSphereVisible: { radius: 10 } } as Scene, {});
      // as the canvas sets it: the frame (`Camera.getFrame`)
      Camera.getFrame(controls.viewport, camera);
      const drag = { dx: 0, dy: 0, pageX: 0, pageY: 0, x: 0, y: 0, useDelta: false };
      const base = { ...drag, buttons: ButtonsType.Flag.Primary, button: ButtonsType.Flag.Primary };
      const modifiers = ModifiersKeys.None;
      input.drag.next({ ...base, modifiers, isStart: true, pageX: 200, pageY: 200, x: 200, y: 200 } as never);
      input.drag.next({ ...base, modifiers, isStart: false, pageX: 230, pageY: 210, x: 230, y: 210 } as never);
      controls.update(16);

      // the scene's near side (+Z) as the camera now sees it: x right, y up, z toward the viewer
      const z = Vec3.normalize(Vec3(), camera.state.position);
      const x = Vec3.normalize(Vec3(), Vec3.cross(Vec3(), camera.state.up, z));
      const y = Vec3.cross(Vec3(), z, x);
      const seen = [x[2], y[2], z[2]];

      const rate = trackballTurnRate({ camera, input, props: { trackball: controls.props } } as unknown as Canvas3D);
      const pose: GizmoPose = { position: Vec3(), rotation: Quat.identity(), scale: Vec3.create(1, 1, 1) };
      const frame = { right: Vec3.create(1, 0, 0), up: Vec3.create(0, 1, 0), perPixel: 1, turnRate: rate };
      const turned = turn(pose, frame, 30, 10);
      const near = Vec3.transformQuat(Vec3(), Vec3.create(0, 0, 1), turned.rotation);
      for (let i = 0; i < 3; i++) expect(near[i]).toBeCloseTo(seen[i], 5);
      controls.dispose();
      input.dispose();
    }
  });
});
