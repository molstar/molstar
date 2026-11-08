/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Scene } from '@molstar/graphics/gl/scene';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { DebugRegistry } from './debug-registry.js';
import { CameraAxes, CameraHelperParams } from './camera-helper.js';
import { OverlayHelper } from './overlay-helper.js';
import { PointerHelper, PointerHelperParams } from './pointer-helper.js';
import { ScaleLine } from './scale-line.js';
import { ScreenHelper } from './screen-helper.js';

export const HelperParams = {
  camera: PD.Group({
    helper: PD.Group(CameraHelperParams),
  }),
  pointer: PD.Group(PointerHelperParams),
};
export const DefaultHelperProps = PD.getDefaultValues(HelperParams);
export type HelperProps = PD.Values<typeof HelperParams>;

export class Helper {
  readonly debug: DebugRegistry;
  /** the overlay layer, see `OverlayHelper` */
  readonly overlay: OverlayHelper;
  /** the screen layer, see `ScreenHelper` */
  readonly screen: ScreenHelper;
  /** the corner axes, in the screen layer */
  readonly axes: CameraAxes;
  /** the scale line, in the screen layer */
  readonly scaleLine: ScaleLine;
  readonly pointer: PointerHelper;

  constructor(graphics: GraphicsContext, scene: Scene, props: Partial<HelperProps> = {}) {
    const p = { ...DefaultHelperProps, ...props };

    this.debug = new DebugRegistry(graphics, scene);
    this.overlay = new OverlayHelper(graphics);
    this.screen = new ScreenHelper(graphics);
    this.axes = new CameraAxes(this.screen, p.camera.helper);
    this.scaleLine = new ScaleLine(this.screen);
    this.scaleLine.setProps(p.camera.helper.scale);
    this.pointer = new PointerHelper(graphics, p.pointer);
  }
}
