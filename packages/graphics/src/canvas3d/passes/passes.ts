/**
 * Copyright (c) 2020-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { DrawPass } from './draw.js';
import { PickPass } from './pick.js';
import { MultiSamplePass } from './multi-sample.js';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { AssetManager } from '@molstar/core/util/assets';
import { IlluminationPass } from './illumination.js';
import { PassTargets } from './pass-targets.js';

export class Passes {
  /** the targets of the passes, see `PassTargets` */
  readonly targets: PassTargets;
  readonly draw: DrawPass;
  readonly pick: PickPass;
  readonly multiSample: MultiSamplePass;
  readonly illumination: IlluminationPass;

  constructor(
    private graphics: GraphicsContext,
    assetManager: AssetManager,
    attribs: Partial<{ pickScale: number; transparency: 'wboit' | 'dpoit' | 'blended' }> = {},
  ) {
    const drs = this.graphics.getDrawingBufferSize();
    this.targets = new PassTargets(graphics);
    this.draw = new DrawPass(this.targets, assetManager, drs.width, drs.height, attribs.transparency || 'blended');
    this.pick = new PickPass(graphics, drs.width, drs.height, attribs.pickScale || 0.25);
    this.multiSample = new MultiSamplePass(graphics, this.draw);
    this.illumination = new IlluminationPass(graphics, this.draw);
  }

  getByteCount() {
    return (
      this.draw.getByteCount() +
      this.pick.getByteCount() +
      this.multiSample.getByteCount() +
      this.illumination.getByteCount()
    );
  }

  setPickScale(pickScale: number) {
    this.pick.setPickScale(pickScale);
  }

  setTransparency(transparency: 'wboit' | 'dpoit' | 'blended') {
    this.draw.setTransparency(transparency);
  }

  updateSize() {
    const drs = this.graphics.getDrawingBufferSize();
    // Avoid setting dimensions to 0x0 because it causes "empty textures are not allowed" error.
    const width = Math.max(drs.width, 2);
    const height = Math.max(drs.height, 2);
    this.draw.setSize(width, height);
    this.pick.setSize(width, height);
    this.multiSample.syncSize();
    this.illumination.setSize(width, height);
  }
}
