/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { ObservableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';

/** The crop as fractions of the image (0–1), like Mol*'s `Viewport.relativeCrop`. */
export interface CropBox {
  x: number;
  y: number;
  width: number;
  height: number;
}

/**
 * Image export (implemented over `ViewportScreenshotHelper`): a live preview with the crop, the screenshot params
 * (resolution, format, transparency, axes, illumination), and copy or download.
 */
export interface ImageExportModel {
  /** A preview image URL (an object URL of a small render); a CSS background in the mocks. */
  preview: ReadonlyState<string | undefined>;
  /** Width over height of the image the settings produce. */
  aspect: ReadonlyState<number>;
  /** Crop to the visible content automatically, or keep the box the user dragged. */
  autoCrop: ObservableState<boolean>;
  crop: ReadonlyState<CropBox>;
  setCrop(crop: CropBox): void;
  /** The output size in pixels ("1920 × 1080"), for the line under the preview. */
  size: ReadonlyState<string>;
  form: ParamsFormModel<any>;
  busy: ReadonlyState<boolean>;
  copy(): void;
  download(): void;
}
