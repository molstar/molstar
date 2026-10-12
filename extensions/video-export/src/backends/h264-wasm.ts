/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import HME from 'h264-mp4-encoder';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { VideoRendererBackend } from '../backend.js';

/** Explicitly import this backend for environments without WebCodecs. */
export const H264WasmBackend: VideoRendererBackend = {
  name: 'H.264 WASM',
  dimensionAlignment: 2,
  async getParams() {
    return {
      quantization: PD.Numeric(
        18,
        { min: 10, max: 51, step: 1 },
        {
          description: 'Lower values preserve more detail.',
          help: () => ({
            description:
              'Controls compression for the H.264 WASM encoder, from 10 to 51. Lower values preserve more detail and generally produce larger files; higher values increase compression and reduce quality. The default is 18.',
          }),
        },
      ),
    };
  },
  getPresets(params) {
    const defaults = PD.getDefaultValues(params);
    return [
      [{ ...defaults, fps: 30, quantization: 18 }, 'Balanced (default)'],
      [{ ...defaults, fps: 15, quantization: 28 }, 'Fast preview'],
      [{ ...defaults, fps: 30, quantization: 10 }, 'High quality'],
      [{ ...defaults, fps: 60, quantization: 18 }, 'Smooth motion'],
    ];
  },
  async validate(config, values) {
    if (config.width % 2 || config.height % 2) return 'H.264 WASM requires even dimensions.';
    if (!Number.isInteger(values.quantization) || values.quantization < 10 || values.quantization > 51) {
      return 'Quantization must be an integer from 10 to 51.';
    }
  },
  async create(config, values) {
    const encoder = await HME.createH264MP4Encoder();
    try {
      encoder.width = config.width;
      encoder.height = config.height;
      encoder.frameRate = config.fps;
      encoder.quantizationParameter = values.quantization;
      encoder.initialize();
    } catch (e) {
      encoder.delete();
      throw e;
    }
    return {
      addFrame: (data) => encoder.addFrameRgba(data),
      finalize() {
        encoder.finalize();
        // Own the bytes before releasing the WASM heap.
        const data = new Uint8Array(encoder.FS.readFile(encoder.outputFilename));
        return { data, mimeType: 'video/mp4', extension: 'mp4' };
      },
      dispose: () => encoder.delete(),
    };
  },
};
