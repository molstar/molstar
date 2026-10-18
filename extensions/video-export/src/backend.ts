/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

export interface VideoRendererConfig {
  width: number;
  height: number;
  fps: number;
}

export interface VideoResult {
  data: Uint8Array<ArrayBuffer>;
  mimeType: string;
  extension: string;
}

export interface VideoRendererSession {
  addFrame(data: Uint8Array | Uint8ClampedArray, timestamp: number, duration: number): Promise<void> | void;
  finalize(): Promise<VideoResult> | VideoResult;
  dispose(): Promise<void> | void;
}

/** Encoding only; animation playback and frame capture are shared by all backends. */
export interface VideoRendererBackend {
  readonly name: string;
  readonly dimensionAlignment: number;
  getParams(config: VideoRendererConfig, values?: Record<string, any>): Promise<PD.Params>;
  /** Optional presets for the combined frame-rate and backend parameter schema. */
  getPresets?(params: PD.Params, values?: Record<string, any>): PD.Select<Record<string, any>>['options'];
  /** Undefined when the configuration can be encoded, otherwise a user-facing reason. */
  validate(config: VideoRendererConfig, values: Record<string, any>): Promise<string | undefined>;
  create(config: VideoRendererConfig, values: Record<string, any>): Promise<VideoRendererSession>;
}

/** Reject supplied choices before PD normalization can replace them with defaults. */
function validateChoices(params: PD.Params, values: Record<string, any>, prefix = '') {
  for (const [key, value] of Object.entries(values)) {
    if (value === undefined) continue;
    const path = prefix ? `${prefix}.${key}` : key;
    const param = params[key];
    // FPS is shared by the renderer and may accompany backend settings from the UI/presets.
    if (!prefix && key === 'fps' && !param) continue;
    if (!param) throw new Error(`Unknown video parameter: ${path}`);
    if (param.type === 'select') {
      if (!param.options.some(([option]) => option === value)) throw new Error(`Unsupported video ${path}: ${value}`);
    } else if (param.type === 'group') {
      if (!value || typeof value !== 'object' || Array.isArray(value)) throw new Error(`Invalid video ${path}`);
      validateChoices(param.params, value, path);
    } else if (param.type === 'mapped') {
      if (!value || !param.select.options.some(([option]) => option === value.name)) {
        throw new Error(`Unsupported video ${path}: ${value?.name}`);
      }
      const child = param.map(value.name);
      if (child.type === 'group' && value.params !== undefined) {
        if (!value.params || typeof value.params !== 'object' || Array.isArray(value.params)) {
          throw new Error(`Invalid video ${path}.params`);
        }
        validateChoices(child.params, value.params, `${path}.params`);
      }
    }
  }
}

/** Fill defaults while rejecting explicitly supplied settings which cannot be used. */
export async function resolveVideoParams(backend: VideoRendererBackend, config: VideoRendererConfig, values = {}) {
  const params = await backend.getParams(config, values);
  validateChoices(params, values);
  const normalized = PD.normalizeParams(params, values, 'all');
  for (const key of ['container', 'codec']) {
    const requested = (values as Record<string, any>)[key];
    if (requested !== undefined && (requested?.name ?? requested) !== (normalized[key]?.name ?? normalized[key])) {
      throw new Error(`Unsupported video ${key}: ${requested?.name ?? requested}`);
    }
  }
  const reason = await backend.validate(config, normalized);
  if (reason) throw new Error(reason);
  return normalized;
}
