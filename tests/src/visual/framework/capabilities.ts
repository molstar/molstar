/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { DpoitPass } from '@molstar/graphics/canvas3d/passes/dpoit';
import { IlluminationPass } from '@molstar/graphics/canvas3d/passes/illumination';
import { WboitPass } from '@molstar/graphics/canvas3d/passes/wboit';
import { getGaussianSurfaceGpuUnsupportedReason } from '@molstar/graphics/repr/structure/visual/gaussian-surface-mesh';
import { GraphicsFeatures } from '@molstar/graphics/gl/capabilities';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { isWebGLContext, type WebGLContext } from '@molstar/graphics/gl/webgl/context';

/**
 * Backend-neutral capability keys: the graphics features plus derived rendering features.
 * Features required by Mol* (e.g. float textures, min/max blending) are not listed.
 */
export const CapabilityKeys = [
  ...GraphicsFeatures,
  'wboit',
  'dpoit',
  'illumination',
  'gaussianSurfaceGpu',
  'releaseContext',
] as const;
export type CapabilityKey = (typeof CapabilityKeys)[number];

export type Capabilities = { [K in CapabilityKey]: boolean };

export interface CapabilityReport {
  backend: string;
  version: string;
  renderer: string;
  vendor: string;
  userAgent: string;
  capabilities: Capabilities;
  limits: Record<string, number>;
  extensions: string[];
}

/** Why derived capabilities are missing, per capabilities object of `getCapabilities` */
const MissingReasons = new WeakMap<Capabilities, { [K in CapabilityKey]?: string }>();

export function getCapabilities(graphics: GraphicsContext): Capabilities {
  const capabilities = {} as Capabilities;
  for (const k of GraphicsFeatures) capabilities[k] = graphics.capabilities[k];
  const reasons: { [K in CapabilityKey]?: string } = {
    wboit: WboitPass.getUnsupportedReason(graphics),
    dpoit: DpoitPass.getUnsupportedReason(graphics),
    illumination: IlluminationPass.getUnsupportedReason(graphics),
    gaussianSurfaceGpu: getGaussianSurfaceGpuUnsupportedReason(graphics),
    // probed without releasing: WebGPU destroys its device, WebGL needs the extension
    releaseContext:
      isWebGLContext(graphics) && !graphics.gl.getExtension('WEBGL_lose_context') ? 'no WEBGL_lose_context' : undefined,
  };
  capabilities.wboit = !reasons.wboit;
  capabilities.dpoit = !reasons.dpoit;
  capabilities.illumination = !reasons.illumination;
  capabilities.gaussianSurfaceGpu = !reasons.gaussianSurfaceGpu;
  capabilities.releaseContext = !reasons.releaseContext;
  MissingReasons.set(capabilities, reasons);
  return capabilities;
}

/** `keys` with the reasons of derived capabilities, e.g. `dpoit (missing floatBlend)` */
export function describeMissing(capabilities: Capabilities, keys: readonly CapabilityKey[]) {
  const reasons = MissingReasons.get(capabilities);
  return keys.map((k) => (reasons?.[k] ? `${k} (${reasons[k]})` : k));
}

export function getWebGLCapabilityReport(webgl: WebGLContext): CapabilityReport {
  const { gl } = webgl;
  const debugInfo = gl.getExtension('WEBGL_debug_renderer_info');
  const renderer = debugInfo ? gl.getParameter(debugInfo.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER);
  const vendor = debugInfo ? gl.getParameter(debugInfo.UNMASKED_VENDOR_WEBGL) : gl.getParameter(gl.VENDOR);
  return {
    backend: webgl.capabilities.backend,
    version: `${gl.getParameter(gl.VERSION)}`,
    renderer: `${renderer}`,
    vendor: `${vendor}`,
    userAgent: navigator.userAgent,
    capabilities: getCapabilities(webgl),
    limits: { ...webgl.resources.limits },
    extensions: gl.getSupportedExtensions() ?? [],
  };
}

export function getMissingCapabilities(
  capabilities: Capabilities,
  required: readonly CapabilityKey[],
): CapabilityKey[] {
  return required.filter((k) => !capabilities[k]);
}
