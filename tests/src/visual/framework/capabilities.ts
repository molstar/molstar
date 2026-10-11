/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { DpoitPass } from '@molstar/graphics/canvas3d/passes/dpoit';
import { IlluminationPass } from '@molstar/graphics/canvas3d/passes/illumination';
import { WboitPass } from '@molstar/graphics/canvas3d/passes/wboit';
import { GraphicsFeatures } from '@molstar/graphics/gl/capabilities';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { type WebGLContext } from '@molstar/graphics/gl/webgl/context';

/**
 * Backend-neutral capability keys: the graphics features plus derived rendering features.
 * Features required by Mol* (e.g. float textures, min/max blending) are not listed.
 */
export const CapabilityKeys = [...GraphicsFeatures, 'wboit', 'dpoit', 'illumination'] as const;
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

export function getCapabilities(graphics: GraphicsContext): Capabilities {
  const capabilities = {} as Capabilities;
  for (const k of GraphicsFeatures) capabilities[k] = graphics.capabilities[k];
  capabilities.wboit = WboitPass.isSupported(graphics);
  capabilities.dpoit = DpoitPass.isSupported(graphics);
  capabilities.illumination = IlluminationPass.isSupported(graphics);
  return capabilities;
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
