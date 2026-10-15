/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Canvas3DParams, type Canvas3DProps, type PartialCanvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import { type ApiContext, ApiResult, type OperationDefinition } from '../../api/area.js';

// The `canvas-settings` feature's operations (operations.md §2, the `canvas` row; viewer-spec.md §3.1): added to the
// core's `canvas` area, as only built-in features may. The card's sections edit through `canvas.settings` (one section's
// values, as the classic simple settings map them); tools and commands take the named ones.

/** The card's sections, as the classic simple settings group Canvas3D's params. */
export const CanvasSections = {
  background: PD.Group({
    color: PD.Color(Canvas3DParams.renderer.params.backgroundColor.defaultValue, { label: 'Color' }),
    transparent: PD.Boolean(false),
    style: Canvas3DParams.postprocessing.params.background,
  }),
  lighting: PD.Group({
    occlusion: Canvas3DParams.postprocessing.params.occlusion,
    shadow: Canvas3DParams.postprocessing.params.shadow,
    outline: Canvas3DParams.postprocessing.params.outline,
    dof: { ...Canvas3DParams.postprocessing.params.dof, label: 'Depth of field' },
    fog: Canvas3DParams.cameraFog,
  }),
  camera: Canvas3DParams.camera,
  clipping: PD.Group({ ...Canvas3DParams.cameraClipping.params }),
  advanced: PD.Group({
    illumination: Canvas3DParams.illumination,
    multiSample: Canvas3DParams.multiSample,
    sharpening: Canvas3DParams.postprocessing.params.sharpening,
    bloom: Canvas3DParams.postprocessing.params.bloom,
  }),
};

export type CanvasSection = keyof typeof CanvasSections;

/** A section's values, read from the canvas's props. */
export function sectionValues(section: CanvasSection, p: Canvas3DProps): Record<string, unknown> {
  switch (section) {
    case 'background':
      return {
        color: p.renderer.backgroundColor,
        transparent: p.transparentBackground,
        style: p.postprocessing.background,
      };
    case 'lighting':
      return {
        occlusion: p.postprocessing.occlusion,
        shadow: p.postprocessing.shadow,
        outline: p.postprocessing.outline,
        dof: p.postprocessing.dof,
        fog: p.cameraFog,
      };
    case 'camera':
      return p.camera as unknown as Record<string, unknown>;
    case 'clipping':
      return p.cameraClipping as unknown as Record<string, unknown>;
    case 'advanced':
      return {
        illumination: p.illumination,
        multiSample: p.multiSample,
        sharpening: p.postprocessing.sharpening,
        bloom: p.postprocessing.bloom,
      };
  }
}

/** The canvas props a section's values set. */
export function sectionProps(section: CanvasSection, v: Record<string, any>): PartialCanvas3DProps {
  switch (section) {
    case 'background':
      return {
        renderer: { backgroundColor: v.color },
        transparentBackground: v.transparent,
        postprocessing: { background: v.style },
      };
    case 'lighting':
      return {
        postprocessing: { occlusion: v.occlusion, shadow: v.shadow, outline: v.outline, dof: v.dof },
        cameraFog: v.fog,
      };
    case 'camera':
      return { camera: v as never };
    case 'clipping':
      return { cameraClipping: v as never };
    case 'advanced':
      return {
        illumination: v.illumination,
        multiSample: v.multiSample,
        postprocessing: { sharpening: v.sharpening, bloom: v.bloom },
      };
  }
}

const OnOff = () => PD.Select<'keep' | 'on' | 'off'>('keep', PD.arrayToOptions(['keep', 'on', 'off'] as const));

export const SettingsParams = {
  section: PD.Select<CanvasSection>('background', PD.arrayToOptions(Object.keys(CanvasSections) as CanvasSection[])),
  values: PD.Value<Record<string, unknown>>({}),
};
export const LightingParams = { occlusion: OnOff(), shadow: OnOff(), outline: OnOff(), fog: OnOff() };
export const IlluminationParams = { on: PD.Boolean(true, { description: 'Path-traced global illumination' }) };
export const CameraParams = {
  mode: PD.Select<'keep' | 'perspective' | 'orthographic'>(
    'keep',
    PD.arrayToOptions(['keep', 'perspective', 'orthographic'] as const),
  ),
  fov: PD.Numeric(0, { min: 0, max: 130, step: 1 }, { description: 'Field of view in degrees; 0: keep' }),
};
export const ClippingParams = {
  radius: PD.Numeric(100, { min: 0, max: 100, step: 1 }, { description: 'How much of the scene to show (100: all)' }),
  far: PD.Boolean(true, { description: 'Hide the scene in the distance' }),
};

const canvasModule = ['canvas'] as const;

function props(ctx: ApiContext): Canvas3DProps | undefined {
  return ctx.plugin.canvas3d?.props;
}

function set(ctx: ApiContext, change: PartialCanvas3DProps, did: string) {
  ctx.plugin.setCanvas3DProps(change);
  return ApiResult.ok(did, { view: true });
}

/** A mapped param's 'on' case's defaults. */
const defaultsOf = (m: PD.Any) => (m.type === 'group' ? PD.getDefaultValues(m.params) : {});

/** On or off, as Canvas3D's mapped params say it; turning on keeps the current settings, or takes the defaults. */
const onOff = (on: boolean, current: { name: string; params: unknown }, defaults: { name: string; params: unknown }) =>
  on ? (current.name === 'on' ? current : defaults) : { name: 'off', params: {} };

/** Checks a definition against its own params, then lists it with the others. */
const typed = <P extends PD.Params>(d: OperationDefinition<P, unknown>) => d as OperationDefinition<any, any>;

export const canvasSettingsOperations = (): {
  area: string;
  name: string;
  operation: OperationDefinition<any, any>;
}[] => [
  {
    area: 'canvas',
    name: 'settings',
    operation: typed({
      id: 'canvas.settings',
      title: 'Change canvas settings',
      params: SettingsParams,
      modules: canvasModule,
      surface: {
        kind: 'ui-only',
        reason: "The Canvas card's sections; tools use lighting, illumination, camera, clipping",
      },
      async run({ section, values }, ctx) {
        const p = props(ctx);
        if (!p) return ApiResult.fail('unavailable', 'There is no canvas.');
        const current = sectionValues(section, p);
        const next = { ...current, ...values };
        if (PD.areEqual(CanvasSections[section].params, current, next)) {
          return ApiResult.already(`The ${section} settings are already so.`);
        }
        return set(ctx, sectionProps(section, next), `Changed the ${section} settings.`);
      },
    }),
  },
  {
    area: 'canvas',
    name: 'lighting',
    operation: typed({
      id: 'canvas.lighting',
      title: 'Change the lighting',
      params: LightingParams,
      modules: canvasModule,
      async run(v, ctx) {
        const p = props(ctx);
        if (!p) return ApiResult.fail('unavailable', 'There is no canvas.');
        const pp = Canvas3DParams.postprocessing.params as unknown as Record<string, PD.Mapped<any>>;
        const wanted = (key: 'occlusion' | 'shadow' | 'outline') =>
          v[key] === 'keep'
            ? p.postprocessing[key]
            : onOff(v[key] === 'on', p.postprocessing[key], { name: 'on', params: defaultsOf(pp[key].map('on')) });
        const change = {
          postprocessing: { occlusion: wanted('occlusion'), shadow: wanted('shadow'), outline: wanted('outline') },
          cameraFog:
            v.fog === 'keep'
              ? p.cameraFog
              : onOff(v.fog === 'on', p.cameraFog, {
                  name: 'on',
                  params: defaultsOf(Canvas3DParams.cameraFog.map('on')),
                }),
        };
        const said = (['occlusion', 'shadow', 'outline', 'fog'] as const)
          .filter((k) => v[k] !== 'keep')
          .map((k) => `${k} ${v[k]}`);
        if (!said.length)
          return ApiResult.fail('invalid', 'Change what?', 'occlusion, shadow, outline or fog: on or off.');
        const same =
          change.postprocessing.occlusion.name === p.postprocessing.occlusion.name &&
          change.postprocessing.shadow.name === p.postprocessing.shadow.name &&
          change.postprocessing.outline.name === p.postprocessing.outline.name &&
          change.cameraFog.name === p.cameraFog.name;
        if (same) return ApiResult.already('The lighting is already so.');
        return set(ctx, change as PartialCanvas3DProps, `Lighting: ${said.join(', ')}.`);
      },
    }),
  },
  {
    area: 'canvas',
    name: 'illumination',
    operation: typed({
      id: 'canvas.illumination',
      title: 'Turn illumination on or off',
      params: IlluminationParams,
      modules: canvasModule,
      async run({ on }, ctx) {
        const p = props(ctx);
        if (!p) return ApiResult.fail('unavailable', 'There is no canvas.');
        if (p.illumination.enabled === on) return ApiResult.already(`Illumination is already ${on ? 'on' : 'off'}.`);
        return set(ctx, { illumination: { enabled: on } }, `Illumination ${on ? 'on' : 'off'}.`);
      },
    }),
  },
  {
    area: 'canvas',
    name: 'camera',
    operation: typed({
      id: 'canvas.camera',
      title: 'Set the camera',
      params: CameraParams,
      modules: canvasModule,
      async run({ mode, fov }, ctx) {
        const p = props(ctx);
        if (!p) return ApiResult.fail('unavailable', 'There is no canvas.');
        const nextMode = mode === 'keep' ? p.camera.mode : mode;
        const nextFov = fov > 0 ? Math.max(10, fov) : p.camera.fov;
        if (nextMode === p.camera.mode && nextFov === p.camera.fov)
          return ApiResult.already('The camera is already so.');
        return set(ctx, { camera: { mode: nextMode, fov: nextFov } }, `Camera: ${nextMode}, ${nextFov}°.`);
      },
    }),
  },
  {
    area: 'canvas',
    name: 'clipping',
    operation: typed({
      id: 'canvas.clipping',
      title: 'Clip the scene',
      params: ClippingParams,
      modules: canvasModule,
      async run({ radius, far }, ctx) {
        const p = props(ctx);
        if (!p) return ApiResult.fail('unavailable', 'There is no canvas.');
        const r = Math.min(99, radius);
        if (p.cameraClipping.radius === r && p.cameraClipping.far === far) {
          return ApiResult.already('The clipping is already so.');
        }
        return set(ctx, { cameraClipping: { radius: r, far } }, `Clipping: ${radius}%${far ? ', far plane' : ''}.`);
      },
    }),
  },
];
