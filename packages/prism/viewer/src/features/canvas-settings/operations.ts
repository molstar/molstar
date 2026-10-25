/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Canvas3DParams, type Canvas3DProps, type PartialCanvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import { type ApiContext, ApiResult, type ApiFailure, type OperationDefinition } from '../../api/area.js';
import { OptionalParam, SwitchParam, type SwitchValue } from '../../api/params.js';
import { checkValues } from '../../api/run.js';
import { changeWords } from '../../api/change-words.js';

// The `canvas-settings` feature's operations (api.md §2, the `canvas` row; viewer-spec.md §3.1): added to the
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
  // WebXR's settings (Canvas3D's `xr`, the classic XRManager's params): used the next time the headset view is entered.
  xr: Canvas3DParams.xr,
  advanced: PD.Group({
    illumination: Canvas3DParams.illumination,
    multiSample: Canvas3DParams.multiSample,
    sharpening: Canvas3DParams.postprocessing.params.sharpening,
    bloom: Canvas3DParams.postprocessing.params.bloom,
    volumeResolutionScale: { ...Canvas3DParams.volumeResolutionScale, label: 'Volume resolution' },
  }),
};

export type CanvasSection = keyof typeof CanvasSections;

/** How a change of a field is said in the History ("Background white", "Shadows on"). */
const SectionNames: Readonly<Record<string, string>> = {
  color: 'Background',
  transparent: 'Transparent background',
  style: 'Background style',
  occlusion: 'Occlusion',
  shadow: 'Shadows',
  outline: 'Outline',
  fog: 'Fog',
  multiSample: 'Multi-sampling',
  volumeResolutionScale: 'Volume resolution',
};

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
    case 'xr':
      return p.xr as unknown as Record<string, unknown>;
    case 'advanced':
      return {
        illumination: p.illumination,
        multiSample: p.multiSample,
        sharpening: p.postprocessing.sharpening,
        bloom: p.postprocessing.bloom,
        volumeResolutionScale: p.volumeResolutionScale,
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
    case 'xr':
      return { xr: v as never };
    case 'advanced':
      return {
        illumination: v.illumination,
        multiSample: v.multiSample,
        postprocessing: { sharpening: v.sharpening, bloom: v.bloom },
        volumeResolutionScale: v.volumeResolutionScale,
      };
  }
}

const pp = Canvas3DParams.postprocessing.params;
/** A mapped param's 'on' case's params (Mol*'s own on/off props). */
const onParams = (m: PD.Any): PD.Params =>
  m.type === 'mapped' && m.map('on').type === 'group' ? (m.map('on') as PD.Group<any>).params : {};
const { enabled: _enabled, ...illuminationSettings } = Canvas3DParams.illumination.params;

export const SettingsParams = {
  section: PD.Select<CanvasSection>('background', PD.arrayToOptions(Object.keys(CanvasSections) as CanvasSection[])),
  values: PD.Value<Record<string, unknown>>({}),
};
/** The lighting's switches (api.md §1): each off, or on with its settings (given: set; absent: kept). */
export const LightingParams = {
  occlusion: SwitchParam(onParams(pp.occlusion), { description: 'Ambient occlusion: { name: "on", params } or off' }),
  shadow: SwitchParam(onParams(pp.shadow), { description: 'Shadows: { name: "on", params } or off' }),
  outline: SwitchParam(onParams(pp.outline), { description: 'Outline: { name: "on", params } or off' }),
  fog: SwitchParam(onParams(Canvas3DParams.cameraFog), { description: 'Fog: { name: "on", params } or off' }),
  illumination: SwitchParam(illuminationSettings, {
    description: 'Path-traced global illumination: { name: "on", params } or off',
  }),
};
export const CameraParams = {
  mode: OptionalParam(
    PD.Select<'perspective' | 'orthographic'>(
      'perspective',
      PD.arrayToOptions(['perspective', 'orthographic'] as const),
    ),
  ),
  fov: OptionalParam(
    PD.Numeric(45, { min: 10, max: 130, step: 1 }, { label: 'Field of view', description: 'Field of view in degrees' }),
  ),
};
export const XrParams = {
  on: PD.Boolean(true, { description: 'Enter (true) or leave (false) the headset view' }),
};
export const ClippingParams = {
  // As Canvas3D's `cameraClipping.radius`: how much is clipped away, around the camera's target.
  clip: OptionalParam(
    PD.Numeric(
      0,
      { min: 0, max: 99, step: 1 },
      { label: 'Clip (%)', description: 'Percent of the scene clipped away (0: none)' },
    ),
  ),
  far: OptionalParam(PD.Boolean(true, { label: 'Hide distant scene', description: 'Hide the scene in the distance' })),
};

const canvasModule = ['canvas'] as const;

/** A section's values given, checked against its params: the fields it has, a mapped one with the name of an option. */
function checkSection(section: CanvasSection, given: Record<string, unknown>): ApiFailure | undefined {
  const params = CanvasSections[section].params as PD.Params;
  const checked = checkValues(params, given);
  if (ApiResult.isFailure(checked)) return checked;
  for (const [key, value] of Object.entries(given)) {
    const p = params[key];
    if (!p || p.type !== 'mapped') continue;
    const names = p.select.options.map((o) => o[0]);
    const name = value && typeof value === 'object' ? (value as { name?: unknown }).name : undefined;
    if (typeof name !== 'string' || !names.includes(name)) {
      return ApiResult.fail('invalid', `'${key}' takes a name and params.`, `Names: ${names.join(', ')}.`);
    }
  }
  return undefined;
}

function props(ctx: ApiContext): Canvas3DProps | undefined {
  return ctx.plugin.canvas3d?.props;
}

function set(ctx: ApiContext, change: PartialCanvas3DProps, did: string) {
  ctx.plugin.setCanvas3DProps(change);
  return ApiResult.ok(did, { view: true });
}

/** A mapped param's 'on' case's defaults. */
const defaultsOf = (m: PD.Any) => (m.type === 'group' ? structuredClone(PD.getDefaultValues(m.params)) : {});

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
        const problem = checkSection(section, values);
        if (problem) return problem;
        const current = sectionValues(section, p);
        const next = { ...current, ...values };
        if (PD.areEqual(CanvasSections[section].params, current, next)) {
          return ApiResult.already(`The ${section} settings are already so.`);
        }
        const said = changeWords(CanvasSections[section].params, current, next, SectionNames);
        return set(
          ctx,
          sectionProps(section, next),
          `Changed the ${section} settings${said.length ? `: ${said.join(', ')}` : ''}.`,
        );
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
        const keys = (['occlusion', 'shadow', 'outline', 'fog', 'illumination'] as const).filter(
          (k) => v[k] !== undefined,
        );
        if (!keys.length) {
          return ApiResult.fail(
            'invalid',
            'Change what?',
            'occlusion, shadow, outline, fog or illumination: { name: "on" } (with params to set) or { name: "off" }.',
          );
        }
        // Absent: kept; off; on: its current settings (or the defaults), the given ones over them (api.md §1).
        const switched = (current: { name: string; params: any }, given: SwitchValue | undefined, m: PD.Any) =>
          given === undefined
            ? current
            : given.name === 'off'
              ? { name: 'off', params: {} }
              : {
                  name: 'on',
                  params: {
                    ...(current.name === 'on' ? current.params : defaultsOf(m.type === 'mapped' ? m.map('on') : m)),
                    ...(given.params ?? {}),
                  },
                };
        const illumination =
          v.illumination === undefined
            ? p.illumination
            : v.illumination.name === 'off'
              ? { ...p.illumination, enabled: false }
              : { ...p.illumination, ...(v.illumination.params ?? {}), enabled: true };
        const change = {
          postprocessing: {
            occlusion: switched(p.postprocessing.occlusion, v.occlusion, pp.occlusion),
            shadow: switched(p.postprocessing.shadow, v.shadow, pp.shadow),
            outline: switched(p.postprocessing.outline, v.outline, pp.outline),
          },
          cameraFog: switched(p.cameraFog, v.fog, Canvas3DParams.cameraFog),
          illumination,
        };
        const before = {
          postprocessing: {
            occlusion: p.postprocessing.occlusion,
            shadow: p.postprocessing.shadow,
            outline: p.postprocessing.outline,
          },
          cameraFog: p.cameraFog,
          illumination: p.illumination,
        };
        if (JSON.stringify(change) === JSON.stringify(before)) return ApiResult.already('The lighting is already so.');
        const said = keys.map((k) => `${k} ${v[k]!.name}`);
        return set(ctx, change as PartialCanvas3DProps, `Lighting: ${said.join(', ')}.`);
      },
    }),
  },
  // WebXR (Canvas3D's `xr`): the headset takes the view until it is left; the browser asks the person first.
  {
    area: 'canvas',
    name: 'xr',
    operation: typed({
      id: 'canvas.xr',
      title: 'Enter or leave XR',
      params: XrParams,
      modules: canvasModule,
      history: false,
      // Without a canvas the `canvas` module says so; with one, whether WebXR is there.
      when: (ctx) => {
        const xr = ctx.plugin.canvas3d?.xr;
        return !xr || xr.isSupported.value ? true : 'XR is not supported by this browser or device.';
      },
      async run({ on }, ctx) {
        const xr = ctx.plugin.canvas3d?.xr;
        if (!xr) return ApiResult.fail('unavailable', 'There is no canvas.');
        if (!xr.isSupported.value)
          return ApiResult.fail('unavailable', 'XR is not supported by this browser or device.');
        const presenting = xr.isPresenting.value;
        const next = on;
        if (next === presenting) return ApiResult.already(next ? 'Already in XR.' : 'Not in XR.');
        if (!next) {
          await xr.end();
          return ApiResult.ok('Left XR.', { view: true });
        }
        // `request` reports a failure on `requestFailed` rather than throwing.
        let failed: unknown;
        const sub = xr.requestFailed.subscribe((e) => {
          failed = e;
        });
        try {
          await xr.request();
        } finally {
          sub.unsubscribe();
        }
        if (failed !== undefined) {
          return ApiResult.fail(
            'failed',
            `Could not enter XR: ${failed instanceof Error ? failed.message : String(failed)}`,
          );
        }
        return ApiResult.ok('Entered XR.', { view: true });
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
      async run(v, ctx) {
        const { mode, fov } = v as { mode?: 'perspective' | 'orthographic'; fov?: number };
        const p = props(ctx);
        if (!p) return ApiResult.fail('unavailable', 'There is no canvas.');
        if (mode === undefined && fov === undefined) {
          return ApiResult.fail('invalid', 'Change what?', 'mode (perspective, orthographic) or fov (degrees).');
        }
        const nextMode = mode ?? p.camera.mode;
        const nextFov = fov ?? p.camera.fov;
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
      async run(v, ctx) {
        const { clip, far: farGiven } = v as { clip?: number; far?: boolean };
        const p = props(ctx);
        if (!p) return ApiResult.fail('unavailable', 'There is no canvas.');
        if (clip === undefined && farGiven === undefined) {
          return ApiResult.fail('invalid', 'Change what?', 'clip (percent clipped away) or far (true or false).');
        }
        // Absent: kept (api.md §1).
        const radius = clip === undefined ? p.cameraClipping.radius : Math.max(0, Math.min(99, clip));
        const far = farGiven ?? p.cameraClipping.far;
        if (p.cameraClipping.radius === radius && p.cameraClipping.far === far) {
          return ApiResult.already('The clipping is already so.');
        }
        const said = radius ? `${radius}% clipped away` : 'nothing clipped';
        return set(ctx, { cameraClipping: { radius, far } }, `Clipping: ${said}${far ? ', far plane on' : ''}.`);
      },
    }),
  },
];
