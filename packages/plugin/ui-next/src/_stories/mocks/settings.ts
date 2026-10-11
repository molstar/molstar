/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ParamsFormModel } from '../../forms/params-form.js';
import type { SettingsGroup, SettingsModel } from '../../model/settings.js';

const onOff = (on: boolean, params: PD.Params = {}, info: PD.Info = {}) =>
  PD.MappedStatic(on ? 'on' : 'off', { on: PD.Group(params), off: PD.Group({}) }, info);

// Canvas3D / Canvas3DContext params (a representative subset), mapped into groups by intent.
const Appearance = {
  backgroundColor: PD.Color(0xffffff as any, { label: 'Background' }),
  transparentBackground: PD.Boolean(false, { label: 'Transparent background' }),
  checkeredTransparentBackground: PD.Boolean(true, { label: 'Checkered when transparent' }),
  occlusion: onOff(
    true,
    {
      samples: PD.Numeric(32, { min: 1, max: 256, step: 1 }),
      radius: PD.Numeric(5, { min: 0, max: 20, step: 0.1 }),
      bias: PD.Numeric(0.8, { min: 0, max: 3, step: 0.1 }),
    },
    { label: 'Ambient occlusion' },
  ),
  shadow: onOff(false, { steps: PD.Numeric(1, { min: 1, max: 64, step: 1 }) }, { label: 'Shadows' }),
  outline: onOff(false, { scale: PD.Numeric(1, { min: 1, max: 5, step: 1 }) }, { label: 'Outline' }),
  cameraFog: onOff(true, { intensity: PD.Numeric(15, { min: 1, max: 100, step: 1 }) }, { label: 'Fog' }),
  marking: PD.Group(
    {
      enabled: PD.Boolean(true),
      highlightEdgeColor: PD.Color(0x999999 as any, { label: 'Highlight edge' }),
      selectEdgeColor: PD.Color(0xffff00 as any, { label: 'Select edge' }),
      ghostEdgeStrength: PD.Numeric(0.3, { min: 0, max: 1, step: 0.01 }),
    },
    { label: 'Marking' },
  ),
};

const Camera = {
  mode: PD.Select('perspective', [
    ['perspective', 'Perspective'],
    ['orthographic', 'Orthographic'],
  ]),
  fov: PD.Numeric(45, { min: 10, max: 130, step: 1 }, { label: 'Field of view' }),
  clipping: PD.Group(
    {
      radius: PD.Numeric(100, { min: 0, max: 100, step: 1 }),
      far: PD.Boolean(true),
      minNear: PD.Numeric(5, { min: 0.1, max: 100, step: 0.1 }, { label: 'Min near' }),
    },
    { label: 'Clipping' },
  ),
  cameraResetDurationMs: PD.Numeric(250, { min: 0, max: 1000, step: 1 }, { label: 'Reset duration (ms)' }),
  sceneRadiusFactor: PD.Numeric(1, { min: 1, max: 10, step: 0.1 }, { label: 'Scene radius factor' }),
  manualReset: PD.Boolean(false, {
    label: 'Manual reset',
    description: 'Do not reset the camera when the scene changes.',
  }),
};

const Interaction = {
  pickPadding: PD.Numeric(3, { min: 0, max: 10, step: 1 }, { label: 'Pick padding (px)' }),
  userInteractionReleaseMs: PD.Numeric(250, { min: 0, max: 1000, step: 1 }, { label: 'Interaction release (ms)' }),
  trackball: PD.Group(
    {
      rotateSpeed: PD.Numeric(5, { min: 1, max: 10, step: 1 }),
      zoomSpeed: PD.Numeric(7, { min: 1, max: 15, step: 1 }),
      panSpeed: PD.Numeric(1, { min: 0.1, max: 5, step: 0.1 }),
      staticMoving: PD.Boolean(true, { label: 'No inertia' }),
    },
    { label: 'Trackball' },
  ),
  pointer: PD.Group({ enabled: PD.Boolean(true) }, { label: 'Pointer' }),
};

const Performance = {
  resolutionMode: PD.Select('auto', [
    ['auto', 'Auto'],
    ['scaled', 'Scaled'],
    ['native', 'Native'],
  ]),
  pixelScale: PD.Numeric(1, { min: 0.1, max: 2, step: 0.05 }),
  pickScale: PD.Numeric(0.25, { min: 0.1, max: 1, step: 0.05 }),
  transparency: PD.Select('blended', [
    ['blended', 'Weighted, blended'],
    ['wboit', 'Weighted, order-independent'],
    ['dpoit', 'Depth peeling'],
  ]),
  dpoitIterations: PD.Numeric(2, { min: 1, max: 10, step: 1 }, { label: 'Depth peeling iterations' }),
  multiSample: onOff(true, { sampleLevel: PD.Numeric(2, { min: 0, max: 5, step: 1 }) }, { label: 'Multi-sampling' }),
  hiZ: PD.Boolean(false, { label: 'Hi-Z culling' }),
};

const Interface = {
  theme: PD.Select('auto', [
    ['auto', 'Follow the system'],
    ['light', 'Light'],
    ['dark', 'Dark'],
    ['blue', 'Blue'],
  ]),
  layout: PD.Select('auto', [
    ['auto', 'Auto (compact when narrow)'],
    ['desktop', 'Desktop'],
    ['compact', 'Compact'],
  ]),
  pinPanel: PD.Boolean(false, { label: 'Pin the panel', description: 'Keep the panel beside the canvas.' }),
  sequence: PD.Boolean(true, { label: 'Sequence', description: 'Show the sequence above the canvas.' }),
  tooltips: PD.Boolean(true, { label: 'Tooltips' }),
};

// plugin.state.behaviors, by category: each a switch, and its params behind the sliders.
const Behaviors = {
  representation: PD.Group(
    {
      chemicalComponentDictionary: onOff(true, {}, { label: 'wwPDB Chemical Component Dictionary' }),
    },
    { label: 'Representation', isExpanded: true },
  ),
  interaction: PD.Group(
    {
      highlightLoci: onOff(true, { mark: PD.Boolean(true) }, { label: 'Highlight Loci on Canvas' }),
      selectLoci: onOff(true, { mark: PD.Boolean(true) }, { label: 'Select Loci on Canvas' }),
      defaultLociLabel: onOff(true, {}, { label: 'Provide Default Loci Label' }),
      focusLoci: onOff(true, {}, { label: 'Representation Focus Loci' }),
      cameraFocusLoci: onOff(
        true,
        {
          minRadius: PD.Numeric(8, { min: 1, max: 50, step: 1 }),
          extraRadius: PD.Numeric(4, { min: 0, max: 20, step: 0.5 }),
          durationMs: PD.Numeric(250, { min: 0, max: 1000, step: 1 }),
        },
        { label: 'Camera Focus Loci on Canvas' },
      ),
      axisHelper: onOff(true, {}, { label: 'Camera Axis Helper' }),
      cameraControls: onOff(true, {}, { label: 'Camera Controls on Canvas' }),
      snapshotControls: onOff(true, {}, { label: 'Snapshot Controls' }),
      structureFocus: onOff(
        true,
        { expandRadius: PD.Numeric(5, { min: 1, max: 10, step: 1 }) },
        { label: 'Structure Focus Representation' },
      ),
    },
    { label: 'Interaction', isExpanded: true },
  ),
  customProperties: PD.Group(
    {
      structureInfo: onOff(true, {}, { label: 'Structure Info' }),
      accessibleSurfaceArea: onOff(true, {}, { label: 'Accessible Surface Area' }),
      sifts: onOff(true, {}, { label: 'SIFTS Mapping' }),
      interactions: onOff(true, {}, { label: 'Interactions' }),
      validation: onOff(false, {}, { label: 'wwPDB Validation Report' }),
    },
    { label: 'Custom Properties', isExpanded: true },
  ),
};

const Advanced = {
  renderer: PD.Group(
    {
      exposure: PD.Numeric(1, { min: 0, max: 3, step: 0.01 }),
      ambientIntensity: PD.Numeric(0.4, { min: 0, max: 2, step: 0.01 }),
      lightIntensity: PD.Numeric(0.6, { min: 0, max: 2, step: 0.01 }),
      backgroundOpacity: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
    },
    { label: 'Renderer' },
  ),
  handle: onOff(false, {}, { label: 'Handle' }),
  xr: PD.Group(
    {
      minTargetDistance: PD.Numeric(0.4, { min: 0, max: 1, step: 0.01 }),
      disablePostprocessing: PD.Boolean(true),
    },
    { label: 'XR' },
  ),
};

/** The Settings tab's groups for the stories. */
export function mockSettings(): SettingsModel {
  const group = (
    id: string,
    title: string,
    icon: SettingsGroup['icon'],
    params: PD.Params,
    wideLabels?: boolean,
  ): SettingsGroup => ({ id, title, icon, form: new ParamsFormModel(params), wideLabels });
  const groups = [
    group('appearance', 'Appearance', 'palette', Appearance),
    group('camera', 'Camera', 'camera', Camera),
    group('interaction', 'Interaction', 'mouse-pointer-2', Interaction),
    group('performance', 'Performance', 'gauge', Performance),
    group('interface', 'Interface', 'app-window', Interface),
    group('behaviors', 'Behaviors', 'toggle-right', Behaviors, true),
    group('advanced', 'Advanced', 'wrench', Advanced),
  ];
  return { groups };
}
