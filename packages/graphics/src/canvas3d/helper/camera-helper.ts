import * as ShapeUtils from '@molstar/graphics/geo/shape/shape';
/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { produce } from '@molstar/core/util/produce';
import { ValueCell } from '@molstar/core/util/value-cell';
import { Interval } from '@molstar/core/data/int/interval';
import { addSimpleCylinder } from '@molstar/graphics/geo/geometry/mesh/builder/cylinder';
import { addSphere } from '@molstar/graphics/geo/geometry/mesh/builder/sphere';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { Text } from '@molstar/graphics/geo/geometry/text/text';
import { TextBuilder } from '@molstar/graphics/geo/geometry/text/text-builder';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { DataLoci, EmptyLoci, isEveryLoci, type Loci } from '@molstar/model/model/loci';
import type { Shape } from '@molstar/model/model/shape';
import { Visual } from '@molstar/graphics/repr/visual';
import { ColorNames } from '@molstar/core/util/color/names';
import { type MarkerAction, MarkerActions } from '@molstar/core/util/marker-action';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Representation } from '@molstar/graphics/repr/representation';
import { ScaleLineParams } from './scale-line.js';
import type { ScreenHelper } from './screen-helper.js';

// TODO add a grid

const AxesParams = {
  alpha: PD.Numeric(0.51, { min: 0, max: 1, step: 0.01 }, { isEssential: true, label: 'Opacity' }),
  colorX: PD.Color(ColorNames.red, { isEssential: true }),
  colorY: PD.Color(ColorNames.green, { isEssential: true }),
  colorZ: PD.Color(ColorNames.blue, { isEssential: true }),
  scale: PD.Numeric(0.33, { min: 0.1, max: 2, step: 0.1 }, { isEssential: true }),
  location: PD.Select('bottom-left', PD.arrayToOptions(Representation.ScreenEdgeAnchors)),
  locationOffsetX: PD.Numeric(0),
  locationOffsetY: PD.Numeric(0),
  originColor: PD.Color(ColorNames.grey),
  radiusScale: PD.Numeric(0.075, { min: 0.01, max: 0.3, step: 0.001 }),
  showPlanes: PD.Boolean(true),
  planeColorXY: PD.Color(ColorNames.grey, { label: 'Plane Color XY' }),
  planeColorXZ: PD.Color(ColorNames.grey, { label: 'Plane Color XZ' }),
  planeColorYZ: PD.Color(ColorNames.grey, { label: 'Plane Color YZ' }),
  showLabels: PD.Boolean(false),
  labelX: PD.Text('X'),
  labelY: PD.Text('Y'),
  labelZ: PD.Text('Z'),
  labelColorX: PD.Color(ColorNames.grey),
  labelColorY: PD.Color(ColorNames.grey),
  labelColorZ: PD.Color(ColorNames.grey),
  labelOpacity: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
  labelScale: PD.Numeric(0.25, { min: 0.1, max: 1.0, step: 0.01 }),
};
type AxesParams = typeof AxesParams;
type AxesProps = PD.Values<AxesParams>;

export const CameraHelperParams = {
  axes: PD.MappedStatic(
    'on',
    {
      on: PD.Group(AxesParams),
      off: PD.Group({}),
    },
    { cycle: true, description: 'Show camera orientation axes' },
  ),
  scale: PD.MappedStatic(
    'on',
    {
      on: PD.Group(ScaleLineParams),
      off: PD.Group({}),
    },
    { cycle: true, description: 'Show a scale line: a length in Å at the camera target (`ScaleLine`)' },
  ),
};
export type CameraHelperParams = typeof CameraHelperParams;
export type CameraHelperProps = PD.Values<CameraHelperParams>;
type CameraAxesProps = Pick<CameraHelperProps, 'axes'>;

/**
 * The corner axes (`camera.helper.axes` of the canvas): the camera's orientation in a corner of the screen layer
 * (`ScreenHelper`), turning with the camera; a click on an axis or a plane is a `CameraAxesLoci` (the plugin aligns
 * the camera to it). In CSS pixels; three objects, the axes, the planes and the labels (one mesh of overlapping parts
 * would show its insides through its translucency: the transparent pass writes no depth).
 */
export class CameraAxes {
  props: CameraAxesProps = {
    axes: { name: 'off', params: {} },
  };

  private meshRenderObject: GraphicsRenderObject | undefined;
  private planesRenderObject: GraphicsRenderObject | undefined;
  private textRenderObject: GraphicsRenderObject | undefined;

  constructor(
    private readonly screen: ScreenHelper,
    props: Partial<CameraHelperProps> = {},
  ) {
    this.setProps(props);
  }

  setProps(props: Partial<CameraAxesProps>) {
    this.props = produce(this.props, (p) => {
      if (props.axes === undefined) return;
      p.axes.name = props.axes.name;
      if (props.axes.name === 'off') {
        this.screen.release(this);
        this.meshRenderObject = this.planesRenderObject = this.textRenderObject = undefined;
        return;
      }
      const params = props.axes.params;
      this.meshRenderObject = createMeshRenderObject(params);
      this.planesRenderObject = params.showPlanes ? createPlanesRenderObject(params) : undefined;
      this.textRenderObject = params.showLabels ? createTextRenderObject(params) : undefined;

      // inward by the labels' radius, or the axes' (as their place has always been)
      const inset = (this.textRenderObject ?? this.meshRenderObject).values.boundingSphere.ref.value.radius;
      const scene = this.screen.getScene(this, {
        anchor: params.location,
        offsetX: params.locationOffsetX,
        offsetY: params.locationOffsetY,
        followRotation: true,
        inset,
      });
      scene.clear();
      for (const o of this.renderObjects) scene.add(o);
      scene.commit();

      p.axes.params = { ...params };
    });
    // made again at their `alpha`: still hovered, at full
    if (this.hovered) this.applyAlpha();
  }

  get renderObjects(): GraphicsRenderObject[] {
    const objects: GraphicsRenderObject[] = [];
    for (const o of [this.meshRenderObject, this.planesRenderObject, this.textRenderObject]) if (o) objects.push(o);
    return objects;
  }

  get isEnabled() {
    return this.props.axes.name !== 'off';
  }

  private hovered = false;

  /**
   * At full alpha while the pointer is over them, else at their `alpha` (dimmed, as navigation widgets are): the axes
   * and planes' `uAlpha`, no rebuild; true when it changed (a redraw)
   */
  setHovered(hovered: boolean) {
    if (this.hovered === hovered) return false;
    this.hovered = hovered;
    this.applyAlpha();
    return this.isEnabled;
  }

  private applyAlpha() {
    if (this.props.axes.name === 'off') return;
    const alpha = this.hovered ? 1 : this.props.axes.params.alpha;
    for (const o of [this.meshRenderObject, this.planesRenderObject]) {
      if (o) ValueCell.updateIfChanged(o.values.uAlpha, alpha);
    }
  }

  getLoci(pickingId: PickingId) {
    const { objectId, groupId, instanceId } = pickingId;
    if (
      (!this.meshRenderObject || objectId !== this.meshRenderObject.id) &&
      (!this.planesRenderObject || objectId !== this.planesRenderObject.id) &&
      (!this.textRenderObject || objectId !== this.textRenderObject.id)
    )
      return EmptyLoci;
    if (groupId === CameraHelperAxis.None) return EmptyLoci;
    return CameraAxesLoci(this, groupId, instanceId);
  }

  private eachGroup = (loci: Loci, apply: (interval: Interval) => boolean): boolean => {
    if (!isCameraAxesLoci(loci)) return false;
    let changed = false;
    for (const ro of [this.meshRenderObject, this.planesRenderObject, this.textRenderObject]) {
      if (!ro) continue;
      const groupCount = ro.values.uGroupCount.ref.value;
      for (const { groupId, instanceId } of loci.elements) {
        const idx = instanceId * groupCount + groupId;
        if (apply(Interval.ofSingleton(idx))) changed = true;
      }
    }
    return changed;
  };

  mark(loci: Loci, action: MarkerAction) {
    if (!MarkerActions.is(MarkerActions.Highlighting, action)) return false;
    if (!isEveryLoci(loci)) {
      if (!isCameraAxesLoci(loci)) return false;
      if (loci.data !== this) return false;
    }
    // each one marked (not `||`: the first marked would skip the others)
    const axes = Visual.mark(this.meshRenderObject, loci, action, this.eachGroup);
    const planes = Visual.mark(this.planesRenderObject, loci, action, this.eachGroup);
    const text = Visual.mark(this.textRenderObject, loci, action, this.eachGroup);
    return axes || planes || text;
  }
}

export enum CameraHelperAxis {
  None = 0,
  X,
  Y,
  Z,
  XY,
  XZ,
  YZ,
  Origin,
}

function getAxisLabel(axis: number, cameraHelper: CameraAxes) {
  const a = cameraHelper.props.axes;
  const x = a.name === 'on' ? a.params.labelX : 'X';
  const y = a.name === 'on' ? a.params.labelY : 'Y';
  const z = a.name === 'on' ? a.params.labelZ : 'Z';
  switch (axis) {
    case CameraHelperAxis.X:
      return `${x} Axis`;
    case CameraHelperAxis.Y:
      return `${y} Axis`;
    case CameraHelperAxis.Z:
      return `${z} Axis`;
    case CameraHelperAxis.XY:
      return `${x}${y} Plane`;
    case CameraHelperAxis.XZ:
      return `${x}${z} Plane`;
    case CameraHelperAxis.YZ:
      return `${y}${z} Plane`;
    case CameraHelperAxis.Origin:
      return 'Origin';
    default:
      return 'Axes';
  }
}

function CameraAxesLoci(cameraHelper: CameraAxes, groupId: number, instanceId: number) {
  return DataLoci('camera-axes', cameraHelper, [{ groupId, instanceId }], void 0 /** bounding sphere */, () =>
    getAxisLabel(groupId, cameraHelper),
  );
}
export type CameraAxesLoci = ReturnType<typeof CameraAxesLoci>;
export function isCameraAxesLoci(x: Loci): x is CameraAxesLoci {
  return x.kind === 'data-loci' && x.tag === 'camera-axes';
}

/** The axes' length, short of the labels where they are shown */
function axisLength(props: AxesProps) {
  const textScale = props.showLabels ? (100 * props.labelScale) / 3 : 0;
  return 100 * props.scale - textScale;
}

/**
 * The axes: a slim arrow each (a shaft and a small cone head, ending where the old rods' end caps did) from a small
 * origin sphere; no part inside another, so a translucent arrow shows only its front (`transparentBackfaces` off)
 */
function createAxesMesh(props: AxesProps, mesh?: Mesh) {
  const state = MeshBuilder.createState(512, 256, mesh);
  const length = axisLength(props);
  const radius = props.radiusScale * 100 * props.scale * 0.6;
  const head = Math.min(radius * 4.5, length * 0.25);
  // closed at the origin's end: a translucent open tube would show its inside
  const shaft = { radiusTop: radius, radiusBottom: radius, radialSegments: 24, bottomCap: true };
  const cone = { radiusTop: 0, radiusBottom: radius * 2.2, radialSegments: 24, bottomCap: true };

  // the origin: a small sphere the shafts start at, not inside it (an overlap shows through the translucency)
  const origin = radius * 1.5;
  state.currentGroup = CameraHelperAxis.Origin;
  addSphere(state, Vec3.origin, origin, 3);

  const axes = [
    [CameraHelperAxis.X, Vec3.unitX],
    [CameraHelperAxis.Y, Vec3.unitY],
    [CameraHelperAxis.Z, Vec3.unitZ],
  ] as const;
  for (const [group, unit] of axes) {
    state.currentGroup = group;
    const tip = Vec3.scale(Vec3(), unit, length + radius);
    const base = Vec3.scale(Vec3(), unit, length + radius - head);
    // `addSimpleCylinder` takes the start's radius as `radiusBottom` (and caps it): from the origin sphere's surface
    addSimpleCylinder(state, Vec3.scale(Vec3(), unit, origin), base, shaft);
    // `addSimpleCylinder` puts `radiusBottom` at its start: the head's base at the shaft's end
    addSimpleCylinder(state, base, tip, cone);
  }

  return MeshBuilder.getMesh(state);
}

/**
 * The planes, translucent: a square each between its two axes, clear of them; both windings (only the face toward the
 * viewer is drawn, so a plane does not darken itself)
 */
function createPlanesMesh(props: AxesProps, mesh?: Mesh) {
  const state = MeshBuilder.createState(64, 32, mesh);
  const length = axisLength(props);
  const from = length * 0.12;
  const to = length * 0.5;
  const square = (group: number, u: Vec3, v: Vec3) => {
    state.currentGroup = group;
    const p = (a: number, b: number) => Vec3.add(Vec3(), Vec3.scale(Vec3(), u, a), Vec3.scale(Vec3(), v, b));
    MeshBuilder.addTriangle(state, p(from, from), p(to, from), p(to, to));
    MeshBuilder.addTriangle(state, p(from, from), p(to, to), p(from, to));
    MeshBuilder.addTriangle(state, p(from, from), p(to, to), p(to, from));
    MeshBuilder.addTriangle(state, p(from, from), p(from, to), p(to, to));
  };
  square(CameraHelperAxis.XY, Vec3.unitX, Vec3.unitY);
  square(CameraHelperAxis.XZ, Vec3.unitZ, Vec3.unitX);
  square(CameraHelperAxis.YZ, Vec3.unitY, Vec3.unitZ);
  return MeshBuilder.getMesh(state);
}

function getAxesMeshShape(props: AxesProps, shape?: Shape<Mesh>, planes = false) {
  const scale = 100 * props.scale;
  const mesh = planes ? createPlanesMesh(props, shape?.geometry) : createAxesMesh(props, shape?.geometry);
  mesh.setBoundingSphere(Sphere3D.create(Vec3.create(scale / 2, scale / 2, scale / 2), scale + scale / 4));
  const getColor = (groupId: number) => {
    switch (groupId) {
      case CameraHelperAxis.X:
        return props.colorX;
      case CameraHelperAxis.Y:
        return props.colorY;
      case CameraHelperAxis.Z:
        return props.colorZ;
      case CameraHelperAxis.XY:
        return props.planeColorXY;
      case CameraHelperAxis.XZ:
        return props.planeColorXZ;
      case CameraHelperAxis.YZ:
        return props.planeColorYZ;
      case CameraHelperAxis.Origin:
        return props.originColor;
      default:
        return ColorNames.grey;
    }
  };
  return ShapeUtils.create(
    'axes-mesh',
    {},
    mesh,
    getColor,
    () => 1,
    () => '',
  );
}

function createMeshRenderObject(props: AxesProps) {
  const shape = getAxesMeshShape(props);
  return ShapeUtils.createRenderObject(shape, {
    ...PD.getDefaultValues(Mesh.Params),
    ...props,
    // shaded: an arrow toward the viewer reads as one, not as a flat disk
    ignoreLight: false,
    cellSize: 0,
  });
}

function createPlanesRenderObject(props: AxesProps) {
  const shape = getAxesMeshShape(props, undefined, true);
  return ShapeUtils.createRenderObject(shape, {
    ...PD.getDefaultValues(Mesh.Params),
    ...props,
    ignoreLight: true,
    cellSize: 0,
  });
}

//

function createAxesText(props: AxesProps, text?: Text) {
  const builder = TextBuilder.create(props, 8, 8, text);
  const scale = 100 * props.scale;

  const x = Vec3.scale(Vec3(), Vec3.unitX, scale);
  const y = Vec3.scale(Vec3(), Vec3.unitY, scale);
  const z = Vec3.scale(Vec3(), Vec3.unitZ, scale);

  const textScale = 100 * props.labelScale;
  builder.add(props.labelX, x[0], x[1], x[2], 0.0, textScale, CameraHelperAxis.X);
  builder.add(props.labelY, y[0], y[1], y[2], 0.0, textScale, CameraHelperAxis.Y);
  builder.add(props.labelZ, z[0], z[1], z[2], 0.0, textScale, CameraHelperAxis.Z);

  return builder.getText();
}

function getAxesTextShape(props: AxesProps, shape?: Shape<Text>) {
  const scale = 100 * props.scale;
  const text = createAxesText(props, shape?.geometry);
  text.setBoundingSphere(Sphere3D.create(Vec3.create(scale / 2, scale / 2, scale / 2), scale));
  const getColor = (groupId: number) => {
    switch (groupId) {
      case CameraHelperAxis.X:
        return props.labelColorX;
      case CameraHelperAxis.Y:
        return props.labelColorY;
      case CameraHelperAxis.Z:
        return props.labelColorZ;
      default:
        return ColorNames.grey;
    }
  };
  return ShapeUtils.create(
    'axes-text',
    {},
    text,
    getColor,
    () => 1,
    () => '',
  );
}

function createTextRenderObject(props: AxesProps) {
  const shape = getAxesTextShape(props);
  return ShapeUtils.createRenderObject(shape, {
    ...PD.getDefaultValues(Text.Params),
    ...props,
    alpha: props.labelOpacity,
    cellSize: 0,
  });
}
