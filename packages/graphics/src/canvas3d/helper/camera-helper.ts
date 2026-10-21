import * as ShapeUtils from '@molstar/graphics/geo/shape/shape';
/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { produce } from '@molstar/core/util/produce';
import { Interval } from '@molstar/core/data/int/interval';
import { addSimpleCylinder } from '@molstar/graphics/geo/geometry/mesh/builder/cylinder';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import type { PickingId } from '@molstar/graphics/geo/geometry/picking';
import { Text } from '@molstar/graphics/geo/geometry/text/text';
import { TextBuilder } from '@molstar/graphics/geo/geometry/text/text-builder';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { Scene } from '@molstar/graphics/gl/scene';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { DataLoci, EmptyLoci, isEveryLoci, type Loci } from '@molstar/model/model/loci';
import type { Shape } from '@molstar/model/model/shape';
import { Visual } from '@molstar/graphics/repr/visual';
import { ColorNames } from '@molstar/core/util/color/names';
import { type MarkerAction, MarkerActions } from '@molstar/core/util/marker-action';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { assertUnreachable } from '@molstar/core/util/type-helpers';
import { Camera, type ICamera } from '../camera.js';
import type { Viewport } from '../camera/util.js';

// TODO add scale line/grid

const AxesParams = {
  alpha: PD.Numeric(0.51, { min: 0, max: 1, step: 0.01 }, { isEssential: true, label: 'Opacity' }),
  colorX: PD.Color(ColorNames.red, { isEssential: true }),
  colorY: PD.Color(ColorNames.green, { isEssential: true }),
  colorZ: PD.Color(ColorNames.blue, { isEssential: true }),
  scale: PD.Numeric(0.33, { min: 0.1, max: 2, step: 0.1 }, { isEssential: true }),
  location: PD.Select(
    'bottom-left',
    PD.arrayToOptions(['bottom-left', 'bottom-right', 'top-left', 'top-right'] as const),
  ),
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
};
export type CameraHelperParams = typeof CameraHelperParams;
export type CameraHelperProps = PD.Values<CameraHelperParams>;

export class CameraHelper {
  scene: Scene;
  camera: Camera;
  props: CameraHelperProps = {
    axes: { name: 'off', params: {} },
  };

  /** The axes and the planes, each its own (one mesh of overlapping parts showed its insides through its translucency) */
  private meshRenderObject: GraphicsRenderObject | undefined;
  private planesRenderObject: GraphicsRenderObject | undefined;
  private textRenderObject: GraphicsRenderObject | undefined;
  private pixelRatio = 1;

  constructor(
    private graphics: GraphicsContext,
    props: Partial<CameraHelperProps> = {},
  ) {
    this.scene = Scene.create(graphics, 'blended');

    this.camera = new Camera();
    Vec3.set(this.camera.up, 0, 1, 0);
    Vec3.set(this.camera.target, 0, 0, 0);

    this.setProps(props);
  }

  setContext(graphics: GraphicsContext) {
    this.graphics = graphics;
    this.scene.setContext(graphics);
  }

  setProps(props: Partial<CameraHelperProps>) {
    this.props = produce(this.props, (p) => {
      if (props.axes !== undefined) {
        p.axes.name = props.axes.name;
        if (props.axes.name === 'on') {
          this.scene.clear();
          this.pixelRatio = this.graphics.pixelRatio;
          const params = {
            ...props.axes.params,
            scale: props.axes.params.scale * this.pixelRatio,
            labelScale: props.axes.params.labelScale * this.pixelRatio,
          };
          this.meshRenderObject = createMeshRenderObject(params);
          this.scene.add(this.meshRenderObject);
          if (props.axes.params.showPlanes) {
            this.planesRenderObject = createPlanesRenderObject(params);
            this.scene.add(this.planesRenderObject);
          } else {
            this.planesRenderObject = undefined;
          }
          if (props.axes.params.showLabels) {
            this.textRenderObject = createTextRenderObject(params);
            this.scene.add(this.textRenderObject);
          } else {
            this.textRenderObject = undefined;
          }
          this.scene.commit();

          Vec3.set(this.camera.position, 0, 0, params.scale * 200);
          Mat4.lookAt(this.camera.view, this.camera.position, this.camera.target, this.camera.up);

          p.axes.params = { ...props.axes.params };
        }
      }
    });
  }

  get isEnabled() {
    return this.props.axes.name === 'on';
  }

  getLoci(pickingId: PickingId) {
    const { objectId, groupId, instanceId } = pickingId;
    if (
      ((!this.meshRenderObject || objectId !== this.meshRenderObject.id) &&
        (!this.planesRenderObject || objectId !== this.planesRenderObject.id) &&
        (!this.textRenderObject || objectId !== this.textRenderObject.id)) ||
      groupId === CameraHelperAxis.None
    )
      return EmptyLoci;
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

  update(camera: ICamera) {
    if (!this.meshRenderObject || this.props.axes.name === 'off') return;

    if (this.pixelRatio !== this.graphics.pixelRatio) {
      this.setProps(this.props);
    }

    updateCamera(this.camera, camera.viewport, camera.viewOffset);
    Mat4.extractRotation(this.scene.view, camera.view);

    const r = this.textRenderObject
      ? this.textRenderObject.values.boundingSphere.ref.value.radius
      : this.meshRenderObject.values.boundingSphere.ref.value.radius;
    const l = this.props.axes.params.location;
    const ox = this.props.axes.params.locationOffsetX * this.pixelRatio;
    const oy = this.props.axes.params.locationOffsetY * this.pixelRatio;
    if (l === 'bottom-left') {
      Mat4.setTranslation(
        this.scene.view,
        Vec3.create(-camera.viewport.width / 2 + r + ox, -camera.viewport.height / 2 + r + oy, 0),
      );
    } else if (l === 'bottom-right') {
      Mat4.setTranslation(
        this.scene.view,
        Vec3.create(camera.viewport.width / 2 - r - ox, -camera.viewport.height / 2 + r + oy, 0),
      );
    } else if (l === 'top-left') {
      Mat4.setTranslation(
        this.scene.view,
        Vec3.create(-camera.viewport.width / 2 + r + ox, camera.viewport.height / 2 - r - oy, 0),
      );
    } else if (l === 'top-right') {
      Mat4.setTranslation(
        this.scene.view,
        Vec3.create(camera.viewport.width / 2 - r - ox, camera.viewport.height / 2 - r - oy, 0),
      );
    } else {
      assertUnreachable(l);
    }
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

function getAxisLabel(axis: number, cameraHelper: CameraHelper) {
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

function CameraAxesLoci(cameraHelper: CameraHelper, groupId: number, instanceId: number) {
  return DataLoci('camera-axes', cameraHelper, [{ groupId, instanceId }], void 0 /** bounding sphere */, () =>
    getAxisLabel(groupId, cameraHelper),
  );
}
export type CameraAxesLoci = ReturnType<typeof CameraAxesLoci>;
export function isCameraAxesLoci(x: Loci): x is CameraAxesLoci {
  return x.kind === 'data-loci' && x.tag === 'camera-axes';
}

function updateCamera(camera: Camera, viewport: Viewport, viewOffset: Camera.ViewOffset) {
  const { near, far } = camera;

  const fullLeft = -viewport.width / 2;
  const fullRight = viewport.width / 2;
  const fullTop = viewport.height / 2;
  const fullBottom = -viewport.height / 2;

  const dx = (fullRight - fullLeft) / 2;
  const dy = (fullTop - fullBottom) / 2;
  const cx = (fullRight + fullLeft) / 2;
  const cy = (fullTop + fullBottom) / 2;

  let left = cx - dx;
  let right = cx + dx;
  let top = cy + dy;
  let bottom = cy - dy;

  if (viewOffset.enabled) {
    const scaleW = (fullRight - fullLeft) / viewOffset.width;
    const scaleH = (fullTop - fullBottom) / viewOffset.height;
    left += scaleW * viewOffset.offsetX;
    right = left + scaleW * viewOffset.width;
    top -= scaleH * viewOffset.offsetY;
    bottom = top - scaleH * viewOffset.height;
  }

  Mat4.ortho(camera.projection, left, right, top, bottom, near, far);
}

/** The axes' length, short of the labels where they are shown */
function axisLength(props: AxesProps) {
  const textScale = props.showLabels ? (100 * props.labelScale) / 3 : 0;
  return 100 * props.scale - textScale;
}

/** Where the axes and the planes start from the origin: the corner between them left open */
const Inset = 0.12;

/**
 * The axes: a slim arrow each (a shaft and a small cone head, ending where the old rods' end caps did), starting as far
 * from the origin as the planes do (`Inset`); no part inside another, so a translucent arrow shows only its front
 * (`transparentBackfaces` off)
 */
function createAxesMesh(props: AxesProps, mesh?: Mesh) {
  const state = MeshBuilder.createState(512, 256, mesh);
  const length = axisLength(props);
  const radius = props.radiusScale * 100 * props.scale * 0.6;
  const head = Math.min(radius * 4.5, length * 0.25);
  // closed at the origin's end: a translucent open tube would show its inside
  const shaft = { radiusTop: radius, radiusBottom: radius, radialSegments: 24, bottomCap: true };
  const cone = { radiusTop: 0, radiusBottom: radius * 2.2, radialSegments: 24, bottomCap: true };

  const start = length * Inset;

  const axes = [
    [CameraHelperAxis.X, Vec3.unitX],
    [CameraHelperAxis.Y, Vec3.unitY],
    [CameraHelperAxis.Z, Vec3.unitZ],
  ] as const;
  for (const [group, unit] of axes) {
    state.currentGroup = group;
    const tip = Vec3.scale(Vec3(), unit, length + radius);
    const base = Vec3.scale(Vec3(), unit, length + radius - head);
    // `addSimpleCylinder` takes the start's radius as `radiusBottom` (and caps it)
    addSimpleCylinder(state, Vec3.scale(Vec3(), unit, start), base, shaft);
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
  const from = length * Inset;
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
