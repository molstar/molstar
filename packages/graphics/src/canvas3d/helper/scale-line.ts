/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Representation } from '@molstar/graphics/repr/representation';
import { ValueCell } from '@molstar/core/util/value-cell';
import { Geometry } from '@molstar/graphics/geo/geometry/geometry';
import { Lines } from '@molstar/graphics/geo/geometry/lines/lines';
import { LinesBuilder } from '@molstar/graphics/geo/geometry/lines/lines-builder';
import { Text } from '@molstar/graphics/geo/geometry/text/text';
import { TextBuilder } from '@molstar/graphics/geo/geometry/text/text-builder';
import * as ShapeUtils from '@molstar/graphics/geo/shape/shape';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import type { TextValues } from '@molstar/graphics/gl/renderable/text';
import { Visual } from '@molstar/graphics/repr/visual';
import { Camera, type ICamera } from '../camera.js';
import { Viewport } from '../camera/util.js';
import type { ScreenHelper } from './screen-helper.js';

const tmpFrame = Viewport.zero();

export const ScaleLineParams = {
  location: PD.Select('bottom-middle', PD.arrayToOptions(Representation.ScreenEdgeAnchors)),
  locationOffsetX: PD.Numeric(0),
  locationOffsetY: PD.Numeric(0),
  maxWidth: PD.Numeric(120, { min: 40, max: 400, step: 10 }, { description: 'The longest the line is (CSS pixels)' }),
  color: PD.Color(ColorNames.grey),
  lineWidth: PD.Numeric(2, { min: 1, max: 6, step: 0.5 }, { description: 'CSS pixels' }),
  labelSize: PD.Numeric(20, { min: 8, max: 32, step: 1 }, { description: 'CSS pixels' }),
  labelOutline: PD.Numeric(
    0.2,
    { min: 0, max: 0.5, step: 0.01 },
    { description: "The label's outline, of its glyphs' size (0: none): readable over what has its color" },
  ),
  labelOutlineColor: PD.Color(ColorNames.white),
};
export type ScaleLineProps = PD.Values<typeof ScaleLineParams>;

/** From the corner, inward (CSS pixels): the line's end at the corner */
const Margin = 16;
/** The ends' ticks, up from the line (CSS pixels) */
const TickHeight = 6;

/**
 * The scale line (`camera.helper.scale`): a line on the screen layer (`ScreenHelper`, upright) as long as the largest
 * of 1, 2, 5 × 10ⁿ Å that fits `maxWidth` at the camera's target depth (exact orthographic, at the focus in
 * perspective), its length above it. A line of unit length stretched along x by its transform (lines are as wide on
 * screen at any length), the label made again only when its value changes.
 */
export class ScaleLine {
  props: { name: 'on'; params: ScaleLineProps } | { name: 'off'; params: {} } = { name: 'off', params: {} };

  private lines: GraphicsRenderObject | undefined;
  private label: { readonly renderObject: GraphicsRenderObject; readonly text: Text } | undefined;
  /** The value shown (Å) and the line's length (CSS pixels) */
  private value = 0;
  private length = 0;

  constructor(private readonly screen: ScreenHelper) {}

  get isEnabled() {
    return this.props.name === 'on';
  }

  /** Made here, not in a frame: a frame draws only what was compiled before it (`update` changes values only) */
  setProps(props: ScaleLine['props'] | undefined) {
    // none: off (props and states from before the scale line, e.g. `camera.helper` with `axes` only)
    props ??= { name: 'off', params: {} };
    this.props = props;
    this.lines = this.label = undefined;
    this.value = this.length = 0;
    if (props.name !== 'on') {
      this.screen.release(this);
      return;
    }
    this.lines = createLinesRenderObject(props.params);
    this.label = createLabel(props.params);
    const scene = this.scene();
    scene.clear();
    scene.add(this.lines);
    scene.add(this.label.renderObject);
    scene.commit();
  }

  /** The line's length and label for `camera`, in place (the objects' values): before a frame's passes */
  update(camera: ICamera, pixelRatio: number) {
    if (this.props.name !== 'on' || !this.lines || !this.label) return;
    const params = this.props.params;
    const perPixel = worldPerPixel(camera, pixelRatio);
    if (!(perPixel > 0) || !Number.isFinite(perPixel)) return;
    const value = niceValue(params.maxWidth * perPixel);
    const length = value / perPixel;
    if (value === this.value && Math.abs(length - this.length) < 0.01) return;

    // from its anchor toward the screen's middle (left or right edge), or centered on it (an edge's middle)
    const side = Representation.screenAnchorSides(params.location).x;
    const start = side < 0 ? 0 : side > 0 ? -length : -length / 2;
    const line = Mat4.fromTranslation(Mat4(), Vec3.create(start, 0, 0));
    Visual.setTransform(this.lines, Mat4.mul(line, line, Mat4.fromScaling(Mat4(), Vec3.create(length, 1, 1))));
    if (value !== this.value) setLabel(this.label, params, formatLength(value));
    const at = Vec3.create(start + length / 2, TickHeight + 3, 0);
    Visual.setTransform(this.label.renderObject, Mat4.fromTranslation(Mat4(), at));
    this.value = value;
    this.length = length;
    this.scene().update(void 0, false);
  }

  private scene() {
    const params = this.props.params as ScaleLineProps;
    return this.screen.getScene(this, {
      anchor: params.location,
      offsetX: params.locationOffsetX,
      offsetY: params.locationOffsetY,
      followRotation: false,
      inset: Margin,
    });
  }
}

/** The world size of a CSS pixel at the camera's target */
export function worldPerPixel(camera: ICamera, pixelRatio: number) {
  const { position, target, fov } = camera.state;
  // the frame the field of view spans, not just its drawn part
  const cssHeight = Camera.getFrame(tmpFrame, camera).height / pixelRatio;
  return (2 * Math.tan(fov / 2) * Vec3.distance(position, target)) / cssHeight;
}

/** The largest of 1, 2, 5 × 10ⁿ at most `max` */
export function niceValue(max: number) {
  const exponent = Math.floor(Math.log10(max));
  const base = 10 ** exponent;
  const f = max / base;
  return (f >= 5 ? 5 : f >= 2 ? 2 : 1) * base;
}

/** A length in Å as the label shows it: Å below 1000 Å, nm below 10000 Å, µm above */
export function formatLength(angstrom: number) {
  const [value, unit] =
    angstrom < 1000 ? [angstrom, 'Å'] : angstrom < 10000 ? [angstrom / 10, 'nm'] : [angstrom / 1e4, 'µm'];
  // as written, without float noise (0.30000000000000004)
  return `${Number.parseFloat(value.toPrecision(6))} ${unit}`;
}

/** The line, unit long from its origin along x (placed and stretched by its transform), a tick up at each end */
function createLinesRenderObject(props: ScaleLineProps) {
  const builder = LinesBuilder.create(4, 4);
  builder.add(0, 0, 0, 1, 0, 0, 0);
  builder.add(0, 0, 0, 0, TickHeight, 0, 0);
  builder.add(1, 0, 0, 1, TickHeight, 0, 0);
  const lines = builder.getLines();
  lines.setBoundingSphere(Sphere3D.create(Vec3.create(0.5, TickHeight / 2, 0), TickHeight));
  const shape = ShapeUtils.create(
    'scale-line',
    {},
    lines,
    () => props.color,
    () => props.lineWidth,
    () => '',
  );
  return ShapeUtils.createRenderObject(shape, {
    ...PD.getDefaultValues(Lines.Params),
    sizeFactor: 1,
    lineSizeAttenuation: false,
  });
}

const LabelProps = {
  ...PD.getDefaultValues(Text.Params),
  attachment: 'bottom-center' as const,
  cellSize: 0,
  fontWeight: 'bold' as const,
};

/** The label's text props: centered above the line, outlined (`borderWidth`, `borderColor`) */
function labelProps(props: ScaleLineProps) {
  return { ...LabelProps, borderWidth: props.labelOutline, borderColor: props.labelOutlineColor };
}

function buildText(props: ScaleLineProps, label: string, text?: Text) {
  const builder = TextBuilder.create(labelProps(props), 16, 16, text);
  builder.add(label, 0, 0, 0, 0, props.labelSize, 0);
  return builder.getText();
}

/** The label's object, its text a placeholder until the first `update` */
function createLabel(props: ScaleLineProps) {
  const text = buildText(props, formatLength(1));
  const shape = ShapeUtils.create(
    'scale-line-label',
    {},
    text,
    () => props.color,
    () => 1,
    () => '',
  );
  return { renderObject: ShapeUtils.createRenderObject(shape, labelProps(props)), text };
}

/** The label's text built again into its buffers, its counts and bounds after it (as shape representations do) */
function setLabel(
  label: { readonly renderObject: GraphicsRenderObject; readonly text: Text },
  props: ScaleLineProps,
  s: string,
) {
  const text = buildText(props, s, label.text);
  const { values } = label.renderObject;
  ValueCell.updateIfChanged(values.drawCount, Geometry.getDrawCount(text));
  ValueCell.updateIfChanged(values.uVertexCount, Geometry.getVertexCount(text));
  Text.Utils.updateBoundingSphere(values as TextValues, text);
}
