/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualUpdateState } from '../../util';
import { VisualContext } from '../../visual';
import { Structure, StructureElement } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { EPSILON, Mat4, Quat, Vec3 } from '../../../mol-math/linear-algebra';
import { Axes3D, Box3D } from '../../../mol-math/geometry';
import { ComplexImageParams, ComplexImageVisual, ComplexVisual } from '../complex-visual';
import { eachSerialElement, ElementIterator, getSerialElementLoci } from './util/element';
import { Image, InterpolationTypes } from '../../../mol-geo/geometry/image/image';
import { transformPositionArray } from '../../../mol-geo/util';
import { BaseGeometry } from '../../../mol-geo/geometry/base';
import { PositionLocation } from '../../../mol-geo/util/location-iterator';
import { Color } from '../../../mol-util/color/color';
import { clamp } from '../../../mol-math/interpolate';
import { ColorTheme } from '../../../mol-theme/color';
import { packIntToRGBArray } from '../../../mol-util/number-packing';
import { SizeTheme } from '../../../mol-theme/size';
import { Plane3D } from '../../../mol-math/geometry/primitives/plane3d';
import { degToRad } from '../../../mol-math/misc';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3set = Vec3.set;
const v3transformMat4 = Vec3.transformMat4;
const v3squaredDistance = Vec3.squaredDistance;

export const PlaneImageParams = {
    ...ComplexImageParams,
    interpolation: PD.Select('nearest', PD.objectToOptions(InterpolationTypes)),
    imageResolution: PD.Numeric(0.5, { min: 0.01, max: 20, step: 0.01 }, { description: 'Grid resolution/cell spacing.', ...BaseGeometry.CustomQualityParamInfo }),
    mode: PD.Select('frame', PD.arrayToOptions(['frame', 'plane'] as const), { description: 'Frame: slice through the structure along arbitrary axes in any step size. Plane: an arbitrary plane defined by point and normal.' }),
    offset: PD.Numeric(0, { min: -1, max: 1, step: 0.01 }, { isEssential: true, immediateUpdate: true, hideIf: p => p.mode !== 'frame', description: 'Relative offset from center.' }),
    axis: PD.Select('c', PD.arrayToOptions(['a', 'b', 'c'] as const), { isEssential: true, hideIf: p => p.mode !== 'frame' }),
    rotation: PD.Group({
        axis: PD.Vec3(Vec3.create(1, 0, 0), {}, { description: 'Axis of rotation' }),
        angle: PD.Numeric(0, { min: -180, max: 180, step: 1 }, { immediateUpdate: true, description: 'Axis rotation angle in Degrees' }),
    }, { isExpanded: true, hideIf: p => p.mode !== 'frame' }),
    plane: PD.Group({
        point: PD.Vec3(Vec3.create(0, 0, 0), {}, { description: 'Plane point' }),
        normal: PD.Vec3(Vec3.create(1, 0, 0), {}, { description: 'Plane normal' }),
    }, { isExpanded: true, hideIf: p => p.mode !== 'plane' }),
    extent: PD.Select('frame', PD.arrayToOptions(['frame', 'sphere'] as const), { description: 'Extent of the plane, either box (frame) or sphere.' }),
    margin: PD.Numeric(4, { min: 0, max: 50, step: 1 }, { immediateUpdate: true, description: 'Margin around the structure in Angstrom' }),
    frame: PD.Select('principalAxes', PD.arrayToOptions(['principalAxes', 'boundingBox'] as const)),
    antialias: PD.Boolean(true, { description: 'Antialiasing of structure edges.' }),
    cutout: PD.Boolean(false, { description: 'Cutout the structure from the image.' }),
    defaultColor: PD.Color(Color(0xCCCCCC), { description: 'Default color for parts of the image that are not covered by the color theme.' }),
    includeParent: PD.Boolean(false, { description: 'Show parent structure (but within extent of this structure).' }),
};
export type PlaneImageParams = typeof PlaneImageParams

export function PlaneImageVisual(materialId: number): ComplexVisual<PlaneImageParams> {
    return ComplexImageVisual<PlaneImageParams>({
        defaultProps: PD.getDefaultValues(PlaneImageParams),
        createGeometry: createPlaneImage,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<PlaneImageParams>, currentProps: PD.Values<PlaneImageParams>, newTheme: Theme, currentTheme: Theme) => {
            state.createGeometry = (
                newProps.imageResolution !== currentProps.imageResolution ||
                newProps.mode !== currentProps.mode ||
                newProps.margin !== currentProps.margin ||
                newProps.frame !== currentProps.frame ||
                newProps.extent !== currentProps.extent ||
                !Vec3.equals(newProps.rotation.axis, currentProps.rotation.axis) ||
                newProps.rotation.angle !== currentProps.rotation.angle ||
                newProps.offset !== currentProps.offset ||
                newProps.axis !== currentProps.axis ||
                !Vec3.equals(newProps.plane.point, currentProps.plane.point) ||
                !Vec3.equals(newProps.plane.normal, currentProps.plane.normal) ||
                newProps.antialias !== currentProps.antialias ||
                newProps.cutout !== currentProps.cutout ||
                newProps.defaultColor !== currentProps.defaultColor ||
                !ColorTheme.areEqual(newTheme.color, currentTheme.color) ||
                !SizeTheme.areEqual(newTheme.size, currentTheme.size)
            );
        }
    }, materialId);
}

//

export interface PlaneImageProps {
    imageResolution: number,
    mode: 'frame' | 'plane',
    offset: number,
    axis: 'a' | 'b' | 'c',
    margin: number,
    frame: 'principalAxes' | 'boundingBox',
    extent: 'frame' | 'sphere',
    rotation: { axis: Vec3, angle: number },
    plane: { point: Vec3, normal: Vec3 },
    antialias: boolean,
    cutout: boolean,
    defaultColor: Color,
    includeParent: boolean,
}

function getFrame(structure: Structure, props: PlaneImageProps) {
    const { mode, axis, frame, extent, margin, rotation, plane, includeParent } = props;

    if (includeParent && structure.child) {
        structure = structure.child;
    }

    const size = Vec3();
    const scale = Vec3();
    const major = Vec3();
    const minor = Vec3();
    const normal = Vec3();
    const center = Vec3();

    let a = 0, b = 0, c = 0;
    let dirA: Vec3, dirB: Vec3, dirC: Vec3;

    if (frame === 'principalAxes') {
        const axes = Structure.getPrincipalAxes(structure).boxAxes;
        [a, b, c] = Axes3D.size(Vec3(), axes);
        dirA = axes.dirA;
        dirB = axes.dirB;
        dirC = axes.dirC;
        Vec3.copy(center, axes.origin);
    } else {
        [a, b, c] = Box3D.size(Vec3(), structure.boundary.box);
        dirA = Vec3.create(1, 0, 0);
        dirB = Vec3.create(0, 1, 0);
        dirC = Vec3.create(0, 0, 1);
        Vec3.copy(center, structure.boundary.sphere.center);
    }

    Vec3.set(scale, a, b, c);

    if (axis === 'c') {
        Vec3.set(size, a, b, c);
        Vec3.copy(major, dirA);
        Vec3.copy(minor, dirB);
        Vec3.copy(normal, dirC);
    } else if (axis === 'b') {
        Vec3.set(size, a, c, b);
        Vec3.copy(major, dirA);
        Vec3.copy(normal, dirB);
        Vec3.copy(minor, dirC);
    } else {
        Vec3.set(size, b, c, a);
        Vec3.copy(normal, dirA);
        Vec3.copy(major, dirB);
        Vec3.copy(minor, dirC);
    }

    if (rotation.angle !== 0) {
        const ra = Vec3();
        Vec3.scaleAndAdd(ra, ra, dirA, rotation.axis[0]);
        Vec3.scaleAndAdd(ra, ra, dirB, rotation.axis[1]);
        Vec3.scaleAndAdd(ra, ra, dirC, rotation.axis[2]);
        Vec3.normalize(ra, ra);

        const rm = Mat4.fromRotation(Mat4(), degToRad(rotation.angle), ra);
        Vec3.transformDirection(major, major, rm);
        Vec3.transformDirection(minor, minor, rm);
        Vec3.transformDirection(normal, normal, rm);
    }

    if (extent === 'sphere' || rotation.angle !== 0) {
        const r = structure.boundary.sphere.radius * 2;
        const s = Vec3.magnitude(Box3D.size(Vec3(), Box3D.fromSphere3D(Box3D(), structure.boundary.sphere)));
        Vec3.set(size, s, s, r);
        if (extent === 'sphere') {
            Vec3.set(scale, r, r, r);
        }
    }

    Vec3.addScalar(size, size, margin * 2);
    Vec3.addScalar(scale, scale, margin * 2);

    const trimRotation = Quat.identity();
    if (frame === 'principalAxes') {
        Quat.fromBasis(trimRotation,
            Vec3.normalize(Vec3(), dirA),
            Vec3.normalize(Vec3(), dirB),
            Vec3.normalize(Vec3(), dirC)
        );
    }

    if (mode === 'plane') {
        Vec3.copy(center, plane.point);
        Vec3.copy(normal, plane.normal);

        Vec3.cross(major, normal, Vec3.unitX);
        if (Vec3.dot(major, major) < EPSILON) {
            Vec3.cross(major, normal, Vec3.unitY);
        }
        Vec3.normalize(major, major);

        Vec3.cross(minor, normal, major);
        Vec3.normalize(minor, minor);
    }

    const trim: Image.Trim = {
        type: extent === 'sphere' ? 2 : 3,
        center,
        scale,
        rotation: trimRotation,
        transform: Mat4.identity(),
    };

    return { size, major, minor, normal, center, trim };
}

export function createPlaneImage(ctx: VisualContext, structure: Structure, theme: Theme, props: PlaneImageProps, image?: Image): Image {
    const { imageResolution, offset, antialias, cutout, defaultColor } = props;
    const scaleFactor = 1 / imageResolution;

    const color = 'color' in theme.color && theme.color.color
        ? theme.color.color
        : () => Color(0xffffff);

    const { size, major, minor, normal, center, trim } = getFrame(structure, props);

    const scale = Vec3.create(size[0], size[1], 1);
    const offsetDir = Vec3.setMagnitude(Vec3(), normal, size[2] / 2);

    const width = Math.floor(size[1] * scaleFactor);
    const height = Math.floor(size[0] * scaleFactor);

    const m = Mat4.identity();
    const v = Vec3();
    const anchor = Vec3();

    Vec3.add(v, center, major);
    Mat4.targetTo(m, center, v, minor);
    Vec3.scaleAndAdd(anchor, center, offsetDir, offset);
    Mat4.setTranslation(m, anchor);
    Mat4.mul(m, m, Mat4.rotY90);
    Mat4.scale(m, m, scale);

    const { getSerialIndex } = structure.serialMapping;
    const isVertex = theme.color.granularity.startsWith('vertex');

    const plane = Plane3D.fromNormalAndCoplanarPoint(Plane3D(), Vec3.normalize(Vec3(), normal), anchor);
    const invM = Mat4.invert(Mat4(), m);

    const pl = PositionLocation(Vec3(), Vec3());
    const el = StructureElement.Location.create(structure);

    const { units } = structure;
    let maxRadius = 0;
    for (let i = 0, il = units.length; i < il; ++i) {
        const { elements } = units[i];
        el.unit = units[i];
        for (let j = 0, jl = elements.length; j < jl; ++j) {
            el.element = elements[j];
            const r = theme.size.size(el);
            if (r > maxRadius) maxRadius = r;
        }
    }

    const imageArray = new Uint8Array(width * height * 4);
    const groupArray = new Uint8Array(width * height * 4);

    const distArray = new Float32Array(width * height);
    distArray.fill(Number.MAX_VALUE);

    const p = Vec3();
    const pp = Vec3();
    const pn = Vec3();

    if (isVertex) {
        let i = 0;
        for (let ih = 0; ih < height; ++ih) {
            for (let iw = 0; iw < width; ++iw) {
                const y = (clamp(iw + 0.5, 0, width - 1) / width) - 0.5;
                const x = (clamp(ih + 0.5, 0, height - 1) / height) - 0.5;
                Vec3.set(v, x, -y, 0);
                Vec3.transformMat4(v, v, m);

                Vec3.copy(pl.position, v);
                const c = color(pl, false);
                Color.toArray(c, imageArray, i);
                imageArray[i + 3] = cutout ? 0 : 255;

                i += 4;
            }
        }
    } else {
        for (let i = 0, il = width * height * 4; i < il; i += 4) {
            Color.toArray(defaultColor, imageArray, i);
            imageArray[i + 3] = cutout ? 0 : 255;
        }
    }

    for (let i = 0, il = units.length; i < il; ++i) {
        const unit = units[i];
        const { elements, conformation: c } = unit;
        el.unit = units[i];

        for (let j = 0, jl = elements.length; j < jl; ++j) {
            const eI = elements[j];
            el.element = eI;

            c.position(eI, p);
            const dist = Plane3D.distanceToPoint(plane, p);
            if (Math.abs(dist) > maxRadius) continue;

            const r = theme.size.size(el);
            if (Math.abs(dist) > r) continue;

            const rf = Math.cos(Math.abs(dist) / r);
            const tol = antialias ? imageResolution * rf : 0;
            const rTol = r + (tol / 2);
            const rTolSq = rTol * rTol;

            Vec3.scaleAndAdd(pp, p, plane.normal, -dist);
            Vec3.transformMat4(pn, pp, invM);
            Vec3.addScalar(pn, pn, 0.5);

            const x = Math.floor(pn[0] * height);
            const y = width - Math.ceil(pn[1] * width);

            // Number of grid points, round this up...
            const ng = Math.ceil(r * scaleFactor);

            // Extents of grid to consider for this atom
            const begX = Math.max(0, x - ng);
            const begY = Math.max(0, y - ng);

            // Add two to these points:
            // - x,y are floor'd values so this ensures coverage
            // - these are loop limits (exclusive)
            const endX = Math.min(height, x + ng + 2);
            const endY = Math.min(width, y + ng + 2);

            Vec3.copy(pl.position, pp);
            const col = isVertex ? defaultColor : color(el, false);
            const idx = getSerialIndex(el.unit, el.element);

            for (let xi = begX; xi < endX; ++xi) {
                for (let yi = begY; yi < endY; ++yi) {
                    const xx = (clamp(xi + 0.5, 0, height - 1) / height) - 0.5;
                    const yy = (clamp(yi + 0.5, 0, width - 1) / width) - 0.5;
                    v3set(v, xx, -yy, 0);
                    v3transformMat4(v, v, m);
                    const distSq = v3squaredDistance(v, p);
                    if (distSq > rTolSq) continue;

                    const k = xi * width + yi;
                    if (distSq < distArray[k]) {
                        const k4 = k * 4;
                        const d = Math.sqrt(distSq) - r + tol / 2;
                        let f = d > 0 ? 1 - d / tol : 1;

                        if (isVertex) {
                            if (f === 1) {
                                distArray[k] = distSq;
                            } else {
                                if (cutout) {
                                    if (groupArray[k4] !== 0 || groupArray[k4 + 1] !== 0 || groupArray[k4 + 2] !== 0) {
                                        f = 1;
                                    }
                                }
                            }
                        } else {
                            if (f === 1) {
                                distArray[k] = distSq;
                                Color.toArray(col, imageArray, k4);
                            } else {
                                if (cutout) {
                                    Color.toArray(col, imageArray, k4);
                                    if (groupArray[k4] !== 0 || groupArray[k4 + 1] !== 0 || groupArray[k4 + 2] !== 0) {
                                        f = 1;
                                    }
                                } else {
                                    Color.toArray(Color.interpolate(Color.fromArray(imageArray, k4), col, f), imageArray, k4);
                                }
                            }
                        }
                        packIntToRGBArray(idx, groupArray, k4);
                        groupArray[k4 + 3] = antialias ? Math.round(255 * f) : 255;

                        if (cutout) {
                            imageArray[k * 4 + 3] = antialias ? Math.round(255 * f) : 255;
                        }
                    }
                }
            }
        }
    }

    const imageTexture = { width, height, array: imageArray, flipY: true };
    const groupTexture = { width, height, array: groupArray, flipY: true };
    const valueTexture = { width: 1, height: 1, array: new Float32Array(1), flipY: true };

    const corners = new Float32Array([
        -0.5, 0.5, 0,
        0.5, 0.5, 0,
        -0.5, -0.5, 0,
        0.5, -0.5, 0
    ]);
    transformPositionArray(m, corners, 0, 4);

    return Image.create(imageTexture, corners, groupTexture, valueTexture, trim, -1, image);
}
