/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Zhenyu Zhang <jump2cn@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Herman Bergwerf <post@hbergwerf.nl>
 */

import { Vec3 } from '../../../../mol-math/linear-algebra';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { Mesh } from '../../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../../mol-geo/geometry/mesh/mesh-builder';
import { CylinderProps } from '../../../../mol-geo/primitive/cylinder';
import { addFixedCountDashedCylinder, addCylinder, addDoubleCylinder } from '../../../../mol-geo/geometry/mesh/builder/cylinder';
import { VisualContext } from '../../../visual';
import { BaseGeometry } from '../../../../mol-geo/geometry/base';
import { Lines } from '../../../../mol-geo/geometry/lines/lines';
import { LinesBuilder } from '../../../../mol-geo/geometry/lines/lines-builder';
import { Cylinders } from '../../../../mol-geo/geometry/cylinders/cylinders';
import { CylindersBuilder } from '../../../../mol-geo/geometry/cylinders/cylinders-builder';
import { Sphere3D } from '../../../../mol-math/geometry/primitives/sphere3d';

export const LinkCylinderParams = {
    linkScale: PD.Numeric(0.45, { min: 0, max: 1, step: 0.01 }),
    linkSpacing: PD.Numeric(1, { min: 0, max: 2, step: 0.01 }),
    linkCap: PD.Boolean(false),
    aromaticScale: PD.Numeric(0.3, { min: 0, max: 1, step: 0.01 }),
    aromaticSpacing: PD.Numeric(1.5, { min: 0, max: 3, step: 0.01 }),
    aromaticDashCount: PD.Numeric(2, { min: 1, max: 6, step: 1 }),
    dashCount: PD.Numeric(4, { min: 0, max: 10, step: 1 }),
    dashScale: PD.Numeric(0.8, { min: 0, max: 2, step: 0.1 }),
    dashCap: PD.Boolean(true),
    stubCap: PD.Boolean(true),
    radialSegments: PD.Numeric(16, { min: 2, max: 56, step: 2 }, BaseGeometry.CustomQualityParamInfo),
    colorMode: PD.Select('default', PD.arrayToOptions(['default', 'interpolate'] as const), BaseGeometry.ShadingCategory)
};
export const DefaultLinkCylinderProps = PD.getDefaultValues(LinkCylinderParams);
export type LinkCylinderProps = typeof DefaultLinkCylinderProps

export const LinkLineParams = {
    linkScale: PD.Numeric(0.5, { min: 0, max: 1, step: 0.1 }),
    linkSpacing: PD.Numeric(0.1, { min: 0, max: 2, step: 0.01 }),
    aromaticDashCount: PD.Numeric(2, { min: 1, max: 6, step: 1 }),
    dashCount: PD.Numeric(4, { min: 0, max: 10, step: 1 }),
};
export const DefaultLinkLineProps = PD.getDefaultValues(LinkLineParams);
export type LinkLineProps = typeof DefaultLinkLineProps

const tmpV12 = Vec3();
const tmpShiftV12 = Vec3();
const tmpShiftV13 = Vec3();
const up = Vec3.create(0, 1, 0);

/** Calculate 'shift' direction that is perpendiculat to v1 - v2 and goes through v3 */
export function calculateShiftDir(out: Vec3, v1: Vec3, v2: Vec3, v3: Vec3 | null) {
    Vec3.normalize(tmpShiftV12, Vec3.sub(tmpShiftV12, v1, v2));
    if (v3 !== null) {
        Vec3.sub(tmpShiftV13, v1, v3);
    } else {
        Vec3.copy(tmpShiftV13, v1); // no reference point, use v1
    }
    Vec3.normalize(tmpShiftV13, tmpShiftV13);

    // ensure v13 and v12 are not colinear
    let dp = Vec3.dot(tmpShiftV12, tmpShiftV13);
    if (1 - Math.abs(dp) < 1e-5) {
        Vec3.set(tmpShiftV13, 1, 0, 0);
        dp = Vec3.dot(tmpShiftV12, tmpShiftV13);
        if (1 - Math.abs(dp) < 1e-5) {
            Vec3.set(tmpShiftV13, 0, 1, 0);
            dp = Vec3.dot(tmpShiftV12, tmpShiftV13);
        }
    }

    Vec3.setMagnitude(tmpShiftV12, tmpShiftV12, dp);
    Vec3.sub(tmpShiftV13, tmpShiftV13, tmpShiftV12);
    return Vec3.normalize(out, tmpShiftV13);
}

export interface LinkBuilderProps {
    linkCount: number
    position: (posA: Vec3, posB: Vec3, edgeIndex: number, adjust: boolean) => void
    radius: (edgeIndex: number) => number,
    referencePosition?: (edgeIndex: number) => Vec3 | null
    style?: (edgeIndex: number) => LinkStyle
    ignore?: (edgeIndex: number) => boolean
    stub?: (edgeIndex: number) => boolean
}

export const EmptyLinkBuilderProps: LinkBuilderProps = {
    linkCount: 0,
    position: () => { },
    radius: () => 0
};

export const enum LinkStyle {
    Solid = 0,
    Dashed = 1,
    Double = 2,
    OffsetDouble = 3,
    Triple = 4,
    OffsetTriple = 5,
    Disk = 6,
    Aromatic = 7,
    MirroredAromatic = 8,
}

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3scale = Vec3.scale;
const v3add = Vec3.add;
const v3sub = Vec3.sub;
const v3setMagnitude = Vec3.setMagnitude;
const v3dot = Vec3.dot;

/**
 * Each edge is included twice to allow for coloring/picking
 * the half closer to the first vertex, i.e. vertex a.
 */
export function createLinkCylinderMesh(ctx: VisualContext, linkBuilder: LinkBuilderProps, props: LinkCylinderProps, mesh?: Mesh): { mesh: Mesh, boundingSphere?: Sphere3D } {
    const { linkCount, referencePosition, position, style, radius, ignore, stub } = linkBuilder;

    if (!linkCount) return { mesh: Mesh.createEmpty(mesh) };

    const { radialSegments, stubCap } = props;

    const vertexCountEstimate = radialSegments * 2 * linkCount * 2;
    const builderState = MeshBuilder.createState(vertexCountEstimate, vertexCountEstimate / 4, mesh);

    const center = Vec3();
    let count = 0;

    let edgeIndex = 0;
    const addOptions: AddLinkOptions = {
        builderState,
        props,
        assignNonAdjustedPosition: (posA, posB) => {
            position(posA, posB, edgeIndex, false);
        },
        referencePosition: referencePosition ? () => referencePosition(edgeIndex) : undefined,
    };
    const addParams: AddLinkParams = {
        a: Vec3(),
        b: Vec3(),
        group: 0,
        linkStub: false,
        linkStyle: LinkStyle.Solid,
        linkRadius: 0,
    };

    for (let _eI = linkCount; edgeIndex < _eI; ++edgeIndex) {
        if (ignore && ignore(edgeIndex)) continue;

        position(addParams.a, addParams.b, edgeIndex, true);

        v3add(center, center, addParams.a);
        v3add(center, center, addParams.b);
        count += 2;

        addParams.group = edgeIndex;
        addParams.linkStub = stubCap && (stub ? stub(edgeIndex) : false);
        addParams.linkRadius = radius(edgeIndex);
        addParams.linkStyle = style ? style(edgeIndex) : LinkStyle.Solid;
        addLinkCylinderMesh(addOptions, addParams);
    }

    const m = MeshBuilder.getMesh(builderState);
    if (count === 0) return { mesh: m };

    // re-use boundingSphere if it has not changed much
    Vec3.scale(center, center, 1 / count);
    const oldBoundingSphere = mesh ? Sphere3D.clone(mesh.boundingSphere) : undefined;
    if (oldBoundingSphere && Vec3.distance(center, oldBoundingSphere.center) / oldBoundingSphere.radius < 0.1) {
        return { mesh: m, boundingSphere: oldBoundingSphere };
    } else {
        return { mesh: m };
    }
}

export interface AddLinkOptions {
    builderState: MeshBuilder.State,
    props: LinkCylinderProps,
    assignNonAdjustedPosition?: (posA: Vec3, posB: Vec3) => void,
    referencePosition?: () => Vec3 | null,
};

export interface AddLinkParams {
    linkStub: boolean,
    linkRadius: number,
    linkStyle: LinkStyle,
    a: Vec3,
    b: Vec3,
    group: number,
}

const AddLinkCylinderState = {
    cylinderProps: {
        radiusTop: 1,
        radiusBottom: 1,
        radialSegments: 8,
        topCap: false,
        bottomCap: false,
    } as CylinderProps,
    vShift: Vec3(),
    va: Vec3(),
    vb: Vec3(),
};

export function addLinkCylinderMesh(options: AddLinkOptions, params: AddLinkParams) {
    const { builderState, props, assignNonAdjustedPosition, referencePosition } = options;
    const { linkRadius, linkStub, linkStyle, a, b } = params;
    const { linkScale, linkSpacing, radialSegments, linkCap, aromaticScale, aromaticSpacing, aromaticDashCount, dashCount, dashScale, dashCap } = props;

    const { cylinderProps, vShift, va, vb } = AddLinkCylinderState;
    cylinderProps.radialSegments = radialSegments;

    Vec3.copy(va, a);
    Vec3.copy(vb, b);

    v3sub(tmpV12, vb, va);
    const dirFlag = v3dot(tmpV12, up) > 0;

    const [topCap, bottomCap] = dirFlag ? [linkStub, linkCap] : [linkCap, linkStub];
    builderState.currentGroup = params.group;

    if (linkStyle === LinkStyle.Solid) {
        cylinderProps.radiusTop = cylinderProps.radiusBottom = linkRadius;
        cylinderProps.topCap = topCap;
        cylinderProps.bottomCap = bottomCap;

        addCylinder(builderState, va, vb, 0.5, cylinderProps);
    } else if (linkStyle === LinkStyle.Dashed) {
        cylinderProps.radiusTop = cylinderProps.radiusBottom = linkRadius * dashScale;
        cylinderProps.topCap = cylinderProps.bottomCap = dashCap;
        addFixedCountDashedCylinder(builderState, va, vb, 0.5, dashCount, linkStub, cylinderProps);
    } else if (linkStyle === LinkStyle.Double || linkStyle === LinkStyle.OffsetDouble || linkStyle === LinkStyle.Triple || linkStyle === LinkStyle.OffsetTriple || linkStyle === LinkStyle.Aromatic || linkStyle === LinkStyle.MirroredAromatic) {
        const order = (linkStyle === LinkStyle.Double || linkStyle === LinkStyle.OffsetDouble) ? 2 :
            (linkStyle === LinkStyle.Triple || linkStyle === LinkStyle.OffsetTriple) ? 3 : 1.5;
        const multiRadius = linkRadius * (linkScale / (0.5 * order));
        const absOffset = (linkRadius - multiRadius) * linkSpacing;

        calculateShiftDir(vShift, va, vb, referencePosition ? referencePosition() : null);

        cylinderProps.topCap = topCap;
        cylinderProps.bottomCap = bottomCap;

        if (linkStyle === LinkStyle.Aromatic || linkStyle === LinkStyle.MirroredAromatic) {
            cylinderProps.radiusTop = cylinderProps.radiusBottom = linkRadius;
            addCylinder(builderState, va, vb, 0.5, cylinderProps);

            const aromaticOffset = linkRadius + aromaticScale * linkRadius + aromaticScale * linkRadius * aromaticSpacing;

            if (assignNonAdjustedPosition) {
                assignNonAdjustedPosition(va, vb);
            }

            v3setMagnitude(tmpV12, v3sub(tmpV12, vb, va), linkRadius * 0.5);
            v3add(va, va, tmpV12);
            v3sub(vb, vb, tmpV12);

            cylinderProps.radiusTop = cylinderProps.radiusBottom = linkRadius * aromaticScale;
            cylinderProps.topCap = cylinderProps.bottomCap = dashCap;
            v3setMagnitude(vShift, vShift, aromaticOffset);
            v3sub(va, va, vShift);
            v3sub(vb, vb, vShift);
            addFixedCountDashedCylinder(builderState, va, vb, 0.5, aromaticDashCount, linkStub, cylinderProps);

            if (linkStyle === LinkStyle.MirroredAromatic) {
                v3setMagnitude(vShift, vShift, aromaticOffset * 2);
                v3add(va, va, vShift);
                v3add(vb, vb, vShift);
                addFixedCountDashedCylinder(builderState, va, vb, 0.5, aromaticDashCount, linkStub, cylinderProps);
            }
        } else if (linkStyle === LinkStyle.OffsetDouble || linkStyle === LinkStyle.OffsetTriple) {
            const multipleOffset = linkRadius + multiRadius + linkScale * linkRadius * linkSpacing;
            v3setMagnitude(vShift, vShift, multipleOffset);

            cylinderProps.radiusTop = cylinderProps.radiusBottom = linkRadius;
            addCylinder(builderState, va, vb, 0.5, cylinderProps);

            v3scale(tmpV12, tmpV12, linkSpacing * linkScale * 0.2);
            v3add(va, va, tmpV12);
            v3sub(vb, vb, tmpV12);

            cylinderProps.radiusTop = cylinderProps.radiusBottom = multiRadius;
            cylinderProps.topCap = dirFlag ? linkStub : dashCap;
            cylinderProps.bottomCap = dirFlag ? dashCap : linkStub;
            v3setMagnitude(vShift, vShift, multipleOffset);
            v3sub(va, va, vShift);
            v3sub(vb, vb, vShift);
            addCylinder(builderState, va, vb, 0.5, cylinderProps);

            if (order === 3) {
                v3setMagnitude(vShift, vShift, multipleOffset * 2);
                v3add(va, va, vShift);
                v3add(vb, vb, vShift);
                addCylinder(builderState, va, vb, 0.5, cylinderProps);
            }
        } else {
            v3setMagnitude(vShift, vShift, absOffset);

            cylinderProps.radiusTop = cylinderProps.radiusBottom = multiRadius;
            if (order === 3) addCylinder(builderState, va, vb, 0.5, cylinderProps);
            addDoubleCylinder(builderState, va, vb, 0.5, vShift, cylinderProps);
        }
    } else if (linkStyle === LinkStyle.Disk) {
        v3scale(tmpV12, tmpV12, 0.475);
        v3add(va, va, tmpV12);
        v3sub(vb, vb, tmpV12);

        cylinderProps.radiusTop = cylinderProps.radiusBottom = linkRadius;
        cylinderProps.topCap = topCap;
        cylinderProps.bottomCap = bottomCap;

        addCylinder(builderState, va, vb, 0.5, cylinderProps);
    }
};

/**
 * Each edge is included twice to allow for coloring/picking
 * the half closer to the first vertex, i.e. vertex a.
 */
export function createLinkCylinderImpostors(ctx: VisualContext, linkBuilder: LinkBuilderProps, props: LinkCylinderProps, cylinders?: Cylinders): { cylinders: Cylinders, boundingSphere?: Sphere3D } {
    const { linkCount, referencePosition, position, style, radius, ignore, stub } = linkBuilder;

    if (!linkCount) return { cylinders: Cylinders.createEmpty(cylinders) };

    const { linkScale, linkSpacing, linkCap, aromaticScale, aromaticSpacing, aromaticDashCount, dashCount, dashScale, dashCap, stubCap, colorMode } = props;
    const interpolate = colorMode === 'interpolate';
    const colorModeFlag = interpolate === true ? 3 : 2;

    const cylindersCountEstimate = linkCount * 2;
    const builder = CylindersBuilder.create(cylindersCountEstimate, cylindersCountEstimate / 4, cylinders);

    const va = Vec3();
    const vb = Vec3();
    const vm = Vec3();
    const vShift = Vec3();

    const center = Vec3();
    let count = 0;

    for (let edgeIndex = 0, _eI = linkCount; edgeIndex < _eI; ++edgeIndex) {
        if (ignore && ignore(edgeIndex)) continue;

        position(va, vb, edgeIndex, true);

        v3add(center, center, va);
        v3add(center, center, vb);
        count += 2;

        const linkRadius = radius(edgeIndex);
        const linkStyle = style ? style(edgeIndex) : LinkStyle.Solid;
        const linkStub = stubCap && (stub ? stub(edgeIndex) : false);

        if (linkStyle === LinkStyle.Solid) {
            v3scale(vm, v3add(vm, va, vb), 0.5);
            builder.add(va[0], va[1], va[2], vm[0], vm[1], vm[2], 1, linkCap, linkStub, colorModeFlag, edgeIndex);
        } else if (linkStyle === LinkStyle.Dashed) {
            v3scale(vm, v3add(vm, va, vb), 0.5);
            builder.addFixedCountDashes(va, vm, dashCount, dashScale, dashCap, dashCap, linkStub, interpolate, edgeIndex);
        } else if (linkStyle === LinkStyle.Double || linkStyle === LinkStyle.OffsetDouble || linkStyle === LinkStyle.Triple || linkStyle === LinkStyle.OffsetTriple || linkStyle === LinkStyle.Aromatic || linkStyle === LinkStyle.MirroredAromatic) {
            const order = (linkStyle === LinkStyle.Double || linkStyle === LinkStyle.OffsetDouble) ? 2 :
                (linkStyle === LinkStyle.Triple || linkStyle === LinkStyle.OffsetTriple) ? 3 : 1.5;
            const multiScale = linkScale / (0.5 * order);
            const absOffset = (linkRadius - multiScale * linkRadius) * linkSpacing;

            v3scale(vm, v3add(vm, va, vb), 0.5);
            calculateShiftDir(vShift, va, vb, referencePosition ? referencePosition(edgeIndex) : null);

            if (linkStyle === LinkStyle.Aromatic || linkStyle === LinkStyle.MirroredAromatic) {
                builder.add(va[0], va[1], va[2], vm[0], vm[1], vm[2], 1, linkCap, linkStub, colorModeFlag, edgeIndex);

                const aromaticOffset = linkRadius + aromaticScale * linkRadius + aromaticScale * linkRadius * aromaticSpacing;

                position(va, vb, edgeIndex, false);

                v3setMagnitude(tmpV12, v3sub(tmpV12, vb, va), linkRadius * 0.5);
                v3add(va, va, tmpV12);

                v3setMagnitude(vShift, vShift, aromaticOffset);
                v3sub(va, va, vShift);
                v3sub(vm, vm, vShift);
                builder.addFixedCountDashes(va, vm, aromaticDashCount, aromaticScale, dashCap, dashCap, linkStub, interpolate, edgeIndex);

                if (linkStyle === LinkStyle.MirroredAromatic) {
                    v3setMagnitude(vShift, vShift, aromaticOffset * 2);
                    v3add(va, va, vShift);
                    v3add(vm, vm, vShift);
                    builder.addFixedCountDashes(va, vm, aromaticDashCount, aromaticScale, dashCap, dashCap, linkStub, interpolate, edgeIndex);
                }
            } else if (linkStyle === LinkStyle.OffsetDouble || linkStyle === LinkStyle.OffsetTriple) {
                const multipleOffset = linkRadius + multiScale * linkRadius + linkScale * linkRadius * linkSpacing;
                v3setMagnitude(vShift, vShift, multipleOffset);

                builder.add(va[0], va[1], va[2], vm[0], vm[1], vm[2], 1, linkCap, linkStub, colorModeFlag, edgeIndex);

                v3setMagnitude(tmpV12, v3sub(tmpV12, va, vm), linkRadius / 1.5);
                v3sub(va, va, tmpV12);

                if (order === 3) builder.add(va[0] + vShift[0], va[1] + vShift[1], va[2] + vShift[2], vm[0] + vShift[0], vm[1] + vShift[1], vm[2] + vShift[2], multiScale, linkCap, linkStub, colorModeFlag, edgeIndex);
                builder.add(va[0] - vShift[0], va[1] - vShift[1], va[2] - vShift[2], vm[0] - vShift[0], vm[1] - vShift[1], vm[2] - vShift[2], multiScale, dashCap, linkStub, colorModeFlag, edgeIndex);
            } else {
                v3setMagnitude(vShift, vShift, absOffset);

                if (order === 3) builder.add(va[0], va[1], va[2], vm[0], vm[1], vm[2], multiScale, linkCap, linkStub, colorModeFlag, edgeIndex);
                builder.add(va[0] + vShift[0], va[1] + vShift[1], va[2] + vShift[2], vm[0] + vShift[0], vm[1] + vShift[1], vm[2] + vShift[2], multiScale, linkCap, linkStub, colorModeFlag, edgeIndex);
                builder.add(va[0] - vShift[0], va[1] - vShift[1], va[2] - vShift[2], vm[0] - vShift[0], vm[1] - vShift[1], vm[2] - vShift[2], multiScale, linkCap, linkStub, colorModeFlag, edgeIndex);
            }
        } else if (linkStyle === LinkStyle.Disk) {
            v3scale(tmpV12, v3sub(tmpV12, vm, va), 0.475);
            v3add(va, va, tmpV12);
            v3sub(vm, vm, tmpV12);

            builder.add(va[0], va[1], va[2], vm[0], vm[1], vm[2], 1, linkCap, linkStub, colorModeFlag, edgeIndex);
        }
    }

    const c = builder.getCylinders();
    if (count === 0) return { cylinders: c };

    // re-use boundingSphere if it has not changed much
    Vec3.scale(center, center, 1 / count);
    const oldBoundingSphere = cylinders ? Sphere3D.clone(cylinders.boundingSphere) : undefined;
    if (oldBoundingSphere && Vec3.distance(center, oldBoundingSphere.center) / oldBoundingSphere.radius < 0.1) {
        return { cylinders: c, boundingSphere: oldBoundingSphere };
    } else {
        return { cylinders: c };
    }
}

/**
 * Each edge is included twice to allow for coloring/picking
 * the half closer to the first vertex, i.e. vertex a.
 */
export function createLinkLines(ctx: VisualContext, linkBuilder: LinkBuilderProps, props: LinkLineProps, lines?: Lines): { lines: Lines, boundingSphere?: Sphere3D } {
    const { linkCount, referencePosition, position, style, ignore } = linkBuilder;

    if (!linkCount) return { lines: Lines.createEmpty(lines) };

    const { linkScale, linkSpacing, aromaticDashCount, dashCount } = props;

    const linesCountEstimate = linkCount * 2;
    const builder = LinesBuilder.create(linesCountEstimate, linesCountEstimate / 4, lines);

    const va = Vec3();
    const vb = Vec3();
    const vm = Vec3();
    const vShift = Vec3();

    const center = Vec3();
    let count = 0;

    const aromaticOffsetFactor = 4.5;
    const multipleOffsetFactor = 3;

    for (let edgeIndex = 0, _eI = linkCount; edgeIndex < _eI; ++edgeIndex) {
        if (ignore && ignore(edgeIndex)) continue;

        position(va, vb, edgeIndex, true);

        v3add(center, center, va);
        v3add(center, center, vb);
        count += 2;

        const linkStyle = style ? style(edgeIndex) : LinkStyle.Solid;

        if (linkStyle === LinkStyle.Solid) {
            v3scale(vm, v3add(vm, va, vb), 0.5);
            builder.add(va[0], va[1], va[2], vm[0], vm[1], vm[2], edgeIndex);
        } else if (linkStyle === LinkStyle.Dashed) {
            v3scale(vm, v3add(vm, va, vb), 0.5);
            builder.addFixedCountDashes(va, vm, dashCount, edgeIndex);
        } else if (linkStyle === LinkStyle.Double || linkStyle === LinkStyle.OffsetDouble || linkStyle === LinkStyle.Triple || linkStyle === LinkStyle.OffsetTriple || linkStyle === LinkStyle.Aromatic || linkStyle === LinkStyle.MirroredAromatic) {
            const order = linkStyle === LinkStyle.Double || linkStyle === LinkStyle.OffsetDouble ? 2 :
                linkStyle === LinkStyle.Triple || linkStyle === LinkStyle.OffsetTriple ? 3 : 1.5;
            const multiRadius = 1 * (linkScale / (0.5 * order));
            const absOffset = (1 - multiRadius) * linkSpacing;

            v3scale(vm, v3add(vm, va, vb), 0.5);
            calculateShiftDir(vShift, va, vb, referencePosition ? referencePosition(edgeIndex) : null);

            if (linkStyle === LinkStyle.Aromatic || linkStyle === LinkStyle.MirroredAromatic) {
                builder.add(va[0], va[1], va[2], vm[0], vm[1], vm[2], edgeIndex);

                const aromaticOffset = absOffset * aromaticOffsetFactor;

                v3setMagnitude(tmpV12, v3sub(tmpV12, vb, va), aromaticOffset * 0.5);
                v3add(va, va, tmpV12);

                v3setMagnitude(vShift, vShift, aromaticOffset);
                v3sub(va, va, vShift);
                v3sub(vm, vm, vShift);
                builder.addFixedCountDashes(va, vm, aromaticDashCount, edgeIndex);

                if (linkStyle === LinkStyle.MirroredAromatic) {
                    v3setMagnitude(vShift, vShift, aromaticOffset * 2);
                    v3add(va, va, vShift);
                    v3add(vm, vm, vShift);
                    builder.addFixedCountDashes(va, vm, aromaticDashCount, edgeIndex);
                }
            } else if (linkStyle === LinkStyle.OffsetDouble || linkStyle === LinkStyle.OffsetTriple) {
                v3setMagnitude(vShift, vShift, absOffset * multipleOffsetFactor);

                builder.add(va[0], va[1], va[2], vm[0], vm[1], vm[2], edgeIndex);

                v3scale(tmpV12, v3sub(tmpV12, va, vm), linkSpacing * linkScale);
                v3sub(va, va, tmpV12);

                if (order === 3) builder.add(va[0] + vShift[0], va[1] + vShift[1], va[2] + vShift[2], vm[0] + vShift[0], vm[1] + vShift[1], vm[2] + vShift[2], edgeIndex);
                builder.add(va[0] - vShift[0], va[1] - vShift[1], va[2] - vShift[2], vm[0] - vShift[0], vm[1] - vShift[1], vm[2] - vShift[2], edgeIndex);
            } else {
                v3setMagnitude(vShift, vShift, absOffset * 1.5);

                if (order === 3) builder.add(va[0], va[1], va[2], vm[0], vm[1], vm[2], edgeIndex);
                builder.add(va[0] + vShift[0], va[1] + vShift[1], va[2] + vShift[2], vm[0] + vShift[0], vm[1] + vShift[1], vm[2] + vShift[2], edgeIndex);
                builder.add(va[0] - vShift[0], va[1] - vShift[1], va[2] - vShift[2], vm[0] - vShift[0], vm[1] - vShift[1], vm[2] - vShift[2], edgeIndex);
            }
        } else if (linkStyle === LinkStyle.Disk) {
            v3scale(tmpV12, v3sub(tmpV12, vm, va), 0.475);
            v3add(va, va, tmpV12);
            v3sub(vm, vm, tmpV12);

            // TODO what to do here? Line as disk doesn't work well.
            builder.add(va[0], va[1], va[2], vm[0], vm[1], vm[2], edgeIndex);
        }
    }

    const l = builder.getLines();
    if (count === 0) return { lines: l };

    // re-use boundingSphere if it has not changed much
    Vec3.scale(center, center, 1 / count);
    const oldBoundingSphere = lines ? Sphere3D.clone(lines.boundingSphere) : undefined;
    if (oldBoundingSphere && Vec3.distance(center, oldBoundingSphere.center) / oldBoundingSphere.radius < 0.1) {
        return { lines: l, boundingSphere: oldBoundingSphere };
    } else {
        return { lines: l };
    }
}
