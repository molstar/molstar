/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from '../../../../mol-model/structure';
import { Theme } from '../../../../mol-theme/theme';
import { VisualContext } from '../../../visual';
import { WebGLContext } from '../../../../mol-gl/webgl/context';
import { Spheres } from '../../../../mol-geo/geometry/spheres/spheres';
import { Mesh } from '../../../../mol-geo/geometry/mesh/mesh';
import { ElementSphereImpostorProps, ElementSphereMeshProps, createStructureElementSphereImpostor, createStructureElementSphereMesh } from '../../../structure/visual/util/element';
import { checkSphereImpostorSupport } from '../../../structure/visual/util/common';

/** The subset of representation props needed to build a structure target geometry. */
export interface StructureTargetProps {
    sizeFactor: number
    detail: number
    ignoreHydrogens: boolean
    ignoreHydrogensVariant: 'all' | 'non-polar'
    traceOnly: boolean
    stride: number
    tryUseImpostor: boolean
}

export function structureTargetUseImpostor(props: StructureTargetProps, webgl?: WebGLContext): boolean {
    return props.tryUseImpostor && checkSphereImpostorSupport(webgl);
}

/** Build the geometry for a structure target: impostor spheres or a sphere mesh of its atoms. */
export function createStructureTargetGeometry(ctx: VisualContext, structure: Structure, theme: Theme, props: StructureTargetProps, webgl?: WebGLContext, existing?: Spheres | Mesh): Spheres | Mesh | Promise<Spheres | Mesh> {
    if (structureTargetUseImpostor(props, webgl)) {
        const p: ElementSphereImpostorProps = {
            sizeFactor: props.sizeFactor,
            ignoreHydrogens: props.ignoreHydrogens,
            ignoreHydrogensVariant: props.ignoreHydrogensVariant,
            traceOnly: props.traceOnly,
            stride: props.stride,
        };
        return createStructureElementSphereImpostor(ctx, structure, theme, p,
            existing?.kind === 'spheres' ? existing as Spheres : undefined);
    } else {
        const p: ElementSphereMeshProps = {
            detail: props.detail,
            sizeFactor: props.sizeFactor,
            ignoreHydrogens: props.ignoreHydrogens,
            ignoreHydrogensVariant: props.ignoreHydrogensVariant,
            traceOnly: props.traceOnly,
            stride: props.stride,
        };
        return createStructureElementSphereMesh(ctx, structure, theme, p,
            existing?.kind === 'mesh' ? existing as Mesh : undefined);
    }
}

/** Whether the structure geometry must be rebuilt because relevant props changed. */
export function structureTargetGeometryPropsChanged(oldProps: StructureTargetProps, newProps: StructureTargetProps, webgl?: WebGLContext): boolean {
    if (structureTargetUseImpostor(newProps, webgl)) {
        return (
            newProps.ignoreHydrogens !== oldProps.ignoreHydrogens ||
            newProps.ignoreHydrogensVariant !== oldProps.ignoreHydrogensVariant ||
            newProps.traceOnly !== oldProps.traceOnly ||
            newProps.stride !== oldProps.stride
        );
    }
    return (
        newProps.sizeFactor !== oldProps.sizeFactor ||
        newProps.detail !== oldProps.detail ||
        newProps.ignoreHydrogens !== oldProps.ignoreHydrogens ||
        newProps.ignoreHydrogensVariant !== oldProps.ignoreHydrogensVariant ||
        newProps.traceOnly !== oldProps.traceOnly ||
        newProps.stride !== oldProps.stride
    );
}

/** Whether toggling impostor usage requires recreating the render object (geometry kind changes). */
export function structureTargetMustRecreate(oldProps: StructureTargetProps, newProps: StructureTargetProps): boolean {
    return oldProps.tryUseImpostor !== newProps.tryUseImpostor;
}

/** Whether the impostor sphere size factor changed (impostor spheres carry size as a uniform). */
export function structureTargetSizeFactorChanged(oldProps: StructureTargetProps, newProps: StructureTargetProps, webgl?: WebGLContext): boolean {
    return structureTargetUseImpostor(newProps, webgl) && newProps.sizeFactor !== oldProps.sizeFactor;
}
