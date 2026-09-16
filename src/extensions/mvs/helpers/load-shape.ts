/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { Mat4 } from '../../../mol-math/linear-algebra';
import { Color } from '../../../mol-util/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { alphaForNode, clippingForNode, transformFromRotationTranslation } from '../load-helpers';
import { getChildren } from '../tree/generic/tree-schema';
import { MolstarNode, MolstarSubtree } from '../tree/molstar/molstar-tree';
import { DefaultColor } from '../tree/mvs/mvs-tree';
import { decodeColor } from './utils';


/** The parsed-resource kind a `shape` node sits on, taken from the upstream `parse` node. */
export type ShapeFormat = 'vtp' | 'ply' | 'obj'

/** Format-specific extensions to a `shape` node, carried in `custom` rather than as spec params
 * and prefixed by the format they apply to -- the same `<prefix>_<name>` convention as
 * `molstar_mesh_params` and friends. */
interface ShapeCustomProps {
    /** VTP: name of the scalar data array to color by. */
    vtp_attribute?: string,
    /** VTP: whether `vtp_attribute` names a per-point or a per-cell array. Defaults to `'point'`. */
    vtp_attribute_source?: 'point' | 'cell',
    /** VTP: name of a molstar color list to map attribute values through, e.g. `'viridis'`.
     * Custom properties are vendor extensions, so a molstar-internal name is appropriate here.
     * Defaults to the shape provider's own default. */
    vtp_palette?: string,
    /** VTP: `[min, max]` for the color scale. Defaults to the min/max of the attribute values. */
    vtp_domain?: [number, number],
    /** PLY: where colors come from. `'vertex'`/`'material'` read the conventional `red`/`green`/`blue`
     * properties of the respective element. Defaults to `'uniform'` (the `color` children). */
    ply_coloring?: 'uniform' | 'vertex' | 'material',
    /** OBJ: where colors come from. `'vertex'` uses colors embedded in the OBJ; `'custom'` colors
     * each material group via `obj_material_colors`. Defaults to `'uniform'`. (`'given'`, which
     * reads Kd from an MTL file, is not supported -- there is no way to reference the MTL yet.) */
    obj_coloring?: 'uniform' | 'vertex' | 'custom',
    /** OBJ: color per material name, for `obj_coloring: 'custom'`. Material names come from the
     * OBJ's own `usemtl` directives, so no MTL file is needed. Colors have to be listed here
     * because setting `coloring` replaces the whole mapped value -- there is no way to ask the
     * provider for its generated per-material palette from here, and no access to the material
     * names to pre-fill them. Materials not listed render grey (the provider's own fallback). */
    obj_material_colors?: Record<string, string>,
    /** Mesh-level rendering params, as for `primitives`. */
    molstar_mesh_params?: Record<string, any>,
}

/** Color for a `shape` node, taken from its `color` children (last one wins), as `volume` does. */
function shapeColor(node: MolstarSubtree<'shape'>): Color {
    const children = getChildren(node).filter(c => c.kind === 'color') as MolstarNode<'color'>[];
    const color = children.length > 0 ? children[children.length - 1].params.color : DefaultColor;
    return decodeColor(color) ?? decodeColor(DefaultColor)!;
}

/** Create params for the shape representation of a `shape` node.
 *
 * These are the params of the format's own shape provider (`createVtpShapeParams` and friends), so
 * an MVS-loaded shape is the same state object as one opened through the UI and keeps the full
 * interactive coloring controls. Every value MVS controls is supplied explicitly -- in particular
 * `colorTheme`/`coloring` is always set, never left to the provider's file-derived default, so the
 * same MVSJ renders the same way regardless of which arrays the resource happens to contain. */
export function shapeRepresentationProps(node: MolstarSubtree<'shape'>, format: ShapeFormat) {
    const custom = (node.custom ?? {}) as ShapeCustomProps;
    const base = {
        alpha: alphaForNode(node),
        clip: clippingForNode(node),
        // Mesh-level rendering params follow the same custom-prop convention as `primitives`.
        ...custom.molstar_mesh_params,
    };
    const color = shapeColor(node);

    if (format === 'vtp') {
        const colorTheme: PD.NamedParams = custom.vtp_attribute
            ? {
                name: 'attribute',
                params: {
                    name: `${custom.vtp_attribute_source ?? 'point'}:${custom.vtp_attribute}`,
                    // Setting `colorTheme` replaces the whole mapped value, so every param of the
                    // group must be supplied -- an omitted `colors` would be undefined rather than
                    // falling back to the provider's default.
                    colors: custom.vtp_palette ?? 'viridis',
                    domain: custom.vtp_domain
                        ? { name: 'custom', params: custom.vtp_domain }
                        : { name: 'auto', params: { symmetric: false } },
                },
            }
            : { name: 'uniform', params: { color } };
        return { ...base, colorTheme };
    }
    if (format === 'ply') {
        const mode = custom.ply_coloring ?? 'uniform';
        if (mode === 'uniform') {
            // PLY's uniform group also carries saturation/lightness modifiers; MVS pins them to 0
            // so the rendered color is exactly the one the file asked for.
            return { ...base, coloring: { name: 'uniform', params: { color, saturation: 0, lightness: 0 } } };
        }
        // The property names are part of the provider's params, and PLY convention is red/green/blue.
        // Pinning them keeps the result independent of what the provider happens to auto-detect;
        // a file using other names renders grey rather than picking something arbitrary.
        return { ...base, coloring: { name: mode, params: { red: 'red', green: 'green', blue: 'blue' } } };
    }
    const objMode = custom.obj_coloring ?? 'uniform';
    if (objMode === 'vertex') {
        return { ...base, coloring: { name: 'vertex', params: {} } };
    }
    if (objMode === 'custom') {
        const colors: Record<string, Color> = {};
        for (const [name, c] of Object.entries(custom.obj_material_colors ?? {})) {
            const decoded = decodeColor(c);
            if (decoded !== undefined) colors[name] = decoded;
        }
        return { ...base, coloring: { name: 'custom', params: colors } };
    }
    return { ...base, coloring: { name: 'uniform', params: { color } } };
}

/** Build a single matrix from a `transform`/`instance` node's params. */
function matrixFromTransformNode(params: MolstarNode<'transform'>['params']): Mat4 {
    if (params.matrix) return Mat4.fromArray(Mat4(), params.matrix, 0);
    if (params.rotation_center) {
        // `rotation_center` needs a transformer that can resolve it against the transformed
        // object (`TransformStructureConformation`'s `components` variant). Shapes carry a plain
        // matrix array instead, so there is nothing to resolve it against. Fail loudly rather
        // than quietly rotating about the origin.
        throw new Error(`'rotation_center' is not supported on 'transform'/'instance' nodes under a 'shape' node; bake the centering into 'matrix' instead.`);
    }
    return transformFromRotationTranslation(params.rotation, params.translation);
}

/** Compose the `transform` and `instance` children of a `shape` node into instance matrices.
 * There is no transformer for shapes, so instead of applying a transform node the matrices are
 * folded into the Shape's own `transforms` array: one entry per instance, each premultiplied by
 * the accumulated `transform` matrix. */
export function shapeTransforms(node: MolstarSubtree<'shape'>): Mat4[] | undefined {
    let transform: Mat4 | undefined;
    for (const child of getChildren(node)) {
        if (child.kind !== 'transform') continue;
        const m = matrixFromTransformNode((child as MolstarNode<'transform'>).params);
        transform = transform ? Mat4.mul(Mat4(), m, transform) : m;
    }

    const instances = getChildren(node).filter(c => c.kind === 'instance') as MolstarNode<'instance'>[];
    if (instances.length === 0) return transform ? [transform] : undefined;

    return instances.map(i => {
        const m = matrixFromTransformNode(i.params);
        return transform ? Mat4.mul(Mat4(), m, transform) : m;
    });
}
