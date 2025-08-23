/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bool, float, int, literal, nullable, OptionalField, RequiredField } from '../generic/field-schema';
import { SimpleParamsSchema, UnionParamsSchema } from '../generic/params-schema';
import { Matrix, Vector3 } from './param-types';

const Cartoon = {
    /** Scales the corresponding visuals */
    size_factor: OptionalField(float, 1, 'Scales the corresponding visuals.'),
    /** Simplify corkscrew helices to tubes. */
    tubular_helices: OptionalField(bool, false, 'Simplify corkscrew helices to tubes.'),
};

const Backbone = {
    /** Scales the corresponding visuals */
    size_factor: OptionalField(float, 1, 'Scales the corresponding visuals.'),
};

const BallAndStick = {
    /** Scales the corresponding visuals */
    size_factor: OptionalField(float, 1, 'Scales the corresponding visuals.'),
    /** Controls whether hydrogen atoms are drawn. */
    ignore_hydrogens: OptionalField(bool, false, 'Controls whether hydrogen atoms are drawn.'),
};

const Line = {
    /** Scales the corresponding visuals */
    size_factor: OptionalField(float, 1, 'Scales the corresponding visuals.'),
    /** Controls whether hydrogen atoms are drawn. */
    ignore_hydrogens: OptionalField(bool, false, 'Controls whether hydrogen atoms are drawn.'),
};

const Spacefill = {
    /** Scales the corresponding visuals */
    size_factor: OptionalField(float, 1, 'Scales the corresponding visuals.'),
    /** Controls whether hydrogen atoms are drawn. */
    ignore_hydrogens: OptionalField(bool, false, 'Controls whether hydrogen atoms are drawn.'),
};

const Carbohydrate = {
    /** Scales the corresponding visuals */
    size_factor: OptionalField(float, 1, 'Scales the corresponding visuals.'),
};

const Surface = {
    /** Type of surface representation. (Default is 'molecular') */
    surface_type: OptionalField(literal('molecular', 'gaussian'), 'molecular', `Type of surface representation. (Default is 'molecular')`),
    /** Scales the corresponding visuals */
    size_factor: OptionalField(float, 1, 'Scales the corresponding visuals.'),
    /** Controls whether hydrogen atoms are drawn. */
    ignore_hydrogens: OptionalField(bool, false, 'Controls whether hydrogen atoms are drawn.'),
};

export const MVSRepresentationParams = UnionParamsSchema(
    'type',
    'Representation type',
    {
        cartoon: SimpleParamsSchema(Cartoon),
        backbone: SimpleParamsSchema(Backbone),
        ball_and_stick: SimpleParamsSchema(BallAndStick),
        line: SimpleParamsSchema(Line),
        spacefill: SimpleParamsSchema(Spacefill),
        carbohydrate: SimpleParamsSchema(Carbohydrate),
        surface: SimpleParamsSchema(Surface),
    },
);

const VolumeIsoSurface = {
    /** Relative isovalue. */
    relative_isovalue: OptionalField(nullable(float), null, 'Relative isovalue.'),
    /** Absolute isovalue. Overrides `relative_isovalue`. */
    absolute_isovalue: OptionalField(nullable(float), null, 'Absolute isovalue. Overrides `relative_isovalue`.'),
    /** Show mesh wireframe. Defaults to false. */
    show_wireframe: OptionalField(bool, false, 'Show mesh wireframe. Defaults to false.'),
    /** Show mesh faces. Defaults to true. */
    show_faces: OptionalField(bool, true, 'Show mesh faces. Defaults to true.'),
};

const VolumeGridSlice = {
    /** Dimension of the grid slice, i.e. 'x', 'y', or 'z'. */
    dimension: RequiredField(literal('x', 'y', 'z'), 'Dimension of the grid slice, i.e. \'x\', \'y\', or \'z\'.'),
    /** Index of the grid slice in the specified dimension. 0-based index, i.e. 0 is the first slice. */
    absolute_index: OptionalField(nullable(int), null, 'Index of the grid slice in the specified dimension. 0-based index, i.e. 0 is the first slice.'),
    /** Relative index of the grid slice in the specified dimension. 0.0 is the first slice, 1.0 is the last slice. Overrides `absolute_index`. */
    relative_index: OptionalField(nullable(float), null, 'Relative index of the grid slice in the specified dimension. 0.0 is the first slice, 1.0 is the last slice. Overrides `absolute_index`.'),
    /** Relative isovalue. */
    relative_isovalue: OptionalField(nullable(float), null, 'Relative isovalue.'),
    /** Absolute isovalue. Overrides `relative_isovalue`. */
    absolute_isovalue: OptionalField(nullable(float), null, 'Absolute isovalue. Overrides `relative_isovalue`.'),
};

export const MVSVolumeRepresentationParams = UnionParamsSchema(
    'type',
    'Representation type',
    {
        'isosurface': SimpleParamsSchema(VolumeIsoSurface),
        'grid_slice': SimpleParamsSchema(VolumeGridSlice),
    },
);

const ClipParamsBase = {
    /** Transformation matrix to applied to each point before clipping. For example, can be used to clip volumes in the grid/fractional space. Default is null. */
    check_transform: OptionalField(nullable(Matrix), null, 'Transformation matrix to applied to each point before clipping. For example, can be used to clip volumes in the grid/fractional space. Default is null.'),
    /** Inverts the clipping region. Default is false. */
    invert: OptionalField(bool, false, 'Inverts the clipping region. Default is false'),
    /** Variant of the clip node, either "object" or "pixel". */
    variant: OptionalField(literal('object', 'pixel'), 'pixel', 'Variant of the clip node, either "object" or "pixel"'),
};

export const MVSClipParams = UnionParamsSchema(
    'type',
    'Clip type',
    {
        plane: SimpleParamsSchema({
            ...ClipParamsBase,
            /** Normal vector of the clipping plane. */
            normal: RequiredField(Vector3, 'Normal vector of the clipping plane.'),
            /** Point on the clipping plane. */
            point: RequiredField(Vector3, 'Point on the clipping plane.'),
        }),
        sphere: SimpleParamsSchema({
            ...ClipParamsBase,
            /** Center of the clipping sphere. */
            center: RequiredField(Vector3, 'Center of the clipping sphere.'),
            /** Radius of the clipping sphere. */
            radius: OptionalField(float, 1, 'Radius of the clipping sphere.'),
        }),
        box: SimpleParamsSchema({
            ...ClipParamsBase,
            /** Center of the clipping box. */
            center: RequiredField(Vector3, 'Center of the clipping box.'),
            /** Size of the clipping box. */
            size: OptionalField(Vector3, [1, 1, 1], 'Size of the clipping box.'),
            /** Rotation matrix (3x3 matrix flattened in column major format (j*3+i indexing), this is equivalent to Fortran-order in numpy). This matrix will multiply the structure coordinates from the left. The default value is the identity matrix (corresponds to no rotation). */
            rotation: OptionalField(Matrix, [1, 0, 0, 0, 1, 0, 0, 0, 1], 'Rotation matrix (3x3 matrix flattened in column major format (j*3+i indexing), this is equivalent to Fortran-order in numpy). This matrix will multiply the structure coordinates from the left. The default value is the identity matrix (corresponds to no rotation).'),
        }),
    },
);