/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bool, float, nullable, OptionalField } from '../generic/field-schema';
import { SimpleParamsSchema, UnionParamsSchema } from '../generic/params-schema';

const Cartoon = {
    /** Scales the corresponding visuals */
    size_factor: OptionalField(float, 1, 'Scales the corresponding visuals.'),
    /** Simplify corkscrew helices to tubes. */
    tubular_helices: OptionalField(bool, false, 'Simplify corkscrew helices to tubes.'),
};

const BallAndStick = {
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
        ball_and_stick: SimpleParamsSchema(BallAndStick),
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

export const MVSVolumeRepresentationParams = UnionParamsSchema(
    'type',
    'Representation type',
    {
        'isosurface': SimpleParamsSchema(VolumeIsoSurface),
    },
);