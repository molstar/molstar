/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Multi-visual representation for parsed EMDB-SFF data: a single mol*
 * Representation that contains both a Mesh sub-visual and a Wireframe
 * (Lines) sub-visual. Modeled after the gaussian-surface representation
 * (`mol-repr/structure/representation/gaussian-surface.ts`), so that the
 * UI shows one representation node with a `Visuals` MultiSelect to toggle
 * mesh / wireframe / both.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { Lines } from '../../mol-geo/geometry/lines/lines';
import { BaseGeometry } from '../../mol-geo/geometry/base';
import { Shape } from '../../mol-model/shape';
import { Color } from '../../mol-util/color';
import { Material } from '../../mol-util/material';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../mol-repr/representation';
import { ShapeRepresentation } from '../../mol-repr/shape/representation';
import { SffData } from '../../mol-io/reader/hff/schema';
import { BuiltMesh, buildLinesFromMesh, buildMesh, colourToColor, segmentLabel } from './model';


export const SffRepresentationParams = {
    ...Mesh.Params,
    ...Lines.Params,
    // SFF mesh segmentations are typically thin oriented surfaces (membranes,
    // organelles); render both sides by default and keep the back face in the
    // segment colour.
    doubleSided: PD.Boolean(true, BaseGeometry.CustomQualityParamInfo),
    interior: PD.Group({
        color: PD.Color(Color.fromRgb(76, 76, 76)),
        colorStrength: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
        substance: Material.getParam(),
        substanceStrength: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
    }),
    visuals: PD.MultiSelect<'mesh' | 'wire'>(
        ['mesh'],
        [['mesh', 'Mesh'], ['wire', 'Wireframe']],
    ),
};
export type SffRepresentationParams = typeof SffRepresentationParams;


/** Per-(SffData) cache shared between mesh and wire visuals; built lazily. */
function makeBuiltCache() {
    let key: SffData | undefined;
    let built: BuiltMesh | undefined;
    return (data: SffData): BuiltMesh => {
        if (data !== key) {
            built = buildMesh(data);
            key = data;
        }
        return built!;
    };
}


function makeMeshShapeGetter() {
    const cache = makeBuiltCache();
    let _shape: Shape<Mesh> | undefined;
    let _key: SffData | undefined;
    return async (_runtime: any, data: SffData) => {
        if (_shape && _key === data) return _shape;
        const built = cache(data);
        const colors = built.segmentByGroup.map(s => colourToColor(s.colour));
        const labels = built.segmentByGroup.map(s => segmentLabel(s));
        const label = data.name?.trim() || 'EMDB-SFF';
        _shape = Shape.create(
            label,
            data,
            built.mesh,
            (g: number) => colors[g] ?? Color.fromNormalizedRgb(0.7, 0.7, 0.7),
            () => 1,
            (g: number) => labels[g] ?? `Segment ${g}`,
        );
        _key = data;
        return _shape;
    };
}


function makeWireShapeGetter() {
    const cache = makeBuiltCache();
    let _shape: Shape<Lines> | undefined;
    let _key: SffData | undefined;
    return async (_runtime: any, data: SffData) => {
        if (_shape && _key === data) return _shape;
        const built = cache(data);
        const lines = buildLinesFromMesh(built);
        const colors = built.segmentByGroup.map(s => colourToColor(s.colour));
        const labels = built.segmentByGroup.map(s => segmentLabel(s));
        const label = data.name?.trim() || 'EMDB-SFF';
        _shape = Shape.create(
            `${label} (wire)`,
            data,
            lines,
            (g: number) => colors[g] ?? Color.fromNormalizedRgb(0.7, 0.7, 0.7),
            () => 1,
            (g: number) => labels[g] ?? `Segment ${g}`,
        );
        _key = data;
        return _shape;
    };
}


// `ShapeRepresentation`'s generics force a single Geometry per representation,
// but `Representation.createMulti` requires sub-reprs to share the same `P`.
// In practice each sub-rep ignores foreign params (mesh ignores Lines params
// and vice versa), so the cast is safe — same trick as gaussian-surface uses.
const SffVisuals: Representation.Def<SffData, SffRepresentationParams> = {
    'mesh': (_ctx: RepresentationContext) =>
        ShapeRepresentation(makeMeshShapeGetter(), Mesh.Utils as any) as unknown as Representation<SffData, SffRepresentationParams>,
    'wire': (_ctx: RepresentationContext) =>
        ShapeRepresentation(makeWireShapeGetter(), Lines.Utils as any) as unknown as Representation<SffData, SffRepresentationParams>,
};


export type SffRepresentation = Representation<SffData, SffRepresentationParams>;
export function SffRepresentation(
    ctx: RepresentationContext,
    getParams: RepresentationParamsGetter<SffData, SffRepresentationParams>,
): SffRepresentation {
    return Representation.createMulti(
        'SFF',
        ctx,
        getParams,
        Representation.StateBuilder,
        SffVisuals,
    );
}

export function getSffRepresentationParams(_: unknown, _data: SffData) {
    return SffRepresentationParams;
}
