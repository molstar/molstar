/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject } from '../objects';
import { DistanceData } from '../../mol-repr/shape/loci/distance';
import { LabelData } from '../../mol-repr/shape/loci/label';
import { OrientationData } from '../../mol-repr/shape/loci/orientation';
import { AngleData } from '../../mol-repr/shape/loci/angle';
import { DihedralData } from '../../mol-repr/shape/loci/dihedral';
import { PlaneData } from '../../mol-repr/shape/loci/plane';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { EPSILON, Mat4, Vec3 } from '../../mol-math/linear-algebra';

export function getDistanceDataFromStructureSelections(s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>): DistanceData {
    const lociA = s[0].loci;
    const lociB = s[1].loci;
    return { pairs: [{ loci: [lociA, lociB] as const }] };
}

export function getAngleDataFromStructureSelections(s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>): AngleData {
    const lociA = s[0].loci;
    const lociB = s[1].loci;
    const lociC = s[2].loci;
    return { triples: [{ loci: [lociA, lociB, lociC] as const }] };
}

export function getDihedralDataFromStructureSelections(s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>): DihedralData {
    const lociA = s[0].loci;
    const lociB = s[1].loci;
    const lociC = s[2].loci;
    const lociD = s[3].loci;
    return { quads: [{ loci: [lociA, lociB, lociC, lociD] as const }] };
}

export function getLabelDataFromStructureSelections(s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>): LabelData {
    const loci = s[0].loci;
    return { infos: [{ loci }] };
}

export function getOrientationDataFromStructureSelections(s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>): OrientationData {
    return { locis: s.map(v => v.loci) };
}

export function getPlaneDataFromStructureSelections(s: ReadonlyArray<PluginStateObject.Molecule.Structure.SelectionEntry>): PlaneData {
    return { locis: s.map(v => v.loci) };
}

const GetTransformState = {
    center: Vec3(),
    rotation: Mat4(),
    translation1: Mat4(),
    translation2: Mat4(),
    local: Mat4(),
};

type TransformParams =
    | { name: 'matrix', params: { data: Mat4, transpose?: boolean } }
    | { name: 'components', params: { translation: Vec3, axis: Vec3, angle: number, localAxis: Vec3, localAngle: number } }

export function transformParamsNeedCenter(src: TransformParams) {
    if (src.name === 'components' && src.params.localAngle > EPSILON) {
        return true;
    }
    return false;
}

export function getTransformFromParams(src: TransformParams, center: Vec3) {
    if (src.name === 'matrix') {
        const transform = Mat4();
        Mat4.copy(transform, src.params.data);
        if (src.params.transpose) Mat4.transpose(transform, transform);
        return transform;
    } else {
        if (src.params.localAngle) {
            Mat4.fromTranslation(GetTransformState.translation1, center);
            Mat4.fromRotation(GetTransformState.rotation, src.params.localAngle * Math.PI / 180, src.params.localAxis);
            Mat4.fromTranslation(GetTransformState.translation2, Vec3.negate(GetTransformState.center, center));
            Mat4.mul3(
                GetTransformState.local,
                GetTransformState.translation1,
                GetTransformState.rotation,
                GetTransformState.translation2
            );
        } else {
            Mat4.setIdentity(GetTransformState.local);
        }

        const transform = Mat4.fromRotation(Mat4(), src.params.angle * Math.PI / 180, src.params.axis);
        Mat4.setTranslation(transform, src.params.translation);
        Mat4.mul(transform, transform, GetTransformState.local);
        return transform;
    }
}

export const TransformParam = PD.MappedStatic(
    'matrix',
    {
        matrix: PD.Group(
            {
                data: PD.Mat4(Mat4.identity()),
                transpose: PD.Boolean(false),
            },
            { isFlat: true }
        ),
        components: PD.Group(
            {
                translation: PD.Vec3(Vec3.create(0, 0, 0)),
                axis: PD.Vec3(Vec3.create(1, 0, 0)),
                angle: PD.Numeric(0, { min: -360, max: 360, step: 1 }, { description: 'Angle in Degrees' }),
                localAxis: PD.Vec3(Vec3.create(1, 0, 0)),
                localAngle: PD.Numeric(0, { min: -360, max: 360, step: 1 }, { description: 'Angle in Degrees' }),
            },
            { isFlat: true }
        ),
    },
    { label: 'Kind' },
);