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
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';

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

export function getTransformFromParams(src:
    | { name: 'matrix', params: { data: Mat4, transpose?: boolean } }
    | { name: 'components', params: { translation: Vec3, axis: Vec3, angle: number } }
) {
    if (src.name === 'matrix') {
        const transform = Mat4();
        Mat4.copy(transform, src.params.data);
        if (src.params.transpose) Mat4.transpose(transform, transform);
        return transform;
    } else {
        const transform = Mat4.fromRotation(Mat4(), src.params.angle * Math.PI / 180, src.params.axis);
        Mat4.setTranslation(transform, src.params.translation);
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
            },
            { isFlat: true }
        ),
    },
    { label: 'Kind' },
);