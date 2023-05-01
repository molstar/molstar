/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { NtCTubeProvider } from './property';
import { NtCTubeSegmentsIterator } from './util';
import { NtCTubeTypes as NTT } from './types';
import { Dnatco } from '../property';
import { DnatcoTypes } from '../types';
import { DnatcoUtil } from '../util';
import { Interval } from '../../../mol-data/int';
import { BaseGeometry, VisualQuality } from '../../../mol-geo/geometry/base';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { addFixedCountDashedCylinder } from '../../../mol-geo/geometry/mesh/builder/cylinder';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { addTube } from '../../../mol-geo/geometry/mesh/builder/tube';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { CylinderProps } from '../../../mol-geo/primitive/cylinder';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { Sphere3D } from '../../../mol-math/geometry/primitives/sphere3d';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { smoothstep } from '../../../mol-math/interpolate';
import { NullLocation } from '../../../mol-model/location';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { Structure, StructureElement, Unit } from '../../../mol-model/structure';
import { structureUnion } from '../../../mol-model/structure/query/utils/structure-set';
import { CustomProperty } from '../../../mol-model-props/common/custom-property';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../../mol-repr/representation';
import { StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder, UnitsRepresentation } from '../../../mol-repr/structure/representation';
import { UnitsMeshParams, UnitsMeshVisual, UnitsVisual } from '../../../mol-repr/structure/units-visual';
import { createCurveSegmentState, CurveSegmentState } from '../../../mol-repr/structure/visual/util/polymer';
import { getStructureQuality, VisualUpdateState } from '../../../mol-repr/util';
import { VisualContext } from '../../../mol-repr/visual';
import { StructureGroup } from '../../../mol-repr/structure/visual/util/common';
import { Theme, ThemeRegistryContext } from '../../../mol-theme/theme';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

const v3add = Vec3.add;
const v3copy = Vec3.copy;
const v3cross = Vec3.cross;
const v3fromArray = Vec3.fromArray;
const v3matchDirection = Vec3.matchDirection;
const v3normalize = Vec3.normalize;
const v3orthogonalize = Vec3.orthogonalize;
const v3scale = Vec3.scale;
const v3slerp = Vec3.slerp;
const v3spline = Vec3.spline;
const v3sub = Vec3.sub;
const v3toArray = Vec3.toArray;

const NtCTubeMeshParams = {
    ...UnitsMeshParams,
    linearSegments: PD.Numeric(4, { min: 2, max: 8, step: 1 }, BaseGeometry.CustomQualityParamInfo),
    radialSegments: PD.Numeric(22, { min: 4, max: 56, step: 2 }, BaseGeometry.CustomQualityParamInfo),
    residueMarkerWidth: PD.Numeric(0.05, { min: 0.01, max: 0.25, step: 0.01 }),
    segmentBoundaryWidth: PD.Numeric(0.05, { min: 0.01, max: 0.25, step: 0.01 }),
};
type NtCTubeMeshParams = typeof NtCTubeMeshParams;

type QualityOptions = Exclude<VisualQuality, 'auto' | 'custom'>;
const LinearSegmentCount: Record<QualityOptions, number> = {
    highest: 6,
    higher: 6,
    high: 4,
    medium: 4,
    low: 3,
    lower: 3,
    lowest: 2,
};
const RadialSegmentCount: Record<QualityOptions, number> = {
    highest: 32,
    higher: 26,
    high: 22,
    medium: 18,
    low: 14,
    lower: 10,
    lowest: 6,
};

const _curvePoint = Vec3();
const _tanA = Vec3();
const _tanB = Vec3();
const _firstTangentVec = Vec3();
const _lastTangentVec = Vec3();
const _firstNormalVec = Vec3();
const _lastNormalVec = Vec3();

const _tmpNormal = Vec3();
const _tangentVec = Vec3();
const _normalVec = Vec3();
const _binormalVec = Vec3();
const _prevNormal = Vec3();
const _nextNormal = Vec3();

function interpolatePointsAndTangents(state: CurveSegmentState, p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, tRange: number[]) {
    const { curvePoints, tangentVectors, linearSegments } = state;
    const tension = 0.5;
    const r = tRange[1] - tRange[0];

    for (let j = 0; j <= linearSegments; ++j) {
        const t = j * r / linearSegments + tRange[0];

        v3spline(_curvePoint, p0, p1, p2, p3, t, tension);
        v3spline(_tanA, p0, p1, p2, p3, t - 0.01, tension);
        v3spline(_tanB, p0, p1, p2, p3, t + 0.01, tension);

        v3toArray(_curvePoint, curvePoints, j * 3);
        v3normalize(_tangentVec, v3sub(_tangentVec, _tanA, _tanB));
        v3toArray(_tangentVec, tangentVectors, j * 3);
    }
}

function interpolateNormals(state: CurveSegmentState, firstDirection: Vec3, lastDirection: Vec3) {
    const { curvePoints, tangentVectors, normalVectors, binormalVectors } = state;

    const n = curvePoints.length / 3;

    v3fromArray(_firstTangentVec, tangentVectors, 0);
    v3fromArray(_lastTangentVec, tangentVectors, (n - 1) * 3);

    v3orthogonalize(_firstNormalVec, _firstTangentVec, firstDirection);
    v3orthogonalize(_lastNormalVec, _lastTangentVec, lastDirection);
    v3matchDirection(_lastNormalVec, _lastNormalVec, _firstNormalVec);

    v3copy(_prevNormal, _firstNormalVec);

    const n1 = n - 1;
    for (let i = 0; i < n; ++i) {
        const j = smoothstep(0, n1, i) * n1;
        const t = i === 0 ? 0 : 1 / (n - j);

        v3fromArray(_tangentVec, tangentVectors, i * 3);

        v3orthogonalize(_normalVec, _tangentVec, v3slerp(_tmpNormal, _prevNormal, _lastNormalVec, t));
        v3toArray(_normalVec, normalVectors, i * 3);

        v3copy(_prevNormal, _normalVec);

        v3normalize(_binormalVec, v3cross(_binormalVec, _tangentVec, _normalVec));
        v3toArray(_binormalVec, binormalVectors, i * 3);
    }

    for (let i = 1; i < n1; ++i) {
        v3fromArray(_prevNormal, normalVectors, (i - 1) * 3);
        v3fromArray(_normalVec, normalVectors, i * 3);
        v3fromArray(_nextNormal, normalVectors, (i + 1) * 3);

        v3scale(_normalVec, v3add(_normalVec, _prevNormal, v3add(_normalVec, _nextNormal, _normalVec)), 1 / 3);
        v3toArray(_normalVec, normalVectors, i * 3);

        v3fromArray(_tangentVec, tangentVectors, i * 3);
        v3normalize(_binormalVec, v3cross(_binormalVec, _tangentVec, _normalVec));
        v3toArray(_binormalVec, binormalVectors, i * 3);
    }
}

function interpolate(state: CurveSegmentState, p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, firstDir: Vec3, lastDir: Vec3, tRange = [0, 1]) {
    interpolatePointsAndTangents(state, p0, p1, p2, p3, tRange);
    interpolateNormals(state, firstDir, lastDir);
}

function createNtCTubeSegmentsIterator(structureGroup: StructureGroup): LocationIterator {
    const { structure, group } = structureGroup;
    const instanceCount = group.units.length;

    const data = NtCTubeProvider.get(structure.model)?.value?.data;
    if (!data) return LocationIterator(0, 1, 1, () => NullLocation);

    const numBlocks = data.data.steps.length * 4;

    const getLocation = (groupId: number, instanceId: number) => {
        if (groupId > numBlocks) return NullLocation;
        const stepIdx = Math.floor(groupId / 4);
        const step = data.data.steps[stepIdx];
        const r = groupId % 4;
        const kind =
            r === 0 ? 'upper' :
                r === 1 ? 'lower' :
                    r === 2 ? 'residue-boundary' : 'segment-boundary';

        return NTT.Location({ step, kind });
    };
    return LocationIterator(totalMeshGroupsCount(data.data.steps) + 1, instanceCount, 1, getLocation);
}

function segmentCount(structure: Structure, props: PD.Values<NtCTubeMeshParams>): { linear: number, radial: number } {
    const quality = props.quality;

    if (quality === 'custom')
        return { linear: props.linearSegments, radial: props.radialSegments };
    else if (quality === 'auto') {
        const autoQuality = getStructureQuality(structure) as QualityOptions;
        return { linear: LinearSegmentCount[autoQuality], radial: RadialSegmentCount[autoQuality] };
    } else
        return { linear: LinearSegmentCount[quality], radial: RadialSegmentCount[quality] };
}

function stepBoundingSphere(step: DnatcoTypes.Step, struLoci: StructureElement.Loci): Sphere3D | undefined {
    const one = DnatcoUtil.residueToLoci(step.auth_asym_id_1, step.auth_seq_id_1, step.label_alt_id_1, step.PDB_ins_code_1, struLoci, 'auth');
    const two = DnatcoUtil.residueToLoci(step.auth_asym_id_2, step.auth_seq_id_2, step.label_alt_id_2, step.PDB_ins_code_2, struLoci, 'auth');

    if (StructureElement.Loci.is(one) && StructureElement.Loci.is(two)) {
        const union = structureUnion(struLoci.structure, [StructureElement.Loci.toStructure(one), StructureElement.Loci.toStructure(two)]);
        return union.boundary.sphere;
    }
    return void 0;
}

function totalMeshGroupsCount(steps: DnatcoTypes.Step[]) {
    // Each segment has two blocks, Residue Boundary marker and a Segment Boundary marker
    return steps.length * 4 - 1; // Subtract one because the last Segment Boundary marker is not drawn
}

function createNtCTubeMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<NtCTubeMeshParams>, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh);

    const prop = NtCTubeProvider.get(structure.model).value;
    if (prop === undefined || prop.data === undefined) return Mesh.createEmpty(mesh);

    const { data } = prop.data;
    if (data.steps.length === 0) return Mesh.createEmpty(mesh);

    const MarkerLinearSegmentCount = 2;
    const segCount = segmentCount(structure, props);
    const vertexCount = Math.floor((segCount.linear * 4 * data.steps.length / structure.model.atomicHierarchy.chains._rowCount) * segCount.radial);
    const chunkSize = Math.floor(vertexCount / 3);
    const diameter = 1.0 * theme.size.props.value;

    const mb = MeshBuilder.createState(vertexCount, chunkSize, mesh);

    const state = createCurveSegmentState(segCount.linear);
    const { curvePoints, normalVectors, binormalVectors, widthValues, heightValues } = state;
    for (let idx = 0; idx <= segCount.linear; idx++) {
        widthValues[idx] = diameter;
        heightValues[idx] = diameter;
    }
    const [normals, binormals] = [binormalVectors, normalVectors]; // Needed so that the tube is not drawn from inside out

    const markerState = createCurveSegmentState(MarkerLinearSegmentCount);
    const { curvePoints: mCurvePoints, normalVectors: mNormalVectors, binormalVectors: mBinormalVectors, widthValues: mWidthValues, heightValues: mHeightValues } = markerState;
    for (let idx = 0; idx <= MarkerLinearSegmentCount; idx++) {
        mWidthValues[idx] = diameter;
        mHeightValues[idx] = diameter;
    }
    const [mNormals, mBinormals] = [mBinormalVectors, mNormalVectors];

    const firstDir = Vec3();
    const lastDir = Vec3();
    const markerDir = Vec3();

    const residueMarkerWidth = props.residueMarkerWidth / 2;
    const it = new NtCTubeSegmentsIterator(structure, unit);
    while (it.hasNext) {
        const segment = it.move();
        if (!segment)
            continue;

        const { p_1, p0, p1, p2, p3, p4, pP } = segment;
        const FirstBlockId = segment.stepIdx * 4;
        const SecondBlockId = FirstBlockId + 1;
        const ResidueMarkerId = FirstBlockId + 2;
        const SegmentBoundaryMarkerId = FirstBlockId + 3;

        const { rmShift, rmPos } = calcResidueMarkerShift(p2, p3, pP);

        if (segment.firstInChain) {
            v3normalize(firstDir, v3sub(firstDir, p2, p1));
            v3normalize(lastDir, v3sub(lastDir, rmPos, p2));
        } else {
            v3copy(firstDir, lastDir);
            v3normalize(lastDir, v3sub(lastDir, rmPos, p2));
        }

        // C5' -> O3' block
        interpolate(state, p0, p1, p2, p3, firstDir, lastDir);
        mb.currentGroup = FirstBlockId;
        addTube(mb, curvePoints, normals, binormals, segCount.linear, segCount.radial, widthValues, heightValues, segment.firstInChain || segment.followsGap, false, 'rounded');

        // O3' -> C5' block
        v3copy(firstDir, lastDir);
        v3normalize(markerDir, v3sub(markerDir, p3, rmPos));
        v3normalize(lastDir, v3sub(lastDir, p4, p3));

        // From O3' to the residue marker
        interpolate(state, p1, p2, p3, p4, firstDir, markerDir, [0, rmShift - residueMarkerWidth]);
        mb.currentGroup = SecondBlockId;
        addTube(mb, curvePoints, normals, binormals, segCount.linear, segCount.radial, widthValues, heightValues, false, false, 'rounded');

        // Residue marker
        interpolate(markerState, p1, p2, p3, p4, markerDir, markerDir, [rmShift - residueMarkerWidth, rmShift + residueMarkerWidth]);
        mb.currentGroup = ResidueMarkerId;
        addTube(mb, mCurvePoints, mNormals, mBinormals, MarkerLinearSegmentCount, segCount.radial, mWidthValues, mHeightValues, false, false, 'rounded');

        if (segment.capEnd) {
            // From the residue marker to C5' of the end
            interpolate(state, p1, p2, p3, p4, markerDir, lastDir, [rmShift + residueMarkerWidth, 1]);
            mb.currentGroup = SecondBlockId;
            addTube(mb, curvePoints, normals, binormals, segCount.linear, segCount.radial, widthValues, heightValues, false, true, 'rounded');
        } else {
            // From the residue marker to C5' of the step boundary marker
            interpolate(state, p1, p2, p3, p4, markerDir, lastDir, [rmShift + residueMarkerWidth, 1 - props.segmentBoundaryWidth]);
            mb.currentGroup = SecondBlockId;
            addTube(mb, curvePoints, normals, binormals, segCount.linear, segCount.radial, widthValues, heightValues, false, false, 'rounded');

            // Step boundary marker
            interpolate(markerState, p1, p2, p3, p4, lastDir, lastDir, [1 - props.segmentBoundaryWidth, 1]);
            mb.currentGroup = SegmentBoundaryMarkerId;
            addTube(mb, mCurvePoints, mNormals, mBinormals, MarkerLinearSegmentCount, segCount.radial, mWidthValues, mHeightValues, false, false, 'rounded');
        }

        if (segment.followsGap) {
            const cylinderProps: CylinderProps = {
                radiusTop: diameter / 2, radiusBottom: diameter / 2, topCap: true, bottomCap: true, radialSegments: segCount.radial,
            };
            mb.currentGroup = FirstBlockId;
            addFixedCountDashedCylinder(mb, p_1, p1, 1, 2 * segCount.linear, false, cylinderProps);
        }
    }

    const boundingSphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, 1.05);

    const m = MeshBuilder.getMesh(mb);
    m.setBoundingSphere(boundingSphere);
    return m;
}

const _rmvCO = Vec3();
const _rmvPO = Vec3();
const _rmPos = Vec3();
const _HalfPi = Math.PI / 2;
function calcResidueMarkerShift(pO: Vec3, pC: Vec3, pP: Vec3): { rmShift: number, rmPos: Vec3 } {
    v3sub(_rmvCO, pC, pO);
    v3sub(_rmvPO, pP, pO);

    // Project position of P atom on the O3' -> C5' vector
    const beta = Vec3.angle(_rmvPO, _rmvCO);
    const alpha = _HalfPi - Math.abs(beta);
    const lengthMO = Math.cos(alpha) * Vec3.magnitude(_rmvPO);
    const shift = lengthMO / Vec3.magnitude(_rmvCO);

    v3scale(_rmvCO, _rmvCO, shift);
    v3add(_rmPos, _rmvCO, pO);

    return { rmShift: shift, rmPos: _rmPos };
}

function getNtCTubeSegmentLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { groupId, objectId, instanceId } = pickingId;
    if (objectId !== id) return EmptyLoci;

    const { structure } = structureGroup;

    const unit = structureGroup.group.units[instanceId];
    if (!Unit.isAtomic(unit)) return EmptyLoci;

    const data = NtCTubeProvider.get(structure.model)?.value?.data ?? undefined;
    if (!data) return EmptyLoci;

    const MeshGroupsCount = totalMeshGroupsCount(data.data.steps);
    if (groupId > MeshGroupsCount) return EmptyLoci;

    const stepIdx = Math.floor(groupId / 4);
    const bs = stepBoundingSphere(data.data.steps[stepIdx], Structure.toStructureElementLoci(structure));

    /*
     * NOTE 1) Each step is drawn with 4 mesh groups. We need to divide/multiply by 4 to convert between steps and mesh groups.
     * NOTE 2) Molstar will create a mesh only for the asymmetric unit. When the entire biological assembly
     *         is displayed, Molstar just copies and transforms the mesh. This means that even though the mesh
     *         might be displayed multiple times, groupIds of the individual blocks in the mesh will be the same.
     *         If there are multiple copies of a mesh, Molstar needs to be able to tell which block belongs to which copy of the mesh.
     *         To do that, Molstar adds an offset to groupIds of the copied meshes. Offset is calculated as follows:
     *
     *         offset = NumberOfBlocks * UnitIndex
     *
     *         "NumberOfBlocks" is the number of valid Location objects got from LocationIterator *or* the greatest groupId set by
     *         the mesh generator - whichever is smaller.
     *
     *         UnitIndex is the index of the Unit the mesh belongs to, starting from 0. (See "unitMap" in the Structure object).
     *         We can also get this index from the value "instanceId" of the "pickingId" object.
     *
     *         If this offset is not applied, picking a piece of one of the copied meshes would actually pick that piece in the original mesh.
     *         This is particularly apparent with highlighting - hovering over items in a copied mesh incorrectly highlights those items in the source mesh.
     *
     *         Molstar can take advantage of the fact that ElementLoci has a reference to the Unit object attached to it. Since we cannot attach ElementLoci
     *         to a step, we need to calculate the offseted groupId here and pass it as part of the DataLoci.
     */
    const offsetGroupId = stepIdx * 4 + (MeshGroupsCount + 1) * instanceId;
    return NTT.Loci(data.data.steps, [stepIdx], [offsetGroupId], bs);
}

function eachNtCTubeSegment(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
    if (NTT.isLoci(loci)) {
        const offsetGroupId = loci.elements[0];
        return apply(Interval.ofBounds(offsetGroupId, offsetGroupId + 4));
    }
    return false;
}

function NtCTubeVisual(materialId: number): UnitsVisual<NtCTubeMeshParams> {
    return UnitsMeshVisual<NtCTubeMeshParams>({
        defaultProps: PD.getDefaultValues(NtCTubeMeshParams),
        createGeometry: createNtCTubeMesh,
        createLocationIterator: createNtCTubeSegmentsIterator,
        getLoci: getNtCTubeSegmentLoci,
        eachLocation: eachNtCTubeSegment,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<NtCTubeMeshParams>, currentProps: PD.Values<NtCTubeMeshParams>) => {
            state.createGeometry = (
                newProps.quality !== currentProps.quality ||
                newProps.residueMarkerWidth !== currentProps.residueMarkerWidth ||
                newProps.segmentBoundaryWidth !== currentProps.segmentBoundaryWidth ||
                newProps.doubleSided !== currentProps.doubleSided ||
                newProps.alpha !== currentProps.alpha ||
                newProps.linearSegments !== currentProps.linearSegments ||
                newProps.radialSegments !== currentProps.radialSegments
            );
        }
    }, materialId);

}
const NtCTubeVisuals = {
    'ntc-tube-symbol': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, NtCTubeMeshParams>) => UnitsRepresentation('NtC Tube Mesh', ctx, getParams, NtCTubeVisual),
};

export const NtCTubeParams = {
    ...NtCTubeMeshParams
};
export type NtCTubeParams = typeof NtCTubeParams;
export function getNtCTubeParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(NtCTubeParams);
}

export type NtCTubeRepresentation = StructureRepresentation<NtCTubeParams>;
export function NtCTubeRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, NtCTubeParams>): NtCTubeRepresentation {
    return Representation.createMulti('NtC Tube', ctx, getParams, StructureRepresentationStateBuilder, NtCTubeVisuals as unknown as Representation.Def<Structure, NtCTubeParams>);
}

export const NtCTubeRepresentationProvider = StructureRepresentationProvider({
    name: 'ntc-tube',
    label: 'NtC Tube',
    description: 'Displays schematic representation of NtC conformers',
    factory: NtCTubeRepresentation,
    getParams: getNtCTubeParams,
    defaultValues: PD.getDefaultValues(NtCTubeParams),
    defaultColorTheme: { name: 'ntc-tube' },
    defaultSizeTheme: { name: 'uniform', props: { value: 2.0 } },
    isApplicable: (structure: Structure) => structure.models.every(m => Dnatco.isApplicable(m)),
    ensureCustomProperties: {
        attach: async (ctx: CustomProperty.Context, structure: Structure) => structure.models.forEach(m => NtCTubeProvider.attach(ctx, m, void 0, true)),
        detach: (data) => data.models.forEach(m => NtCTubeProvider.ref(m, false)),
    },
});
