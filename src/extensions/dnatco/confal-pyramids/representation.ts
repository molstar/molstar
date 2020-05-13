/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { ConfalPyramids, ConfalPyramidsProvider } from './property';
import { ConfalPyramidsUtil } from './util';
import { ConfalPyramidsTypes as CPT } from './types';
import { Interval } from '../../../mol-data/int';
import { Mesh } from '../../../mol-geo/geometry/mesh/mesh';
import { MeshBuilder } from '../../../mol-geo/geometry/mesh/mesh-builder';
import { PickingId } from '../../../mol-geo/geometry/picking';
import { PrimitiveBuilder } from '../../../mol-geo/primitive/primitive';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { EmptyLoci, Loci } from '../../../mol-model/loci';
import { Structure, StructureProperties, Unit } from '../../../mol-model/structure';
import { CustomProperty } from '../../../mol-model-props/common/custom-property';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../../mol-repr/representation';
import { StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder, UnitsRepresentation } from '../../../mol-repr/structure/representation';
import { StructureGroup, UnitsMeshParams, UnitsMeshVisual, UnitsVisual } from '../../../mol-repr/structure/units-visual';
import { VisualUpdateState } from '../../../mol-repr/util';
import { VisualContext } from '../../../mol-repr/visual';
import { getAltResidueLociFromId } from '../../../mol-repr/structure/visual/util/common';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Theme, ThemeRegistryContext } from '../../../mol-theme/theme';
import { NullLocation } from '../../../mol-model/location';

const t = Mat4.identity();
const w = Vec3.zero();
const mp = Vec3.zero();

function calcMidpoint(mp: Vec3, v: Vec3, w: Vec3) {
    Vec3.sub(mp, v, w);
    Vec3.scale(mp, mp, 0.5);
    Vec3.add(mp, mp, w);
}

function shiftVertex(vec: Vec3, ref: Vec3, scale: number) {
    Vec3.sub(w, vec, ref);
    Vec3.scale(w, w, scale);
    Vec3.add(vec, vec, w);
}

const ConfalPyramidsMeshParams = {
    ...UnitsMeshParams
};
type ConfalPyramidsMeshParams = typeof ConfalPyramidsMeshParams;

function createConfalPyramidsIterator(structureGroup: StructureGroup): LocationIterator {
    const { structure, group } = structureGroup;
    const instanceCount = group.units.length;

    const prop = ConfalPyramidsProvider.get(structure.model).value;
    if (prop === undefined || prop.data === undefined) {
        return LocationIterator(0, 1, () => NullLocation);
    }

    const { locations } = prop.data;

    const getLocation = (groupIndex: number, instanceIndex: number) => {
        if (locations.length <= groupIndex) return NullLocation;
        return locations[groupIndex];
    };
    return LocationIterator(locations.length, instanceCount, getLocation);
}

function createConfalPyramidsMesh(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<ConfalPyramidsMeshParams>, mesh?: Mesh) {
    if (!Unit.isAtomic(unit)) return Mesh.createEmpty(mesh);

    const prop = ConfalPyramidsProvider.get(structure.model).value;
    if (prop === undefined || prop.data === undefined) return Mesh.createEmpty(mesh);

    const { pyramids } = prop.data;
    if (pyramids.length === 0) return Mesh.createEmpty(mesh);

    const mb = MeshBuilder.createState(512, 512, mesh);

    const handler = (pyramid: CPT.Pyramid, first: ConfalPyramidsUtil.FirstResidueAtoms, second: ConfalPyramidsUtil.SecondResidueAtoms, firsLocIndex: number, secondLocIndex: number) => {
        if (firsLocIndex === -1 || secondLocIndex === -1)
            throw new Error('Invalid location index');

        const scale = (pyramid.confal_score - 20.0) / 100.0;
        const O3 = first.O3.pos;
        const OP1 = second.OP1.pos; const OP2 = second.OP2.pos; const O5 = second.O5.pos; const P = second.P.pos;

        shiftVertex(O3, P, scale);
        shiftVertex(OP1, P, scale);
        shiftVertex(OP2, P, scale);
        shiftVertex(O5, P, scale);
        calcMidpoint(mp, O3, O5);

        mb.currentGroup = firsLocIndex;
        let pb = PrimitiveBuilder(3);
        /* Upper part (for first residue in step) */
        pb.add(O3, OP1, OP2);
        pb.add(O3, mp, OP1);
        pb.add(O3, OP2, mp);
        MeshBuilder.addPrimitive(mb, t, pb.getPrimitive());

        /* Lower part (for second residue in step */
        mb.currentGroup = secondLocIndex;
        pb = PrimitiveBuilder(3);
        pb.add(mp, O5, OP1);
        pb.add(mp, OP2, O5);
        pb.add(O5, OP2, OP1);
        MeshBuilder.addPrimitive(mb, t, pb.getPrimitive());
    };

    const walker = new ConfalPyramidsUtil.UnitWalker(structure, unit, handler);
    walker.walk();

    return MeshBuilder.getMesh(mb);
}

function getConfalPyramidLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { groupId, objectId, instanceId } = pickingId;
    if (objectId !== id) return EmptyLoci;

    const { structure } = structureGroup;

    const unit = structureGroup.group.units[instanceId];
    if (!Unit.isAtomic(unit)) return EmptyLoci;

    const prop = ConfalPyramidsProvider.get(structure.model).value;
    if (prop === undefined || prop.data === undefined) return EmptyLoci;

    const { locations } = prop.data;

    if (locations.length <= groupId) return EmptyLoci;
    const altId = StructureProperties.atom.label_alt_id(CPT.toElementLocation(locations[groupId]));
    const rI = unit.residueIndex[locations[groupId].element.element];

    return getAltResidueLociFromId(structure, unit, rI, altId);
}

function eachConfalPyramid(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
    return false; // TODO: Implement me
}

function ConfalPyramidsVisual(materialId: number): UnitsVisual<ConfalPyramidsMeshParams> {
    return UnitsMeshVisual<ConfalPyramidsMeshParams>({
        defaultProps: PD.getDefaultValues(ConfalPyramidsMeshParams),
        createGeometry: createConfalPyramidsMesh,
        createLocationIterator: createConfalPyramidsIterator,
        getLoci: getConfalPyramidLoci,
        eachLocation: eachConfalPyramid,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<ConfalPyramidsMeshParams>, currentProps: PD.Values<ConfalPyramidsMeshParams>) => {
        }
    }, materialId);
}
const ConfalPyramidsVisuals = {
    'confal-pyramids-symbol': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, UnitsMeshParams>) => UnitsRepresentation('Confal Pyramids Symbol Mesh', ctx, getParams, ConfalPyramidsVisual),
};

export const ConfalPyramidsParams = {
    ...UnitsMeshParams
};
export type ConfalPyramidsParams = typeof ConfalPyramidsParams;
export function getConfalPyramidsParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(ConfalPyramidsParams);
}

export type ConfalPyramidsRepresentation = StructureRepresentation<ConfalPyramidsParams>;
export function ConfalPyramidsRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, ConfalPyramidsParams>): ConfalPyramidsRepresentation {
    const repr = Representation.createMulti('Confal Pyramids', ctx, getParams, StructureRepresentationStateBuilder, ConfalPyramidsVisuals as unknown as Representation.Def<Structure, ConfalPyramidsParams>);
    return repr;
}

export const ConfalPyramidsRepresentationProvider = StructureRepresentationProvider({
    name: 'confal-pyramids',
    label: 'Confal Pyramids',
    description: 'Displays schematic depiction of conformer classes and confal values',
    factory: ConfalPyramidsRepresentation,
    getParams: getConfalPyramidsParams,
    defaultValues: PD.getDefaultValues(ConfalPyramidsParams),
    defaultColorTheme: { name: 'confal-pyramids' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (structure: Structure) => structure.models.some(m => ConfalPyramids.isApplicable(m)),
    ensureCustomProperties: {
        attach: (ctx: CustomProperty.Context, structure: Structure) => ConfalPyramidsProvider.attach(ctx, structure.model, void 0, true),
        detach: (data) => ConfalPyramidsProvider.ref(data.model, false),
    }
});
