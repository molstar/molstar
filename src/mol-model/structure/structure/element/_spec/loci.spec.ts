/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column, Table } from '../../../../../mol-data/db';
import { RuntimeContext } from '../../../../../mol-task';
import { createModels } from '../../../../../mol-model-formats/structure/basic/parser';
import { BasicSchema, createBasic } from '../../../../../mol-model-formats/structure/basic/schema';
import { EntityBuilder } from '../../../../../mol-model-formats/structure/common/entity';
import { ComponentBuilder } from '../../../../../mol-model-formats/structure/common/component';
import { MoleculeType } from '../../../model/types';
import { Box3D, Sphere3D } from '../../../../../mol-math/geometry';
import { SymmetryOperator } from '../../../../../mol-math/geometry/symmetry-operator';
import { Mat4, Vec3 } from '../../../../../mol-math/linear-algebra';
import { OrderedSet } from '../../../../../mol-data/int';
import { Structure } from '../../structure';
import { Unit } from '../../unit';
import { UnitIndex } from '../element';
import { Loci } from '../loci';

/** Helper: build a single-chain, 4-residue (one CA atom each) `BasicData` for `createModels`. */
function buildSingleChainBasic() {
    const coords: [number, number, number][] = [
        [0, 0, 0],
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 2],
    ];
    const count = coords.length;

    const auth_asym_id = Column.ofConst('A', count, Column.Schema.str);
    const auth_atom_id = Column.ofConst('CA', count, Column.Schema.str);
    const auth_comp_id = Column.ofConst('ALA', count, Column.Schema.str);
    const auth_seq_id = Column.ofIntArray(coords.map((_, i) => i + 1));
    const type_symbol = Column.ofConst('C', count, Column.Schema.str);

    const atom_site = Table.ofPartialColumns(BasicSchema.atom_site, {
        auth_asym_id,
        auth_atom_id,
        auth_comp_id,
        auth_seq_id,
        Cartn_x: Column.ofFloatArray(Float32Array.from(coords.map(c => c[0]))),
        Cartn_y: Column.ofFloatArray(Float32Array.from(coords.map(c => c[1]))),
        Cartn_z: Column.ofFloatArray(Float32Array.from(coords.map(c => c[2]))),
        id: Column.ofIntArray(coords.map((_, i) => i)),

        label_asym_id: auth_asym_id,
        label_atom_id: auth_atom_id,
        label_comp_id: auth_comp_id,
        label_seq_id: auth_seq_id,
        label_entity_id: Column.ofConst('1', count, Column.Schema.str),

        occupancy: Column.ofConst(1, count, Column.Schema.float),
        type_symbol,

        pdbx_PDB_model_num: Column.ofConst(1, count, Column.Schema.int),
    }, count);

    const entityBuilder = new EntityBuilder();
    entityBuilder.getEntityId('ALA', MoleculeType.Protein, 'A');

    const componentBuilder = new ComponentBuilder(auth_seq_id, auth_atom_id);
    componentBuilder.add('ALA', 0);

    return createBasic({
        entity: entityBuilder.getEntityTable(),
        chem_comp: componentBuilder.getChemCompTable(),
        atom_site
    });
}

async function buildSingleChainStructure(): Promise<Structure> {
    const basic = buildSingleChainBasic();
    const trajectory = await createModels(basic, { kind: 'test', name: 'synthetic-test', data: undefined }, RuntimeContext.Synchronous);
    return Structure.ofModel(trajectory.representative);
}

function collectPositions(elements: ReadonlyArray<{ unit: Unit, indices: OrderedSet<UnitIndex> }>, transform?: Mat4): Vec3[] {
    const positions: Vec3[] = [];
    for (const { unit, indices } of elements) {
        const { elements: unitElements, conformation } = unit;
        for (let i = 0, _i = OrderedSet.size(indices); i < _i; i++) {
            const eI = unitElements[OrderedSet.getAt(indices, i)];
            const p = Vec3();
            conformation.position(eI, p);
            if (transform) Vec3.transformMat4(p, p, transform);
            positions.push(p);
        }
    }
    return positions;
}

function expectContains(boundary: { box: Box3D, sphere: Sphere3D }, positions: Vec3[], eps = 1e-4) {
    for (const p of positions) {
        expect(Vec3.distance(boundary.sphere.center, p)).toBeLessThanOrEqual(boundary.sphere.radius + eps);
        expect(p[0]).toBeGreaterThanOrEqual(boundary.box.min[0] - eps);
        expect(p[1]).toBeGreaterThanOrEqual(boundary.box.min[1] - eps);
        expect(p[2]).toBeGreaterThanOrEqual(boundary.box.min[2] - eps);
        expect(p[0]).toBeLessThanOrEqual(boundary.box.max[0] + eps);
        expect(p[1]).toBeLessThanOrEqual(boundary.box.max[1] + eps);
        expect(p[2]).toBeLessThanOrEqual(boundary.box.max[2] + eps);
    }
}

function expectVec3Close(a: Vec3, b: Vec3, precision = 5) {
    expect(a[0]).toBeCloseTo(b[0], precision);
    expect(a[1]).toBeCloseTo(b[1], precision);
    expect(a[2]).toBeCloseTo(b[2], precision);
}

/**
 * The whole-unit/whole-structure fast paths reduce a cached sphere via its stored direction
 * extrema (`BoundaryHelper.includeSphere`/`radiusSphere`), same as the pre-existing
 * `computeStructureBoundary`. This reproduces the sphere center exactly, but can yield a
 * slightly smaller radius than a naive point transform (the true farthest point of a bounded
 * atom sphere isn't necessarily one of the fixed sample directions) - so radius comparisons use
 * a relative tolerance instead of exact/close-to-N-decimals equality.
 */
function expectRadiusApprox(actual: number, expected: number, relTol = 0.05) {
    expect(Math.abs(actual - expected) / expected).toBeLessThanOrEqual(relTol);
}

describe('StructureElement.Loci.getBoundary', () => {
    let structure0: Structure;
    let unit0: Unit;
    let unit1: Unit;
    let multiStructure: Structure;
    let operator: SymmetryOperator;

    let wholeIndices: OrderedSet<UnitIndex>;
    let partialIndices: OrderedSet<UnitIndex>;

    beforeAll(async () => {
        structure0 = await buildSingleChainStructure();

        // Guard: fixture is expected to parse into a single 4-atom unit.
        expect(structure0.units.length).toBe(1);
        unit0 = structure0.units[0];
        expect(unit0.elements.length).toBe(4);

        operator = SymmetryOperator.create('test-op', Mat4.fromTranslation(Mat4(), Vec3.create(10, 0, 0)));
        unit1 = unit0.applyOperator(unit0.id + 1, operator, true);

        multiStructure = Structure.create([unit0, unit1]);

        wholeIndices = OrderedSet.ofBounds<UnitIndex>(0 as UnitIndex, 4 as UnitIndex);
        partialIndices = OrderedSet.ofBounds<UnitIndex>(0 as UnitIndex, 2 as UnitIndex);
    });

    it('whole-structure loci reuses Structure.boundary exactly', () => {
        const loci = Loci.all(structure0);
        const b = Loci.getBoundary(loci);

        expect(b.box.min).toEqual(structure0.boundary.box.min);
        expect(b.box.max).toEqual(structure0.boundary.box.max);
        expect(b.sphere.center).toEqual(structure0.boundary.sphere.center);
        expect(b.sphere.radius).toBeCloseTo(structure0.boundary.sphere.radius, 6);
    });

    it('whole-unit loci (identity operator, not whole structure) matches unit.boundary and contains all atoms', () => {
        const loci = Loci(multiStructure, [{ unit: unit0, indices: wholeIndices }]);
        expect(Loci.isWholeStructure(loci)).toBe(false);

        const b = Loci.getBoundary(loci);
        expectVec3Close(b.sphere.center, unit0.boundary.sphere.center);
        expectRadiusApprox(b.sphere.radius, unit0.boundary.sphere.radius);
        expectContains(b, collectPositions([{ unit: unit0, indices: wholeIndices }]));
    });

    it('partial-unit loci falls back to the per-atom computation and only contains selected atoms', () => {
        const loci = Loci(multiStructure, [{ unit: unit0, indices: partialIndices }]);
        const b = Loci.getBoundary(loci);

        expectContains(b, collectPositions([{ unit: unit0, indices: partialIndices }]));

        const wholeUnitLoci = Loci(multiStructure, [{ unit: unit0, indices: wholeIndices }]);
        const wholeUnitBoundary = Loci.getBoundary(wholeUnitLoci);
        // partial selection's sphere should not be larger than the whole unit's
        expect(b.sphere.radius).toBeLessThanOrEqual(wholeUnitBoundary.sphere.radius + 1e-6);
    });

    it('mixed loci (whole unit + partial unit) combines both fast and slow paths', () => {
        const elements = [
            { unit: unit0, indices: wholeIndices },
            { unit: unit1, indices: partialIndices },
        ];
        const loci = Loci(multiStructure, elements);
        const b = Loci.getBoundary(loci);

        expectContains(b, collectPositions(elements));
    });

    it('whole-unit loci with a non-identity operator applies the operator matrix to the cached unit boundary', () => {
        const loci = Loci(multiStructure, [{ unit: unit1, indices: wholeIndices }]);
        const b = Loci.getBoundary(loci);

        const expectedSphere = Sphere3D.transform(Sphere3D(), unit0.boundary.sphere, operator.matrix);
        expectVec3Close(b.sphere.center, expectedSphere.center);
        expectRadiusApprox(b.sphere.radius, expectedSphere.radius);

        expectContains(b, collectPositions([{ unit: unit1, indices: wholeIndices }]));
    });

    it('boundary param is filled in place', () => {
        const out = { box: Box3D(), sphere: Sphere3D() };
        const loci = Loci.all(structure0);
        const ret = Loci.getBoundary(loci, out);

        expect(ret.box).toBe(out.box);
        expect(ret.sphere).toBe(out.sphere);
        expect(out.box.min).toEqual(structure0.boundary.box.min);
        expect(out.sphere.center).toEqual(structure0.boundary.sphere.center);
    });

    it('does not alias the cached Structure.boundary object', () => {
        const loci = Loci.all(structure0);
        const originalMin = Vec3.clone(structure0.boundary.box.min);

        const fresh = Loci.getBoundary(loci);
        fresh.box.min[0] = 9999;

        expect(structure0.boundary.box.min).toEqual(originalMin);
    });
});
