/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Structure from './structure'
import { StructureSelection, QueryContext } from '../query'
import { ModelSymmetry } from '../model'
import { Task, RuntimeContext } from 'mol-task';
import { SortedArray } from 'mol-data/int';
import Unit from './unit';
import { EquivalenceClasses } from 'mol-data/util';
import { Vec3 } from 'mol-math/linear-algebra';
import { SymmetryOperator, Spacegroup, SpacegroupCell } from 'mol-math/geometry';

namespace StructureSymmetry {
    export function buildAssembly(structure: Structure, asmName: string) {
        return Task.create('Build Assembly', async ctx => {
            const models = structure.models;
            if (models.length !== 1) throw new Error('Can only build assemblies from structures based on 1 model.');

            const assembly = ModelSymmetry.findAssembly(models[0], asmName);
            if (!assembly) throw new Error(`Assembly '${asmName}' is not defined.`);

            const assembler = Structure.Builder();

            const queryCtx = new QueryContext(structure);

            for (const g of assembly.operatorGroups) {
                const selection = g.selector(queryCtx);
                if (StructureSelection.structureCount(selection) === 0) {
                    continue;
                }
                const { units } = StructureSelection.unionStructure(selection);

                for (const oper of g.operators) {
                    for (const unit of units) {
                        assembler.addWithOperator(unit, oper);
                    }
                }
            }

            return assembler.getStructure();
        });
    }

    export function builderSymmetryMates(structure: Structure, radius: number) {
        return Task.create('Find Symmetry Mates', ctx => findMatesRadius(ctx, structure, radius));
    }

    export function buildSymmetryRange(structure: Structure, ijkMin: Vec3, ijkMax: Vec3) {
        return Task.create('Build Symmetry', ctx => findSymmetryRange(ctx, structure, ijkMin, ijkMax));
    }

    /** Builds NCS structure, returns the original if NCS operators are not present. */
    export function buildNcs(structure: Structure) {
        return Task.create('Build NCS', ctx => _buildNCS(ctx, structure));
    }

    export function areUnitsEquivalent(a: Unit, b: Unit) {
        return a.invariantId === b.invariantId && a.model.id === b.model.id && SortedArray.areEqual(a.elements, b.elements);
    }

    export function UnitEquivalenceBuilder() {
        return EquivalenceClasses<number, Unit>(Unit.hashUnit, areUnitsEquivalent);
    }

    export function computeTransformGroups(s: Structure): ReadonlyArray<Unit.SymmetryGroup> {
        const groups = UnitEquivalenceBuilder();
        for (const u of s.units) groups.add(u.id, u);

        const ret: Unit.SymmetryGroup[] = [];
        for (const eqUnits of groups.groups) {
            ret.push(Unit.SymmetryGroup(eqUnits.map(id => s.unitMap.get(id))))
        }

        return ret;
    }

    /** Checks if transform groups are equal up to their unit's transformations */
    export function areTransformGroupsEquivalent(a: ReadonlyArray<Unit.SymmetryGroup>, b: ReadonlyArray<Unit.SymmetryGroup>) {
        if (a.length !== b.length) return false
        for (let i = 0, il = a.length; i < il; ++i) {
            if (a[i].hashCode !== b[i].hashCode) return false
        }
        return true
    }
}

function getOperators(symmetry: ModelSymmetry, ijkMin: Vec3, ijkMax: Vec3) {
    const operators: SymmetryOperator[] = symmetry._operators_333 || [];
    const { spacegroup } = symmetry;
    if (operators.length === 0) {
        operators[0] = Spacegroup.getSymmetryOperator(spacegroup, 0, 0, 0, 0)
        for (let op = 0; op < spacegroup.operators.length; op++) {
            for (let i = ijkMin[0]; i < ijkMax[0]; i++) {
                for (let j = ijkMin[1]; j < ijkMax[1]; j++) {
                    for (let k = ijkMin[2]; k < ijkMax[2]; k++) {
                        // we have added identity as the 1st operator.
                        if (op === 0 && i === 0 && j === 0 && k === 0) continue;
                        operators[operators.length] = Spacegroup.getSymmetryOperator(spacegroup, op, i, j, k);
                    }
                }
            }
        }
        symmetry._operators_333 = operators;
    }
    return operators;
}

function assembleOperators(structure: Structure, operators: ReadonlyArray<SymmetryOperator>) {
    const assembler = Structure.Builder();
    const { units } = structure;
    for (const oper of operators) {
        for (const unit of units) {
            assembler.addWithOperator(unit, oper);
        }
    }
    return assembler.getStructure();
}

async function _buildNCS(ctx: RuntimeContext, structure: Structure) {
    const models = structure.models;
    if (models.length !== 1) throw new Error('Can only build NCS from structures based on 1 model.');

    const operators = models[0].symmetry.ncsOperators;
    if (!operators || !operators.length) return structure;
    return assembleOperators(structure, operators);
}

async function findSymmetryRange(ctx: RuntimeContext, structure: Structure, ijkMin: Vec3, ijkMax: Vec3) {
    const models = structure.models;
    if (models.length !== 1) throw new Error('Can only build symmetries from structures based on 1 model.');

    const { spacegroup } = models[0].symmetry;
    if (SpacegroupCell.isZero(spacegroup.cell)) return structure;

    const operators = getOperators(models[0].symmetry, ijkMin, ijkMax);
    return assembleOperators(structure, operators);
}

async function findMatesRadius(ctx: RuntimeContext, structure: Structure, radius: number) {
    const models = structure.models;
    if (models.length !== 1) throw new Error('Can only build symmetries from structures based on 1 model.');

    const symmetry = models[0].symmetry;
    const { spacegroup } = symmetry;
    if (SpacegroupCell.isZero(spacegroup.cell)) return structure;

    if (ctx.shouldUpdate) await ctx.update('Initialing...');
    const operators = getOperators(symmetry, Vec3.create(-3, -3, -3), Vec3.create(3, 3, 3));
    const lookup = structure.lookup3d;

    const assembler = Structure.Builder();

    const { units } = structure;
    const center = Vec3.zero();
    for (const oper of operators) {
        for (const unit of units) {
            const boundingSphere = unit.lookup3d.boundary.sphere;
            Vec3.transformMat4(center, boundingSphere.center, oper.matrix);

            const closeUnits = lookup.findUnitIndices(center[0], center[1], center[2], boundingSphere.radius + radius);
            for (let uI = 0, _uI = closeUnits.count; uI < _uI; uI++) {
                const closeUnit = units[closeUnits.indices[uI]];
                if (!closeUnit.lookup3d.check(center[0], center[1], center[2], boundingSphere.radius + radius)) continue;
                assembler.addWithOperator(unit, oper);
            }
        }
        if (ctx.shouldUpdate) await ctx.update('Building symmetry...');
    }


    return assembler.getStructure();
}

export default StructureSymmetry;