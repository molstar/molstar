/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4, Tensor } from 'mol-math/linear-algebra'
import SymmetryOperator from 'mol-math/geometry/symmetry-operator'
import Format from '../../format'
import { Assembly, OperatorGroup, OperatorGroups } from '../../properties/symmetry'
import { Queries as Q, Query } from '../../../query'

import mmCIF_Format = Format.mmCIF

export default function create(format: mmCIF_Format): ReadonlyArray<Assembly> {
    const { pdbx_struct_assembly } = format.data;
    if (!pdbx_struct_assembly._rowCount) return [];

    const matrices = getMatrices(format);
    const assemblies: Assembly[] = [];
    for (let i = 0; i < pdbx_struct_assembly._rowCount; i++) {
        assemblies[assemblies.length] = createAssembly(format, i, matrices);
    }
    return assemblies;
}

type Matrices = Map<string, Mat4>
type Generator = { expression: string, asymIds: string[] }

function createAssembly(format: mmCIF_Format, index: number, matrices: Matrices): Assembly {
    const { pdbx_struct_assembly, pdbx_struct_assembly_gen } = format.data;

    const id = pdbx_struct_assembly.id.value(index);
    const details = pdbx_struct_assembly.details.value(index);
    const generators: Generator[] = [];

    const { assembly_id, oper_expression, asym_id_list } = pdbx_struct_assembly_gen;

    for (let i = 0, _i = pdbx_struct_assembly_gen._rowCount; i < _i; i++) {
        if (assembly_id.value(i) !== id) continue;
        generators[generators.length] = {
            expression: oper_expression.value(i),
            asymIds: asym_id_list.value(i)
        };
    }

    return Assembly.create(id, details, operatorGroupsProvider(generators, matrices));
}

function operatorGroupsProvider(generators: Generator[], matrices: Matrices): () => OperatorGroups {
    return () => {
        const groups: OperatorGroup[] = [];

        let operatorOffset = 0;
        for (let i = 0; i < generators.length; i++) {
            const gen = generators[i];
            const operatorList = parseOperatorList(gen.expression);
            const operatorNames = expandOperators(operatorList);
            const operators = getAssemblyOperators(matrices, operatorNames, operatorOffset);
            const selector = Query(Q.generators.atoms({ chainTest: Q.pred.and(
                Q.pred.eq(Q.props.unit.operator_name, SymmetryOperator.DefaultName),
                Q.pred.inSet(Q.props.chain.label_asym_id, gen.asymIds)
            )}));
            groups[groups.length] = { selector, operators };
            operatorOffset += operators.length;
        }

        return groups;
    }
}

function getMatrices({ data }: mmCIF_Format): Matrices {
    const { pdbx_struct_oper_list } = data;
    const { id, matrix, vector, _schema } = pdbx_struct_oper_list;
    const matrices = new Map<string, Mat4>();

    for (let i = 0, _i = pdbx_struct_oper_list._rowCount; i < _i; i++) {
        const m = Tensor.toMat4(_schema.matrix.space, matrix.value(i));
        const t = Tensor.toVec3(_schema.vector.space, vector.value(i));
        Mat4.setTranslation(m, t);
        Mat4.setValue(m, 3, 3, 1);
        matrices.set(id.value(i), m);
    }

    return matrices;
}

function expandOperators(operatorList: string[][]) {
    const ops: string[][] = [];
    const currentOp: string[] = [];
    for (let i = 0; i < operatorList.length; i++) currentOp[i] = '';
    expandOperators1(operatorList, ops, operatorList.length - 1, currentOp);
    return ops;
}

function expandOperators1(operatorNames: string[][], list: string[][], i: number, current: string[]) {
    if (i < 0) {
        list[list.length] = current.slice(0);
        return;
    }

    let ops = operatorNames[i], len = ops.length;
    for (let j = 0; j < len; j++) {
        current[i] = ops[j];
        expandOperators1(operatorNames, list, i - 1, current);
    }
}

function getAssemblyOperators(matrices: Matrices, operatorNames: string[][], startIndex: number) {
    const operators: SymmetryOperator[] = [];

    let index = startIndex;
    for (let op of operatorNames) {
        let m = Mat4.identity();
        for (let i = 0; i < op.length; i++) {
            Mat4.mul(m, m, matrices.get(op[i])!);
        }
        index++;
        operators[operators.length] = SymmetryOperator.create(`A-${index}`, m);
    }

    return operators;
}

function parseOperatorList(value: string): string[][] {
    // '(X0)(1-5)' becomes [['X0'], ['1', '2', '3', '4', '5']]
    // kudos to Glen van Ginkel.

    const oeRegex = /\(?([^\(\)]+)\)?]*/g, groups: string[] = [], ret: string[][] = [];

    let g: any;
    while (g = oeRegex.exec(value)) groups[groups.length] = g[1];

    groups.forEach(g => {
        const group: string[] = [];
        g.split(',').forEach(e => {
            const dashIndex = e.indexOf('-');
            if (dashIndex > 0) {
                const from = parseInt(e.substring(0, dashIndex)), to = parseInt(e.substr(dashIndex + 1));
                for (let i = from; i <= to; i++) group[group.length] = i.toString();
            } else {
                group[group.length] = e.trim();
            }
        });
        ret[ret.length] = group;
    });

    return ret;
}