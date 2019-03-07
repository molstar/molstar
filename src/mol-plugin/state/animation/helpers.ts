/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */


import { SymmetryOperator } from 'mol-math/geometry';
import { Mat4, Vec3 } from 'mol-math/linear-algebra';
import { Structure, StructureSelection, QueryContext } from 'mol-model/structure';
import { StructureUnitTransforms } from 'mol-model/structure/structure/util/unit-transforms';
import { Color } from 'mol-util/color';
import { Overpaint } from 'mol-theme/overpaint';
import { parseMolScript } from 'mol-script/language/parser';
import { transpileMolScript } from 'mol-script/script/mol-script/symbols';
import { compile } from 'mol-script/runtime/query/compiler';

const _unwindMatrix = Mat4.zero();
export function unwindStructureAssembly(structure: Structure, unitTransforms: StructureUnitTransforms, t: number) {
    for (let i = 0, _i = structure.units.length; i < _i; i++) {
        const u = structure.units[i];
        SymmetryOperator.lerpFromIdentity(_unwindMatrix, u.conformation.operator, t);
        unitTransforms.setTransform(_unwindMatrix, u);
    }
}

const _centerVec = Vec3.zero(), _transVec = Vec3.zero(), _transMat = Mat4.zero();
export function explodeStructure(structure: Structure, unitTransforms: StructureUnitTransforms, t: number) {
    const boundary = structure.boundary.sphere;
    const d = boundary.radius * t;

    for (let i = 0, _i = structure.units.length; i < _i; i++) {
        const u = structure.units[i];
        Vec3.transformMat4(_centerVec, u.lookup3d.boundary.sphere.center, u.conformation.operator.matrix);
        Vec3.sub(_transVec, _centerVec, boundary.center);
        Vec3.setMagnitude(_transVec, _transVec, d);
        Mat4.fromTranslation(_transMat, _transVec);

        unitTransforms.setTransform(_transMat, u);
    }
}

type ScriptLayers = { script: { language: string, expression: string }, color: Color }[]
export function getStructureOverpaint(structure: Structure, scriptLayers: ScriptLayers, alpha: number): Overpaint {
    const layers: Overpaint.Layer[] = []
    for (let i = 0, il = scriptLayers.length; i < il; ++i) {
        const { script, color } = scriptLayers[i]
        const parsed = parseMolScript(script.expression)
        if (parsed.length === 0) throw new Error('No query')
        const query = transpileMolScript(parsed[0])

        const compiled = compile<StructureSelection>(query)
        const result = compiled(new QueryContext(structure))
        const loci = StructureSelection.toLoci2(result)

        layers.push({ loci, color })
    }
    return { layers, alpha }
}