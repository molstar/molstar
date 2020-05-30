/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SymmetryOperator } from '../../../../mol-math/geometry/symmetry-operator';
import { arrayFind } from '../../../../mol-data/util';
import { StructureQuery } from '../../query';
import { Model } from '../../model';
import { Spacegroup } from '../../../../mol-math/geometry';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { ModelSymmetry } from '../../../../mol-model-formats/structure/property/symmetry';
import { radToDeg } from '../../../../mol-math/misc';

/** Determine an atom set and a list of operators that should be applied to that set  */
export interface OperatorGroup {
    readonly selector: StructureQuery,
    readonly operators: ReadonlyArray<SymmetryOperator>
}

export type OperatorGroups = ReadonlyArray<OperatorGroup>

export class Assembly {
    readonly id: string;
    readonly details: string;

    private _operators: OperatorGroups;
    get operatorGroups(): OperatorGroups {
        if (this._operators) return this._operators;
        this._operators = this.operatorsProvider();
        return this._operators;
    }

    constructor(id: string, details: string, private operatorsProvider: () => OperatorGroups) {
        this.id = id;
        this.details = details;
    }
}

export namespace Assembly {
    export function create(id: string, details: string, operatorsProvider: () => OperatorGroups): Assembly {
        return new Assembly(id, details, operatorsProvider);
    }
}

interface Symmetry {
    readonly assemblies: ReadonlyArray<Assembly>,
    readonly spacegroup: Spacegroup,
    readonly isNonStandardCrystalFrame: boolean,
    readonly ncsOperators?: ReadonlyArray<SymmetryOperator>,

    /**
     * optionally cached operators from [-3, -3, -3] to [3, 3, 3]
     * around reference point `ref` in fractional coordinates
     */
    _operators_333?: {
        ref: Vec3,
        operators: SymmetryOperator[]
    }
}

namespace Symmetry {
    export const Default: Symmetry = { assemblies: [], spacegroup: Spacegroup.ZeroP1, isNonStandardCrystalFrame: false };

    export function findAssembly(model: Model, id: string): Assembly | undefined {
        const _id = id.toLocaleLowerCase();
        const symmetry = ModelSymmetry.Provider.get(model);
        return symmetry ? arrayFind(symmetry.assemblies, a => a.id.toLowerCase() === _id) : undefined;
    }

    export function getUnitcellLabel(symmetry: Symmetry) {
        const { cell, name, num } = symmetry.spacegroup;
        const { size, anglesInRadians } = cell;
        const a = size[0].toFixed(2);
        const b = size[1].toFixed(2);
        const c = size[2].toFixed(2);
        const alpha = radToDeg(anglesInRadians[0]).toFixed(2);
        const beta = radToDeg(anglesInRadians[1]).toFixed(2);
        const gamma = radToDeg(anglesInRadians[2]).toFixed(2);
        const label: string[] = [];
        // name
        label.push(`Unit Cell <b>${name}</b> #${num}`);
        // sizes
        label.push(`${a}\u00D7${b}\u00D7${c} \u212B`);
        // angles
        label.push(`\u03b1=${alpha}\u00B0 \u03b2=${beta}\u00B0 \u03b3=${gamma}\u00B0`);
        return label.join(' | ');
    }
}

export { Symmetry };