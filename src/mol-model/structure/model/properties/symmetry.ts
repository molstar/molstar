/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SymmetryOperator } from 'mol-math/geometry/symmetry-operator'
import { arrayFind } from 'mol-data/util'
import { Query } from '../../query'
import { Model } from '../../model'
import { Spacegroup } from 'mol-math/geometry';

/** Determine an atom set and a list of operators that should be applied to that set  */
export interface OperatorGroup {
    readonly selector: Query,
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

interface ModelSymmetry {
    readonly assemblies: ReadonlyArray<Assembly>,
    readonly spacegroup: Spacegroup,
    readonly isNonStandardCrytalFrame: boolean,
    readonly ncsOperators?: ReadonlyArray<SymmetryOperator>,

    // optionally cached operators from [-3, -3, -3] to [3, 3, 3]
    _operators_333?: SymmetryOperator[]
}

namespace ModelSymmetry {
    export const Default: ModelSymmetry = { assemblies: [], spacegroup: Spacegroup.ZeroP1, isNonStandardCrytalFrame: false };

    export function findAssembly(model: Model, id: string): Assembly | undefined {
        const _id = id.toLocaleLowerCase();
        return arrayFind(model.symmetry.assemblies, a => a.id.toLowerCase() === _id);
    }
}

export { ModelSymmetry }