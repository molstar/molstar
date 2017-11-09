/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import SymmetryOperator from 'mol-math/geometry/symmetry-operator'
import { Query } from '../../query'
import { Model } from '../../model'

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

interface Symmetry {
    readonly assemblies: ReadonlyArray<Assembly>,
}

namespace Symmetry {
    export const Empty: Symmetry = { assemblies: [] };

    export function findAssembly(model: Model, id: string): Assembly | undefined {
        const { assemblies } = model.symmetry;
        const _id = id.toLocaleLowerCase();
        for (let i = 0; i < assemblies.length; i++) {
            if (assemblies[i].id.toLowerCase() === _id) return assemblies[i];
        }
        return void 0;
    }
}

export default Symmetry