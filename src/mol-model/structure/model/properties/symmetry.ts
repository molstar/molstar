/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import SymmetryOperator from 'mol-math/geometry/symmetry-operator'
import { Query } from '../../query'
import { Model } from '../../model'

export interface OperatorGroup {
    readonly query: Query,
    readonly operators: ReadonlyArray<SymmetryOperator>
}

export type OperatorGroups = ReadonlyArray<ReadonlyArray<OperatorGroup>>

export class Assembly {
    readonly id: string;
    readonly details: string;

    private _operators: OperatorGroups;
    get operators(): OperatorGroups {
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
    export function create(id: string, details: string, operatorsProvider: () => OperatorGroups) {
        return new Assembly(id, details, operatorsProvider);
    }
}

interface Symmetry {
    readonly assemblies: ReadonlyArray<Assembly>,
}

namespace Symmetry {
    export const Empty: Symmetry = { assemblies: [] };

    export function findAssembly(model: Model): Assembly | undefined {
        return void 0;
    }
}

export default Symmetry