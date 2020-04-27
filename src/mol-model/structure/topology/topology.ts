/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UUID } from '../../../mol-util';
import { Column } from '../../../mol-data/db';
import { BasicData } from '../../../mol-model-formats/structure/basic/schema';
import { ModelFormat } from '../../../mol-model-formats/format';

export { Topology };

interface Topology {
    readonly id: UUID
    readonly label: string

    readonly basic: BasicData
    readonly sourceData: ModelFormat

    readonly bonds: {
        readonly indexA: Column<number>,
        readonly indexB: Column<number>
        readonly order: Column<number>
    }

    // TODO
    // readonly angles: {
    //     readonly indexA: Column<number>
    //     readonly indexB: Column<number>
    //     readonly indexC: Column<number>
    // }

    // readonly dihedrals: {
    //     readonly indexA: Column<number>
    //     readonly indexB: Column<number>
    //     readonly indexC: Column<number>
    //     readonly indexD: Column<number>
    // }

    // readonly impropers: {
    //     readonly indexA: Column<number>
    //     readonly indexB: Column<number>
    //     readonly indexC: Column<number>
    //     readonly indexD: Column<number>
    // }

    // TODO: add forces for bonds/angles/dihedrals/impropers
}

namespace Topology {
    export function create(label: string, basic: BasicData, bonds: Topology['bonds'], format: ModelFormat): Topology {
        return {
            id: UUID.create22(),
            label,
            basic,
            sourceData: format,
            bonds
        };
    }
}