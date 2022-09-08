/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { DataLocation } from '../../../mol-model/location';
import { DataLoci } from '../../../mol-model/loci';
import { confalPyramidLabel } from './behavior';

export namespace ConfalPyramidsTypes {
    export const DataTag = 'dnatco-confal-half-pyramid';

    export type Step = {
        PDB_model_number: number,
        name: string,
        auth_asym_id_1: string,
        auth_seq_id_1: number,
        label_comp_id_1: string,
        label_alt_id_1: string,
        PDB_ins_code_1: string,
        auth_asym_id_2: string,
        auth_seq_id_2: number,
        label_comp_id_2: string,
        label_alt_id_2: string,
        PDB_ins_code_2: string,
        confal_score: number,
        NtC: string,
        rmsd: number,
    }

    export type MappedChains = Map<string, MappedResidues>;
    export type MappedResidues = Map<number, number[]>;

    export interface Steps {
        steps: Array<Step>,
        mapping: MappedChains[],
    }

    export interface HalfPyramid {
        step: Step,
        isLower: boolean,
    }

    export interface Location extends DataLocation<HalfPyramid, {}> {}

    export function Location(step: Step, isLower: boolean) {
        return DataLocation(DataTag, { step, isLower }, {});
    }

    export function isLocation(x: any): x is Location {
        return !!x && x.kind === 'data-location' && x.tag === DataTag;
    }

    export interface Loci extends DataLoci<HalfPyramid, {}> {}

    export function Loci(data: HalfPyramid, elements: ReadonlyArray<{}>): Loci {
        return DataLoci(DataTag, data, elements, undefined, () => confalPyramidLabel(data));
    }

    export function isLoci(x: any): x is Loci {
        return !!x && x.kind === 'data-loci' && x.tag === DataTag;
    }
}
