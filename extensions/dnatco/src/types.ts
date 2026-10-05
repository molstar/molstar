/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

export namespace DnatcoTypes {
    export const DataTag = 'dnatco-confal-half-step';

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

    export interface HalfStep {
        step: Step,
        isLower: boolean,
    }
}
