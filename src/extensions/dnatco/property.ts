/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { DnatcoTypes } from './types';
import { Column, Table } from '../../mol-data/db';
import { toTable } from '../../mol-io/reader/cif/schema';
import { Model } from '../../mol-model/structure';
import { CustomProperty } from '../../mol-model-props/common/custom-property';
import { PropertyWrapper } from '../../mol-model-props/common/wrapper';
import { MmcifFormat } from '../../mol-model-formats/structure/mmcif';
import { ParamDefinition as PD } from '../../mol-util/param-definition';

export type DnatcoSteps = PropertyWrapper<DnatcoTypes.Steps | undefined>;

export const DnatcoParams = {};
export type DnatcoParams = typeof DnatcoParams;
export type DnatcoProps = PD.Values<DnatcoParams>;

export namespace Dnatco {
    export const Schema = {
        ndb_struct_ntc_step: {
            id: Column.Schema.int,
            name: Column.Schema.str,
            PDB_model_number: Column.Schema.int,
            label_entity_id_1: Column.Schema.int,
            label_asym_id_1: Column.Schema.str,
            label_seq_id_1: Column.Schema.int,
            label_comp_id_1: Column.Schema.str,
            label_alt_id_1: Column.Schema.str,
            label_entity_id_2: Column.Schema.int,
            label_asym_id_2: Column.Schema.str,
            label_seq_id_2: Column.Schema.int,
            label_comp_id_2: Column.Schema.str,
            label_alt_id_2: Column.Schema.str,
            auth_asym_id_1: Column.Schema.str,
            auth_seq_id_1: Column.Schema.int,
            auth_asym_id_2: Column.Schema.str,
            auth_seq_id_2: Column.Schema.int,
            PDB_ins_code_1: Column.Schema.str,
            PDB_ins_code_2: Column.Schema.str,
        },
        ndb_struct_ntc_step_summary: {
            step_id: Column.Schema.int,
            assigned_CANA: Column.Schema.str,
            assigned_NtC: Column.Schema.str,
            confal_score: Column.Schema.int,
            euclidean_distance_NtC_ideal: Column.Schema.float,
            cartesian_rmsd_closest_NtC_representative: Column.Schema.float,
            closest_CANA: Column.Schema.str,
            closest_NtC: Column.Schema.str,
            closest_step_golden: Column.Schema.str
        }
    };
    export type Schema = typeof Schema;

    export function getStepsFromCif(
        model: Model,
        cifSteps: Table<typeof Dnatco.Schema.ndb_struct_ntc_step>,
        stepsSummary: StepsSummaryTable
    ): DnatcoTypes.Steps {
        const steps = new Array<DnatcoTypes.Step>();
        const mapping = new Array<DnatcoTypes.MappedChains>();

        const {
            id, PDB_model_number, name,
            auth_asym_id_1, auth_seq_id_1, label_comp_id_1, label_alt_id_1, PDB_ins_code_1,
            auth_asym_id_2, auth_seq_id_2, label_comp_id_2, label_alt_id_2, PDB_ins_code_2,
            _rowCount
        } = cifSteps;

        if (_rowCount !== stepsSummary._rowCount) throw new Error('Inconsistent mmCIF data');

        for (let i = 0; i < _rowCount; i++) {
            const {
                NtC,
                confal_score,
                rmsd
            } = getSummaryData(id.value(i), i, stepsSummary);
            const modelNum = PDB_model_number.value(i);
            const chainId = auth_asym_id_1.value(i);
            const seqId = auth_seq_id_1.value(i);
            const modelIdx = modelNum - 1;

            if (mapping.length <= modelIdx || !mapping[modelIdx])
                mapping[modelIdx] = new Map<string, DnatcoTypes.MappedResidues>();

            const step = {
                PDB_model_number: modelNum,
                name: name.value(i),
                auth_asym_id_1: chainId,
                auth_seq_id_1: seqId,
                label_comp_id_1: label_comp_id_1.value(i),
                label_alt_id_1: label_alt_id_1.value(i),
                PDB_ins_code_1: PDB_ins_code_1.value(i),
                auth_asym_id_2: auth_asym_id_2.value(i),
                auth_seq_id_2: auth_seq_id_2.value(i),
                label_comp_id_2: label_comp_id_2.value(i),
                label_alt_id_2: label_alt_id_2.value(i),
                PDB_ins_code_2: PDB_ins_code_2.value(i),
                confal_score,
                NtC,
                rmsd,
            };

            steps.push(step);

            const mappedChains = mapping[modelIdx];
            const residuesOnChain = mappedChains.get(chainId) ?? new Map<number, number[]>();
            const stepsForResidue = residuesOnChain.get(seqId) ?? [];
            stepsForResidue.push(steps.length - 1);

            residuesOnChain.set(seqId, stepsForResidue);
            mappedChains.set(chainId, residuesOnChain);
            mapping[modelIdx] = mappedChains;
        }

        return { steps, mapping };
    }

    export async function fromCif(ctx: CustomProperty.Context, model: Model, props: DnatcoProps): Promise<CustomProperty.Data<DnatcoSteps>> {
        const info = PropertyWrapper.createInfo();
        const data = getCifData(model);
        if (data === undefined) return { value: { info, data: undefined } };

        const fromCif = getStepsFromCif(model, data.steps, data.stepsSummary);
        return { value: { info, data: fromCif } };
    }

    export function getCifData(model: Model) {
        if (!MmcifFormat.is(model.sourceData)) throw new Error('Data format must be mmCIF');
        if (!hasNdbStructNtcCategories(model)) return undefined;
        return {
            steps: toTable(Schema.ndb_struct_ntc_step, model.sourceData.data.frame.categories.ndb_struct_ntc_step),
            stepsSummary: toTable(Schema.ndb_struct_ntc_step_summary, model.sourceData.data.frame.categories.ndb_struct_ntc_step_summary)
        };
    }

    function hasNdbStructNtcCategories(model: Model): boolean {
        if (!MmcifFormat.is(model.sourceData)) return false;
        const names = (model.sourceData).data.frame.categoryNames;
        return names.includes('ndb_struct_ntc_step') && names.includes('ndb_struct_ntc_step_summary');
    }

    export function isApplicable(model?: Model): boolean {
        return !!model && hasNdbStructNtcCategories(model);
    }
}

type StepsSummaryTable = Table<typeof Dnatco.Schema.ndb_struct_ntc_step_summary>;

function getSummaryData(id: number, i: number, stepsSummary: StepsSummaryTable) {
    const {
        step_id,
        confal_score,
        assigned_NtC,
        cartesian_rmsd_closest_NtC_representative,
    } = stepsSummary;

    // Assume that step_ids in ntc_step_summary are in the same order as steps in ntc_step
    for (let j = i; j < stepsSummary._rowCount; j++) {
        if (id === step_id.value(j)) return { NtC: assigned_NtC.value(j), confal_score: confal_score.value(j), rmsd: cartesian_rmsd_closest_NtC_representative.value(j) };
    }
    // Safety net for cases where the previous assumption is not met
    for (let j = 0; j < i; j++) {
        if (id === step_id.value(j)) return { NtC: assigned_NtC.value(j), confal_score: confal_score.value(j), rmsd: cartesian_rmsd_closest_NtC_representative.value(j) };
    }
    throw new Error('Inconsistent mmCIF data');
}
