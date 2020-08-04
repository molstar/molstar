/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { Column, Table } from '../../mol-data/db';
import { toTable } from '../../mol-io/reader/cif/schema';
import { Model } from '../../mol-model/structure';
import { MmcifFormat } from '../../mol-model-formats/structure/mmcif';
import { Color } from '../../mol-util/color';

export namespace DnatcoCommon {
    export const AColor = Color(0xFFC1C1);
    export const BColor = Color(0xC8CFFF);
    export const BIIColor = Color(0x0059DA);
    export const miBColor = Color(0x3BE8FB);
    export const ZColor = Color(0x01F60E);
    export const ICColor = Color(0xFA5CFB);
    export const OPNColor = Color(0xE90000);
    export const SYNColor = Color(0xFFFF01);
    export const NColor = Color(0xF2F2F2);

    export interface NtCObject {
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
        NtC: string
    }

    export const CifSchema = {
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

    export function getCifData(model: Model) {
        if (!MmcifFormat.is(model.sourceData)) throw new Error('Data format must be mmCIF');
        if (!hasNdbStructNtcCategories(model)) return undefined;
        return {
            steps: toTable(CifSchema.ndb_struct_ntc_step, model.sourceData.data.frame.categories.ndb_struct_ntc_step),
            stepsSummary: toTable(CifSchema.ndb_struct_ntc_step_summary, model.sourceData.data.frame.categories.ndb_struct_ntc_step_summary)
        };
    }

    export type StepsSummaryTable = Table<typeof CifSchema.ndb_struct_ntc_step_summary>;

    export function getNtCAndConfalScore(id: number, i: number, stepsSummary: StepsSummaryTable) {
        const { step_id, confal_score, assigned_NtC } = stepsSummary;

        // Assume that step_ids in ntc_step_summary are in the same order as steps in ntc_step
        for (let j = i; j < stepsSummary._rowCount; j++) {
            if (id === step_id.value(j)) return { _NtC: assigned_NtC.value(j), _confal_score: confal_score.value(j) };
        }
        // Safety net for cases where the previous assumption is not met
        for (let j = 0; j < i; j++) {
            if (id === step_id.value(j)) return { _NtC: assigned_NtC.value(j), _confal_score: confal_score.value(j) };
        }
        throw new Error('Inconsistent mmCIF data');
    }

    export function hasNdbStructNtcCategories(model: Model): boolean {
        if (!MmcifFormat.is(model.sourceData)) return false;
        const names = (model.sourceData).data.frame.categoryNames;
        return names.includes('ndb_struct_ntc_step') && names.includes('ndb_struct_ntc_step_summary');
    }

    export function isApplicable(model?: Model): boolean {
        return !!model && hasNdbStructNtcCategories(model);
    }
}
