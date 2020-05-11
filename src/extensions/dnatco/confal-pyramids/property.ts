/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Michal Malý <michal.maly@ibt.cas.cz>
 * @author Jiří Černý <jiri.cerny@ibt.cas.cz>
 */

import { ConfalPyramidsTypes as CPT } from './types';
import { Column, Table } from '../../../mol-data/db';
import { toTable } from '../../../mol-io/reader/cif/schema';
import { CustomPropertyDescriptor } from '../../../mol-model/custom-property';
import { Model } from '../../../mol-model/structure';
import { CustomProperty } from '../../../mol-model-props/common/custom-property';
import { CustomModelProperty } from '../../../mol-model-props/common/custom-model-property';
import { PropertyWrapper } from '../../../mol-model-props/common/wrapper';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { MmcifFormat } from '../../../mol-model-formats/structure/mmcif';

export type ConfalPyramids = PropertyWrapper<CPT.PyramidsData | undefined >;

export namespace ConfalPyramids {
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

    export async function fromCif(ctx: CustomProperty.Context, model: Model, props: ConfalPyramidsProps): Promise<CustomProperty.Data<ConfalPyramids>> {
        const info = PropertyWrapper.createInfo();
        const data = getCifData(model);
        if (data === undefined) return { value: { info, data: undefined } };

        const fromCif = createPyramidsFromCif(model, data.steps, data.stepsSummary);
        return { value: { info, data: fromCif } };
    }

    function getCifData(model: Model) {
        if (!MmcifFormat.is(model.sourceData)) throw new Error('Data format must be mmCIF');
        if (!hasNdbStructNtcCategories(model)) return undefined;
        return {
            steps: toTable(Schema.ndb_struct_ntc_step, model.sourceData.data.frame.categories.ndb_struct_ntc_step),
            stepsSummary: toTable(Schema.ndb_struct_ntc_step_summary, model.sourceData.data.frame.categories.ndb_struct_ntc_step_summary)
        };
    }

    function hasNdbStructNtcCategories(model: Model): boolean {
        if (!MmcifFormat.is(model.sourceData)) throw new Error('Data format must be mmCIF');
        const names = (model.sourceData).data.frame.categoryNames;
        return names.includes('ndb_struct_ntc_step') && names.includes('ndb_struct_ntc_step_summary');
    }

    export function isApplicable(model?: Model): boolean {
        return !!model && hasNdbStructNtcCategories(model);
    }
}

export const ConfalPyramidsParams = {};
export type ConfalPyramidsParams = typeof ConfalPyramidsParams;
export type ConfalPyramidsProps = PD.Values<ConfalPyramidsParams>;

export const ConfalPyramidsProvider: CustomModelProperty.Provider<ConfalPyramidsParams, ConfalPyramids> = CustomModelProperty.createProvider({
    label: 'Confal Pyramids',
    descriptor: CustomPropertyDescriptor({
        name: 'confal_pyramids',
    }),
    type: 'static',
    defaultParams: ConfalPyramidsParams,
    getParams: (data: Model) => ConfalPyramidsParams,
    isApplicable: (data: Model) => ConfalPyramids.isApplicable(data),
    obtain: async (ctx: CustomProperty.Context, data: Model, props: Partial<ConfalPyramidsProps>) => {
        const p = { ...PD.getDefaultValues(ConfalPyramidsParams), ...props };
        return ConfalPyramids.fromCif(ctx, data, p);
    }
});

type StepsSummaryTable = Table<typeof ConfalPyramids.Schema.ndb_struct_ntc_step_summary>;

function createPyramidsFromCif(model: Model,
    steps: Table<typeof ConfalPyramids.Schema.ndb_struct_ntc_step>,
    stepsSummary: StepsSummaryTable): CPT.PyramidsData {
    const pyramids = new Array<CPT.Pyramid>();
    const names = new Map<string, number>();
    const locations = new Array<CPT.Location>();
    let hasMultipleModels = false;

    const {
        id, PDB_model_number, name,
        auth_asym_id_1, auth_seq_id_1, label_comp_id_1, label_alt_id_1, PDB_ins_code_1,
        auth_asym_id_2, auth_seq_id_2, label_comp_id_2, label_alt_id_2, PDB_ins_code_2,
        _rowCount } = steps;

    if (_rowCount !== stepsSummary._rowCount) throw new Error('Inconsistent mmCIF data');

    for (let i = 0; i < _rowCount; i++) {
        const model_num = PDB_model_number.value(i);
        if (model_num !== model.modelNum) {
            hasMultipleModels = true;
            continue; // We are only interested in data for the current model
        }

        const { _NtC, _confal_score } = getNtCAndConfalScore(id.value(i), i, stepsSummary);

        const pyramid = {
            PDB_model_number: model_num,
            name: name.value(i),
            auth_asym_id_1: auth_asym_id_1.value(i),
            auth_seq_id_1: auth_seq_id_1.value(i),
            label_comp_id_1: label_comp_id_1.value(i),
            label_alt_id_1: label_alt_id_1.value(i),
            PDB_ins_code_1: PDB_ins_code_1.value(i),
            auth_asym_id_2: auth_asym_id_2.value(i),
            auth_seq_id_2: auth_seq_id_2.value(i),
            label_comp_id_2: label_comp_id_2.value(i),
            label_alt_id_2: label_alt_id_2.value(i),
            PDB_ins_code_2: PDB_ins_code_2.value(i),
            confal_score: _confal_score,
            NtC: _NtC
        };

        pyramids.push(pyramid);
        names.set(pyramid.name, pyramids.length - 1);

        locations.push(CPT.Location(pyramid, false));
        locations.push(CPT.Location(pyramid, true));
    }

    return { pyramids, names, locations, hasMultipleModels };
}

function getNtCAndConfalScore(id: number, i: number, stepsSummary: StepsSummaryTable) {
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
