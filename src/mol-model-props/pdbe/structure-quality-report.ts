/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifWriter } from 'mol-io/writer/cif';
import { Model, ModelPropertyDescriptor, ResidueIndex, Unit, ResidueCustomProperty, StructureProperties as P } from 'mol-model/structure';
import { residueIdFields } from 'mol-model/structure/export/categories/atom_site';
import CifField = CifWriter.Field;
import { mmCIF_residueId_schema } from 'mol-io/reader/cif/schema/mmcif-extras';
import { Column, Table } from 'mol-data/db';
import { toTable } from 'mol-io/reader/cif/schema';
import { StructureElement } from 'mol-model/structure/structure';


import { QuerySymbolRuntime } from 'mol-script/runtime/query/compiler';
import { CustomPropSymbol } from 'mol-script/language/symbol';
import Type from 'mol-script/language/type';

type IssueMap = ResidueCustomProperty<string[]>

const _Descriptor = ModelPropertyDescriptor({
    isStatic: false,
    name: 'structure_quality_report',
    cifExport: {
        prefix: 'pdbe',
        categories: [{
            name: 'pdbe_structure_quality_report',
            instance(ctx) {
                const issues = StructureQualityReport.get(ctx.model);
                if (typeof ctx.globalCache.pdbe_structure_quality_report !== 'undefined' && ctx.globalCache.pdbe_structure_quality_report !== ctx.model.modelNum) return CifWriter.Category.Empty;
                ctx.globalCache.pdbe_structure_quality_report = ctx.model.modelNum;
                return { fields: _structure_quality_report_fields, rowCount: 1, data: issues ? issues.updated : 'n/a' }
            }
        }, {
            name: 'pdbe_structure_quality_report_issues',
            instance(ctx) {
                const issues = StructureQualityReport.get(ctx.model);
                if (!issues || !issues.map) return CifWriter.Category.Empty;
                return ResidueCustomProperty.createCifCategory(ctx, issues.map, _structure_quality_report_issues_fields);
            }
        }]
    },
    symbols: {
        issueCount: QuerySymbolRuntime.Dynamic(CustomPropSymbol('pdbe', 'structure-quality.issue-count', Type.Num),
            ctx => StructureQualityReport.getIssues(ctx.element).length),
        // TODO: add (hasIssue :: IssueType(extends string) -> boolean) symbol
    }
})

type ExportCtx = ResidueCustomProperty.ExportCtx<string[]>
const _structure_quality_report_issues_fields: CifField<number, ExportCtx>[] = [
    CifField.index('id'),
    ...residueIdFields<number, ExportCtx>((i, d) => d.elements[i]),
    CifField.int<number, ExportCtx>('pdbx_PDB_model_num', (i, d) => P.unit.model_num(d.elements[i])),
    CifField.str<number, ExportCtx>('issues', (i, d) => d.property(i).join(','))
];

const _structure_quality_report_fields: CifField<number, string>[] = [
    CifField.str('updated_datetime_utc', (_, date) => date)
];

function createIssueMapFromJson(modelData: Model, data: any): IssueMap | undefined {
    const ret = new Map<ResidueIndex, string[]>();
    if (!data.molecules) return;

    for (const entity of data.molecules) {
        const entity_id = entity.entity_id.toString();
        for (const chain of entity.chains) {
            const asym_id = chain.struct_asym_id.toString();
            for (const model of chain.models) {
                const model_id = model.model_id.toString();
                if (+model_id !== modelData.modelNum) continue;

                for (const residue of model.residues) {
                    const auth_seq_id = residue.author_residue_number, ins_code = residue.author_insertion_code || '';
                    const idx = modelData.atomicHierarchy.index.findResidue(entity_id, asym_id, auth_seq_id, ins_code);
                    ret.set(idx, residue.outlier_types);
                }
            }
        }
    }

    return ResidueCustomProperty.fromMap(ret, Unit.Kind.Atomic);
}

function createIssueMapFromCif(modelData: Model, data: Table<typeof StructureQualityReport.Schema.pdbe_structure_quality_report_issues>): IssueMap | undefined {
    const ret = new Map<ResidueIndex, string[]>();
    const { label_entity_id, label_asym_id, auth_seq_id, pdbx_PDB_ins_code, issues, pdbx_PDB_model_num, _rowCount } = data;

    for (let i = 0; i < _rowCount; i++) {
        if (pdbx_PDB_model_num.value(i) !== modelData.modelNum) continue;
        const idx = modelData.atomicHierarchy.index.findResidue(label_entity_id.value(i), label_asym_id.value(i), auth_seq_id.value(i), pdbx_PDB_ins_code.value(i));
        ret.set(idx, issues.value(i));
    }

    return ResidueCustomProperty.fromMap(ret, Unit.Kind.Atomic);
}

export namespace StructureQualityReport {
    export interface Data {
        updated: string,
        map: IssueMap | undefined
    }

    export const Descriptor = _Descriptor;

    export const Schema = {
        pdbe_structure_quality_report: {
            updated_datetime_utc: Column.Schema.str
        },
        pdbe_structure_quality_report_issues: {
            id: Column.Schema.int,
            ...mmCIF_residueId_schema,
            pdbx_PDB_model_num: Column.Schema.int,
            issues: Column.Schema.List(',', x => x)
        }
    }

    export async function attachFromCifOrApi(model: Model, params: {
        // provide JSON from api
        PDBe_apiSourceJson?: (model: Model) => Promise<any>
    }) {
        if (get(model)) return true;

        let issueMap, updated = `${new Date().toISOString().replace(/T/, ' ').replace(/\..+/, '')}`;
        if (model.sourceData.kind === 'mmCIF' && model.sourceData.frame.categoryNames.includes('pdbe_structure_quality_report')) {
            const data = toTable(Schema.pdbe_structure_quality_report_issues, model.sourceData.frame.categories.pdbe_structure_quality_report_issues);
            const f = model.sourceData.frame.categories['pdbe_structure_quality_report'].getField('updated_datetime_utc');
            updated = f ? f.str(0) : updated;
            issueMap = createIssueMapFromCif(model, data);
        } else if (params.PDBe_apiSourceJson) {
            const data = await params.PDBe_apiSourceJson(model);
            if (!data) return false;
            issueMap = createIssueMapFromJson(model, data);
        } else {
            return false;
        }

        model.customProperties.add(Descriptor);
        (model._dynamicPropertyData.__StructureQualityReport__ as Data) = {
            updated,
            map: issueMap
        };
        return true;
    }

    export function get(model: Model): Data | undefined {
        return model._dynamicPropertyData.__StructureQualityReport__;
    }

    const _emptyArray: string[] = [];
    export function getIssues(e: StructureElement) {
        if (!Unit.isAtomic(e.unit)) return _emptyArray;
        const issues = StructureQualityReport.get(e.unit.model);
        if (!issues || !issues.map) return _emptyArray;
        const rI = e.unit.residueIndex[e.element];
        return issues.map.has(rI) ? issues.map.get(rI)! : _emptyArray;
    }
}