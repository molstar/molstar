/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column, Table } from '../../../mol-data/db';
import { toTable } from '../../../mol-io/reader/cif/schema';
import { mmCIF_residueId_schema } from '../../../mol-io/reader/cif/schema/mmcif-extras';
import { CifWriter } from '../../../mol-io/writer/cif';
import { Model, ResidueIndex, Unit, IndexedCustomProperty } from '../../../mol-model/structure';
import { residueIdFields } from '../../../mol-model/structure/export/categories/atom_site';
import { StructureElement, CifExportContext, Structure } from '../../../mol-model/structure/structure';
import { CustomPropSymbol } from '../../../mol-script/language/symbol';
import Type from '../../../mol-script/language/type';
import { QuerySymbolRuntime } from '../../../mol-script/runtime/query/compiler';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { arraySetAdd } from '../../../mol-util/array';
import { MmcifFormat } from '../../../mol-model-formats/structure/mmcif';
import { PropertyWrapper } from '../../../mol-model-props/common/wrapper';
import { CustomProperty } from '../../../mol-model-props/common/custom-property';
import { CustomModelProperty } from '../../../mol-model-props/common/custom-model-property';
import { Asset } from '../../../mol-util/assets';
import { CustomPropertyDescriptor } from '../../../mol-model/custom-property';

export { StructureQualityReport };

type StructureQualityReport = PropertyWrapper<{
    issues: IndexedCustomProperty.Residue<string[]>,
    issueTypes: string[]
}| undefined>

namespace StructureQualityReport {
    export const DefaultServerUrl = 'https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry/';
    export function getEntryUrl(pdbId: string, serverUrl: string) {
        return `${serverUrl}/${pdbId.toLowerCase()}`;
    }

    export function isApplicable(model?: Model): boolean {
        return !!model && Model.isFromPdbArchive(model);
    }

    export const Schema = {
        pdbe_structure_quality_report: {
            updated_datetime_utc: Column.Schema.str
        },
        pdbe_structure_quality_report_issues: {
            id: Column.Schema.int,
            ...mmCIF_residueId_schema,
            pdbx_PDB_model_num: Column.Schema.int,
            issue_type_group_id: Column.Schema.int
        },
        pdbe_structure_quality_report_issue_types: {
            group_id: Column.Schema.int,
            issue_type: Column.Schema.str
        }
    };
    export type Schema = typeof Schema

    export function fromJson(model: Model, data: any) {
        const info = PropertyWrapper.createInfo();
        const issueMap = createIssueMapFromJson(model, data);
        return { info, data: issueMap };
    }

    export async function fromServer(ctx: CustomProperty.Context, model: Model, props: StructureQualityReportProps): Promise<CustomProperty.Data<StructureQualityReport>> {
        const url = Asset.getUrlAsset(ctx.assetManager, getEntryUrl(model.entryId, props.serverUrl));
        const json = await ctx.assetManager.resolve(url, 'json').runInContext(ctx.runtime);
        const data = json.data[model.entryId.toLowerCase()];
        if (!data) throw new Error('missing data');
        return { value: fromJson(model, data), assets: [json] };
    }

    export function fromCif(ctx: CustomProperty.Context, model: Model, props: StructureQualityReportProps): StructureQualityReport | undefined {
        let info = PropertyWrapper.tryGetInfoFromCif('pdbe_structure_quality_report', model);
        if (!info) return;
        const data = getCifData(model);
        const issueMap = createIssueMapFromCif(model, data.residues, data.groups);
        return { info, data: issueMap };
    }

    export async function fromCifOrServer(ctx: CustomProperty.Context, model: Model, props: StructureQualityReportProps): Promise<CustomProperty.Data<StructureQualityReport>> {
        const cif = fromCif(ctx, model, props);
        return cif ? { value: cif } : fromServer(ctx, model, props);
    }

    const _emptyArray: string[] = [];
    export function getIssues(e: StructureElement.Location) {
        if (!Unit.isAtomic(e.unit)) return _emptyArray;
        const prop = StructureQualityReportProvider.get(e.unit.model).value;
        if (!prop || !prop.data) return _emptyArray;
        const rI = e.unit.residueIndex[e.element];
        return prop.data.issues.has(rI) ? prop.data.issues.get(rI)! : _emptyArray;
    }

    export function getIssueTypes(structure?: Structure) {
        if (!structure) return _emptyArray;
        const prop = StructureQualityReportProvider.get(structure.models[0]).value;
        if (!prop || !prop.data) return _emptyArray;
        return prop.data.issueTypes;
    }

    function getCifData(model: Model) {
        if (!MmcifFormat.is(model.sourceData)) throw new Error('Data format must be mmCIF.');
        return {
            residues: toTable(Schema.pdbe_structure_quality_report_issues, model.sourceData.data.frame.categories.pdbe_structure_quality_report_issues),
            groups: toTable(Schema.pdbe_structure_quality_report_issue_types, model.sourceData.data.frame.categories.pdbe_structure_quality_report_issue_types),
        };
    }
}

export const StructureQualityReportParams = {
    serverUrl: PD.Text(StructureQualityReport.DefaultServerUrl, { description: 'JSON API Server URL' })
};
export type StructureQualityReportParams = typeof StructureQualityReportParams
export type StructureQualityReportProps = PD.Values<StructureQualityReportParams>

export const StructureQualityReportProvider: CustomModelProperty.Provider<StructureQualityReportParams, StructureQualityReport> = CustomModelProperty.createProvider({
    label: 'Structure Quality Report',
    descriptor: CustomPropertyDescriptor<ReportExportContext, any>({
        name: 'pdbe_structure_quality_report',
        cifExport: {
            prefix: 'pdbe',
            context(ctx: CifExportContext) {
                return createExportContext(ctx);
            },
            categories: [
                PropertyWrapper.defaultInfoCategory<ReportExportContext>('pdbe_structure_quality_report', ctx => ctx.info),
                {
                    name: 'pdbe_structure_quality_report_issues',
                    instance(ctx: ReportExportContext) {
                        return {
                            fields: _structure_quality_report_issues_fields,
                            source: ctx.models.map(data => ({ data, rowCount: data.elements.length }))
                        };
                    }
                }, {
                    name: 'pdbe_structure_quality_report_issue_types',
                    instance(ctx: ReportExportContext) {
                        return CifWriter.Category.ofTable(ctx.issueTypes);
                    }
                }]
        },
        symbols: {
            issueCount: QuerySymbolRuntime.Dynamic(CustomPropSymbol('pdbe', 'structure-quality.issue-count', Type.Num),
                ctx => StructureQualityReport.getIssues(ctx.element).length),
            // TODO: add (hasIssue :: IssueType(extends string) -> boolean) symbol
        }
    }),
    type: 'static',
    defaultParams: StructureQualityReportParams,
    getParams: (data: Model) => StructureQualityReportParams,
    isApplicable: (data: Model) => StructureQualityReport.isApplicable(data),
    obtain: async (ctx: CustomProperty.Context, data: Model, props: Partial<StructureQualityReportProps>) => {
        const p = { ...PD.getDefaultValues(StructureQualityReportParams), ...props };
        return await StructureQualityReport.fromCifOrServer(ctx, data, p);
    }
});

const _structure_quality_report_issues_fields = CifWriter.fields<number, ReportExportContext['models'][0]>()
    .index('id')
    .many(residueIdFields((i, d) => d.elements[i], { includeModelNum: true }))
    .int('issue_type_group_id', (i, d) => d.groupId[i])
    .getFields();

interface ReportExportContext {
    models: {
        elements: StructureElement.Location[],
        groupId: number[]
    }[],
    info: PropertyWrapper.Info,
    issueTypes: Table<StructureQualityReport.Schema['pdbe_structure_quality_report_issue_types']>,
}

function createExportContext(ctx: CifExportContext): ReportExportContext {
    const groupMap = new Map<string, number>();
    const models: ReportExportContext['models'] = [];
    const group_id: number[] = [], issue_type: string[] = [];
    let info: PropertyWrapper.Info = PropertyWrapper.createInfo();

    for (const s of ctx.structures) {
        const prop = StructureQualityReportProvider.get(s.model).value;
        if (prop) info = prop.info;
        if (!prop || !prop.data) continue;

        const { elements, property } = prop.data.issues.getElements(s);
        if (elements.length === 0) continue;

        const elementGroupId: number[] = [];
        for (let i = 0; i < elements.length; i++) {
            const issues = property(i);
            const key = issues.join(',');
            if (!groupMap.has(key)) {
                const idx = groupMap.size + 1;
                groupMap.set(key, idx);
                for (const issue of issues) {
                    group_id.push(idx);
                    issue_type.push(issue);
                }
            }
            elementGroupId[i] = groupMap.get(key)!;
        }
        models.push({ elements, groupId: elementGroupId });
    }

    return {
        info,
        models,
        issueTypes: Table.ofArrays(StructureQualityReport.Schema.pdbe_structure_quality_report_issue_types, { group_id, issue_type })
    };
}

function createIssueMapFromJson(modelData: Model, data: any): StructureQualityReport['data'] | undefined {
    const ret = new Map<ResidueIndex, string[]>();
    if (!data.molecules) return;

    const issueTypes: string[] = [];

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

                    for (const t of residue.outlier_types) {
                        arraySetAdd(issueTypes, t);
                    }
                }
            }
        }
    }

    return {
        issues: IndexedCustomProperty.fromResidueMap(ret),
        issueTypes
    };
}

function createIssueMapFromCif(modelData: Model,
    residueData: Table<typeof StructureQualityReport.Schema.pdbe_structure_quality_report_issues>,
    groupData: Table<typeof StructureQualityReport.Schema.pdbe_structure_quality_report_issue_types>): StructureQualityReport['data'] | undefined {

    const ret = new Map<ResidueIndex, string[]>();
    const { label_entity_id, label_asym_id, auth_seq_id, pdbx_PDB_ins_code, issue_type_group_id, pdbx_PDB_model_num, _rowCount } = residueData;

    const groups = parseIssueTypes(groupData);

    for (let i = 0; i < _rowCount; i++) {
        if (pdbx_PDB_model_num.value(i) !== modelData.modelNum) continue;
        const idx = modelData.atomicHierarchy.index.findResidue(label_entity_id.value(i), label_asym_id.value(i), auth_seq_id.value(i), pdbx_PDB_ins_code.value(i));
        ret.set(idx, groups.get(issue_type_group_id.value(i))!);
    }

    const issueTypes: string[] = [];
    groups.forEach(issues => {
        for (const t of issues) {
            arraySetAdd(issueTypes, t);
        }
    });

    return {
        issues: IndexedCustomProperty.fromResidueMap(ret),
        issueTypes
    };
}

function parseIssueTypes(groupData: Table<typeof StructureQualityReport.Schema.pdbe_structure_quality_report_issue_types>): Map<number, string[]> {
    const ret = new Map<number, string[]>();
    const { group_id, issue_type } = groupData;
    for (let i = 0; i < groupData._rowCount; i++) {
        let group: string[];
        const id = group_id.value(i);
        if (ret.has(id)) group = ret.get(id)!;
        else {
            group = [];
            ret.set(id, group);
        }
        group.push(issue_type.value(i));
    }
    return ret;
}