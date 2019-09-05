/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, Table } from '../../mol-data/db';
import { toTable } from '../../mol-io/reader/cif/schema';
import { mmCIF_residueId_schema } from '../../mol-io/reader/cif/schema/mmcif-extras';
import { CifWriter } from '../../mol-io/writer/cif';
import { Model, CustomPropertyDescriptor, ResidueIndex, Unit, IndexedCustomProperty } from '../../mol-model/structure';
import { residueIdFields } from '../../mol-model/structure/export/categories/atom_site';
import { StructureElement, CifExportContext } from '../../mol-model/structure/structure';
import { CustomPropSymbol } from '../../mol-script/language/symbol';
import Type from '../../mol-script/language/type';
import { QuerySymbolRuntime } from '../../mol-script/runtime/query/compiler';
import { PropertyWrapper } from '../common/wrapper';
import { Task } from '../../mol-task';

export namespace StructureQualityReport {
    export type IssueMap = IndexedCustomProperty.Residue<string[]>
    export type Property = PropertyWrapper<IssueMap | undefined>

    export function get(model: Model): Property | undefined {
        // must be defined before the descriptor so it's not undefined.
        return model._dynamicPropertyData.__PDBeStructureQualityReport__;
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

    export const Descriptor = CustomPropertyDescriptor({
        isStatic: false,
        name: 'pdbe_structure_quality_report',
        cifExport: {
            prefix: 'pdbe',
            context(ctx) {
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
                        }
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
    });

    function getCifData(model: Model) {
        if (model.sourceData.kind !== 'mmCIF') throw new Error('Data format must be mmCIF.');
        return {
            residues: toTable(Schema.pdbe_structure_quality_report_issues, model.sourceData.frame.categories.pdbe_structure_quality_report_issues),
            groups: toTable(Schema.pdbe_structure_quality_report_issue_types, model.sourceData.frame.categories.pdbe_structure_quality_report_issue_types),
        }
    }

    export function createAttachTask(mapUrl: (model: Model) => string, fetch: import('../../mol-util/data-source').AjaxTask) {
        return (model: Model) => Task.create('PDBe Structure Quality Report', async ctx => {
            if (get(model)) return true;

            let issueMap: IssueMap | undefined;
            let info;
            // TODO: return from CIF support once the data is recomputed
            //  = PropertyWrapper.tryGetInfoFromCif('pdbe_structure_quality_report', model);
            // if (info) {
            //     const data = getCifData(model);
            //     issueMap = createIssueMapFromCif(model, data.residues, data.groups);
            // } else
            {
                const url = mapUrl(model);
                const dataStr = await fetch({ url }).runInContext(ctx) as string;
                const data = JSON.parse(dataStr)[model.entryId.toLowerCase()];
                if (!data) return false;
                info = PropertyWrapper.createInfo();
                issueMap = createIssueMapFromJson(model, data);
            }

            model.customProperties.add(Descriptor);
            set(model, { info, data: issueMap });
            return false;
        });
    }

    export async function attachFromCifOrApi(model: Model, params: {
        // optional JSON source
        PDBe_apiSourceJson?: (model: Model) => Promise<any>
    }) {
        if (get(model)) return true;

        let issueMap: IssueMap | undefined;
        let info = PropertyWrapper.tryGetInfoFromCif('pdbe_structure_quality_report', model);
        if (info) {
            const data = getCifData(model);
            issueMap = createIssueMapFromCif(model, data.residues, data.groups);
        } else if (params.PDBe_apiSourceJson) {
            const data = await params.PDBe_apiSourceJson(model);
            if (!data) return false;
            info = PropertyWrapper.createInfo();
            issueMap = createIssueMapFromJson(model, data);
        } else {
            return false;
        }

        model.customProperties.add(Descriptor);
        set(model, { info, data: issueMap });
        return true;
    }

    function set(model: Model, prop: Property) {
        (model._dynamicPropertyData.__PDBeStructureQualityReport__ as Property) = prop;
    }

    export function getIssueMap(model: Model): IssueMap | undefined {
        const prop = get(model);
        return prop && prop.data;
    }

    const _emptyArray: string[] = [];
    export function getIssues(e: StructureElement.Location) {
        if (!Unit.isAtomic(e.unit)) return _emptyArray;
        const prop = StructureQualityReport.get(e.unit.model);
        if (!prop || !prop.data) return _emptyArray;
        const rI = e.unit.residueIndex[e.element];
        return prop.data.has(rI) ? prop.data.get(rI)! : _emptyArray;
    }
}

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
        const prop = StructureQualityReport.get(s.model);
        if (prop) info = prop.info;
        if (!prop || !prop.data) continue;

        const { elements, property } = prop.data.getElements(s);
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
    }
}

function createIssueMapFromJson(modelData: Model, data: any): StructureQualityReport.IssueMap | undefined {
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

    return IndexedCustomProperty.fromResidueMap(ret);
}

function createIssueMapFromCif(modelData: Model,
    residueData: Table<typeof StructureQualityReport.Schema.pdbe_structure_quality_report_issues>,
    groupData: Table<typeof StructureQualityReport.Schema.pdbe_structure_quality_report_issue_types>): StructureQualityReport.IssueMap | undefined {

    const ret = new Map<ResidueIndex, string[]>();
    const { label_entity_id, label_asym_id, auth_seq_id, pdbx_PDB_ins_code, issue_type_group_id, pdbx_PDB_model_num, _rowCount } = residueData;

    const groups = parseIssueTypes(groupData);

    for (let i = 0; i < _rowCount; i++) {
        if (pdbx_PDB_model_num.value(i) !== modelData.modelNum) continue;
        const idx = modelData.atomicHierarchy.index.findResidue(label_entity_id.value(i), label_asym_id.value(i), auth_seq_id.value(i), pdbx_PDB_ins_code.value(i));
        ret.set(idx, groups.get(issue_type_group_id.value(i))!);
    }

    return IndexedCustomProperty.fromResidueMap(ret);
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