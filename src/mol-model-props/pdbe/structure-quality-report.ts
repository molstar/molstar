/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Segmentation } from 'mol-data/int';
import { CifWriter } from 'mol-io/writer/cif';
import { Model, ModelPropertyDescriptor, ResidueIndex, Structure, StructureElement, Unit } from 'mol-model/structure';
import { residueIdFields } from 'mol-model/structure/export/categories/atom_site';
import CifField = CifWriter.Field;
import { mmCIF_residueId_schema } from 'mol-io/reader/cif/schema/mmcif-extras';
import { Column, Table } from 'mol-data/db';
import { toTable } from 'mol-io/reader/cif/schema';

type IssueMap = Map<ResidueIndex, string[]>

const _Descriptor: ModelPropertyDescriptor = {
    isStatic: true,
    name: 'structure_quality_report',
    cifExport: {
        prefix: 'pdbe',
        categories: [{
            name: 'pdbe_structure_quality_report',
            instance() {
                return { fields: _structure_quality_report_fields, rowCount: 1 }
            }
        }, {
            name: 'pdbe_structure_quality_report_issues',
            instance(ctx) {
                const issues = StructureQualityReport.get(ctx.model);
                if (!issues) return CifWriter.Category.Empty;

                const residues = getResidueLoci(ctx.structure, issues);
                return {
                    fields: _structure_quality_report_issues_fields,
                    data: <ExportCtx>{ model: ctx.model, residues, residueIndex: ctx.model.atomicHierarchy.residueAtomSegments.index, issues },
                    rowCount: residues.length
                };
            }
        }]
    }
}

type ExportCtx = { model: Model, residueIndex: ArrayLike<ResidueIndex>, residues: StructureElement[], issues: IssueMap };

const _structure_quality_report_issues_fields: CifField<ResidueIndex, ExportCtx>[] = [
    CifField.index('id'),
    ...residueIdFields<ResidueIndex, ExportCtx>((k, d) => d.residues[k]),
    CifField.str<ResidueIndex, ExportCtx>('issues', (i, d) => d.issues.get(d.residueIndex[d.residues[i].element])!.join(','))
];

const _structure_quality_report_fields: CifField<ResidueIndex, ExportCtx>[] = [
    CifField.str('updated_datetime_utc', () => `${new Date().toISOString().replace(/T/, ' ').replace(/\..+/, '')}`)
];

function getResidueLoci(structure: Structure, issues: IssueMap) {
    const seenResidues = new Set<ResidueIndex>();
    const unitGroups = structure.unitSymmetryGroups;
    const loci: StructureElement[] = [];

    for (const unitGroup of unitGroups) {
        const unit = unitGroup.units[0];
        if (!Unit.isAtomic(unit)) {
            continue;
        }

        const residues = Segmentation.transientSegments(unit.model.atomicHierarchy.residueAtomSegments, unit.elements);
        while (residues.hasNext) {
            const seg = residues.move();
            if (!issues.has(seg.index) || seenResidues.has(seg.index)) continue;

            seenResidues.add(seg.index);
            loci[loci.length] = StructureElement.create(unit, unit.elements[seg.start]);
        }
    }

    loci.sort((x, y) => x.element - y.element);
    return loci;
}

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
                    const idx = modelData.atomicHierarchy.findResidueKey(entity_id, asym_id, '', auth_seq_id, ins_code);
                    ret.set(idx, residue.outlier_types);
                }
            }
        }
    }

    return ret;
}

function createIssueMapFromCif(modelData: Model, data: Table<typeof StructureQualityReport.Schema.pdbe_structure_quality_report_issues>): IssueMap | undefined {
    const ret = new Map<ResidueIndex, string[]>();
    const { label_entity_id, label_asym_id, auth_seq_id, pdbx_PDB_ins_code, issues, _rowCount } = data;

    for (let i = 0; i < _rowCount; i++) {
        const idx = modelData.atomicHierarchy.findResidueKey(label_entity_id.value(i), label_asym_id.value(i), '', auth_seq_id.value(i), pdbx_PDB_ins_code.value(i));
        ret.set(idx, issues.value(i));
    }

    return ret;
}

export namespace StructureQualityReport {
    export const Descriptor = _Descriptor;

    export const Schema = {
        pdbe_structure_quality_report: {
            updated_datetime_utc: Column.Schema.str
        },
        pdbe_structure_quality_report_issues: {
            id: Column.Schema.int,
            ...mmCIF_residueId_schema,
            issues: Column.Schema.List(',', x => x)
        }
    }

    export async function attachFromCifOrApi(model: Model, params: {
        // provide JSON from api
        PDBe_apiSourceJson?: (model: Model) => Promise<any>
    }) {
        if (model.customProperties.has(Descriptor)) return true;

        let issueMap;

        if (model.sourceData.kind === 'mmCIF' && model.sourceData.frame.categoryNames.includes('pdbe_structure_quality_report')) {
            const data = toTable(Schema.pdbe_structure_quality_report_issues, model.sourceData.frame.categories.pdbe_structure_quality_report);
            issueMap = createIssueMapFromCif(model, data);
        } else if (params.PDBe_apiSourceJson) {
            const id = model.label.toLowerCase();
            const json = await params.PDBe_apiSourceJson(model);
            const data = json[id];
            if (!data) return false;
            issueMap = createIssueMapFromJson(model, data);
        } else {
            return false;
        }

        model.customProperties.add(Descriptor);
        model._staticPropertyData.__StructureQualityReport__ = issueMap;
        return true;
    }

    export function get(model: Model): IssueMap | undefined {
        return model._staticPropertyData.__StructureQualityReport__;
    }
}