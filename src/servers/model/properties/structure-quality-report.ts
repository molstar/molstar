/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ResidueIndex, ModelPropertyDescriptor, Model, Structure, Unit, StructureElement, StructureProperties as P  } from 'mol-model/structure';
import fetch from 'node-fetch';
import { CifWriter } from 'mol-io/writer/cif';
import CifField = CifWriter.Field;
import { Segmentation } from 'mol-data/int';

type IssueMap = Map<ResidueIndex, string[]>

const _Descriptor: ModelPropertyDescriptor = {
    isStatic: true,
    name: 'structure_quality_report',
    cifExport: {
        categories: [{
            name: 'structure_quality_report',
            instance(ctx) {
                const issues = StructureQualityReport.get(ctx.model);
                if (!issues) return CifWriter.Category.Empty;

                const residues = getResidueLoci(ctx.structure, issues);
                return {
                    fields: _structure_quality_report_fields,
                    data: <ExportCtx>{ model: ctx.model, residues, residueIndex: ctx.model.atomicHierarchy.residueAtomSegments.index, issues },
                    rowCount: residues.length
                };
            }
        }]
    }
}

type ExportCtx = { model: Model, residueIndex: ArrayLike<ResidueIndex>, residues: StructureElement[], issues: IssueMap };

const _structure_quality_report_fields: CifField<ResidueIndex, ExportCtx>[] = [
    CifField.str<ResidueIndex, ExportCtx>('label_comp_id', (i, d) => P.residue.label_comp_id(d.residues[i])),
    CifField.int<ResidueIndex, ExportCtx>('label_seq_id', (i, d) => P.residue.label_seq_id(d.residues[i])),
    CifField.str<ResidueIndex, ExportCtx>('pdbx_PDB_ins_code', (i, d) => P.residue.pdbx_PDB_ins_code(d.residues[i])),
    CifField.str<ResidueIndex, ExportCtx>('label_asym_id', (i, d) => P.chain.label_asym_id(d.residues[i])),
    CifField.str<ResidueIndex, ExportCtx>('label_entity_id', (i, d) => P.entity.id(d.residues[i])),

    CifField.str<ResidueIndex, ExportCtx>('auth_comp_id', (i, d) => P.residue.auth_comp_id(d.residues[i])),
    CifField.int<ResidueIndex, ExportCtx>('auth_seq_id', (i, d) => P.residue.auth_seq_id(d.residues[i])),
    CifField.str<ResidueIndex, ExportCtx>('auth_asym_id', (i, d) => P.chain.auth_asym_id(d.residues[i])),

    CifField.str<ResidueIndex, ExportCtx>('issues', (i, d) => d.issues.get(d.residueIndex[d.residues[i].element])!.join(','))
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

function createIssueMap(modelData: Model, data: any): IssueMap | undefined {
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

export namespace StructureQualityReport {
    export const Descriptor = _Descriptor;

    export async function attachFromPDBeApi(model: Model) {
        if (model.customProperties.has(Descriptor)) return true;

        const id = model.label.toLowerCase();
        const rawData = await fetch(`https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry/${model.label.toLowerCase()}`, { timeout: 1500 });
        const json = await rawData.json();
        const data = json[id];
        if (!data) return false;
        const issueMap = createIssueMap(model, data);
        if (!issueMap || issueMap.size === 0) return false;

        model.customProperties.add(Descriptor);
        model._staticPropertyData.__StructureQualityReport__ = issueMap;
        return true;
    }

    export function get(model: Model): IssueMap | undefined {
        return model._staticPropertyData.__StructureQualityReport__;
    }
}