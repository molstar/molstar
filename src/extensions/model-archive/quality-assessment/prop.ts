/**
 * Copyright (c) 2021-24 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifFrame } from '../../../mol-io/reader/cif';
import { toDatabase } from '../../../mol-io/reader/cif/schema';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { MmcifFormat } from '../../../mol-model-formats/structure/mmcif';
import { CustomModelProperty } from '../../../mol-model-props/common/custom-model-property';
import { CustomProperty } from '../../../mol-model-props/common/custom-property';
import { CustomPropertyDescriptor } from '../../../mol-model/custom-property';
import { Unit } from '../../../mol-model/structure';
import { Model, ResidueIndex } from '../../../mol-model/structure/model';
import { AtomicIndex } from '../../../mol-model/structure/model/properties/atomic';
import { CustomPropSymbol } from '../../../mol-script/language/symbol';
import { Type } from '../../../mol-script/language/type';
import { QuerySymbolRuntime } from '../../../mol-script/runtime/query/compiler';
import { ParamDefinition, ParamDefinition as PD } from '../../../mol-util/param-definition';

export { QualityAssessment };

interface QualityAssessment {
    local: QualityAssessment.Local[]
    /** id -> metric info */
    localMap: Map<number, QualityAssessment.Local>

    /** default pLDDT metric */
    pLDDT?: Map<ResidueIndex, number>
    /** default qmean metric */
    qmean?: Map<ResidueIndex, number>

    /**
     * @deprecated
     * NOTE: Keeping this around in case someone is using it
     * TODO: Remove in Mol* 5.0
     */
    localMetrics: Map<string, Map<ResidueIndex, number>>
}

namespace QualityAssessment {
    export interface Local {
        id: number
        kind?: 'pLDDT' | 'qmean';
        type: mmCIF_Schema['ma_qa_metric']['type']['T']
        name: string
        domain?: [number, number]
        valueRange: [number, number]
        values: Map<ResidueIndex, number>
    }

    export interface Pairwise {
        id: number
        name: string

        residueRange: [ResidueIndex, ResidueIndex]
        valueRange: [number, number]
        values: Record<ResidueIndex, Record<ResidueIndex, number | undefined> | undefined>
    }

    const Empty = {
        value: {
            local: [],
            localMap: new Map<number, Local>(),
            localMetrics: new Map(),
        } satisfies QualityAssessment
    };

    export function isApplicable(model?: Model, localMetricName?: 'pLDDT' | 'qmean'): boolean {
        if (!model || !MmcifFormat.is(model.sourceData)) return false;
        const { db } = model.sourceData.data;
        const hasLocalMetric = (
            db.ma_qa_metric.id.isDefined &&
            db.ma_qa_metric_local.ordinal_id.isDefined
        );
        if (localMetricName && hasLocalMetric) {
            const nameQuery = localMetricName.toLowerCase();
            for (let i = 0, il = db.ma_qa_metric._rowCount; i < il; i++) {
                if (db.ma_qa_metric.mode.value(i) !== 'local') continue;
                const name = db.ma_qa_metric.name.value(i).toLowerCase();
                if (name.includes(nameQuery)) return true;
            }
            return false;
        } else {
            return hasLocalMetric;
        }
    }

    export function getLocalOptions(model: Model | undefined, kind: 'pLDDT' | 'qmean') {
        if (!model) return ParamDefinition.Select(undefined, [], { label: 'Metric', isHidden: true });
        const local = QualityAssessmentProvider.get(model).value?.local;
        if (!local) return ParamDefinition.Select(undefined, [], { label: 'Metric', isHidden: true });
        const options = local.filter(m => m.kind === kind).map(m => [m.id, `${m.name} (${m.type})`] as [number, string]);
        return ParamDefinition.Select(options[0]?.[0], options, { label: 'Metric ' });
    }

    export async function obtain(ctx: CustomProperty.Context, model: Model, props: QualityAssessmentProps): Promise<CustomProperty.Data<QualityAssessment>> {
        if (!model || !MmcifFormat.is(model.sourceData)) return Empty;
        const { ma_qa_metric, ma_qa_metric_local } = model.sourceData.data.db;
        const { model_id, label_asym_id, label_seq_id, metric_id, metric_value } = ma_qa_metric_local;
        const { index } = model.atomicHierarchy;

        const locals: Local[] = [];
        const localMetrics = new Map<number, Local>();

        // for simplicity we assume names in ma_qa_metric for mode 'local' are unique
        const localMetricValues = new Map<string, Map<ResidueIndex, number>>();
        const localNames = new Map<number, string>();

        for (let i = 0, il = ma_qa_metric._rowCount; i < il; i++) {
            if (ma_qa_metric.mode.value(i) !== 'local') continue;

            const id = ma_qa_metric.id.value(i);
            const name = ma_qa_metric.name.value(i);
            const type = ma_qa_metric.type.value(i);
            const values: Map<ResidueIndex, number> = new Map();
            const ispPLDDType = type.toLowerCase().includes('plddt');
            const has01Range = type.replace(/\s/g, '').includes('[0,1]');

            let domain: [number, number] | undefined;
            if (has01Range) {
                domain = [0, 1];
            } else if (ispPLDDType) {
                domain = [0, 100];
            }

            let kind: 'pLDDT' | 'qmean' | undefined;

            const nameLower = name.toLowerCase();
            if (nameLower.includes('plddt')) kind = 'pLDDT';
            else if (nameLower.includes('qmean')) kind = 'qmean';

            if (!kind && ispPLDDType) kind = 'pLDDT';
            if (!kind && has01Range) kind = 'qmean';

            const metric: Local = { id, kind, type, name, domain, valueRange: [Number.MAX_VALUE, -Number.MAX_VALUE], values };

            localMetrics.set(id, metric);
            locals.push(metric);

            if (localMetricValues.has(name)) {
                console.warn(`local ma_qa_metric with name '${name}' already added`);
                continue;
            }

            localMetricValues.set(name, values);
            localNames.set(id, name);
        }

        const residueKey: AtomicIndex.ResidueLabelKey = {
            label_entity_id: '',
            label_asym_id: '',
            label_seq_id: 0,
            pdbx_PDB_ins_code: undefined,
        };

        for (let i = 0, il = ma_qa_metric_local._rowCount; i < il; i++) {
            if (model_id.value(i) !== model.modelNum) continue;

            const labelAsymId = label_asym_id.value(i);
            const entityIndex = index.findEntity(labelAsymId);

            residueKey.label_entity_id = model.entities.data.id.value(entityIndex);
            residueKey.label_asym_id = labelAsymId;
            residueKey.label_seq_id = label_seq_id.value(i);

            const rI = index.findResidueLabel(residueKey);
            if (rI >= 0) {
                const entry = localMetrics.get(metric_id.value(i));
                if (!entry) continue;

                const value = metric_value.value(i);
                const range = entry.valueRange;
                if (value < range[0]) range[0] = value;
                if (value > range[1]) range[1] = value;
                entry.values.set(rI, value);
            }
        }

        return {
            value: {
                local: Array.from(localMetrics.values()),
                localMap: localMetrics,
                pLDDT: locals.find(m => m.kind === 'pLDDT')?.values,
                qmean: locals.find(m => m.kind === 'qmean')?.values,

                localMetrics: localMetricValues,
            }
        };
    }

    const PairwiseSchema = {
        ma_qa_metric: mmCIF_Schema.ma_qa_metric,
        ma_qa_metric_local_pairwise: mmCIF_Schema.ma_qa_metric_local_pairwise
    };

    export function findModelArchiveCIFPAEMetrics(frame: CifFrame) {
        const { ma_qa_metric, ma_qa_metric_local_pairwise } = toDatabase(PairwiseSchema, frame);
        const result: { id: number, name: string }[] = [];
        if (ma_qa_metric_local_pairwise._rowCount === 0) return result;

        for (let i = 0, il = ma_qa_metric._rowCount; i < il; i++) {
            if (ma_qa_metric.mode.value(i) !== 'local-pairwise') continue;
            const id = ma_qa_metric.id.value(i);
            const name = ma_qa_metric.name.value(i);
            if (!name.toLowerCase().includes('pae')) continue;
            result.push({ id, name });
        }
        return result;
    }

    export function pairwiseMetricFromModelArchiveCIF(model: Model, frame: CifFrame, metricId: number): Pairwise | undefined {
        const db = toDatabase(PairwiseSchema, frame);
        if (!db.ma_qa_metric_local_pairwise._rowCount) return undefined;

        const { ma_qa_metric, ma_qa_metric_local_pairwise } = db;
        const { model_id, label_asym_id_1, label_seq_id_1, label_asym_id_2, label_seq_id_2, metric_id, metric_value } = db.ma_qa_metric_local_pairwise;
        const { index } = model.atomicHierarchy;

        let metric: Pairwise | undefined;

        for (let i = 0, il = ma_qa_metric._rowCount; i < il; i++) {
            if (ma_qa_metric.mode.value(i) !== 'local-pairwise') continue;
            const id = ma_qa_metric.id.value(i);
            if (id !== metricId) continue;

            const name = ma_qa_metric.name.value(i);
            metric = {
                id,
                name,
                residueRange: [Number.MAX_SAFE_INTEGER as ResidueIndex, Number.MIN_SAFE_INTEGER as ResidueIndex],
                valueRange: [Number.MAX_VALUE, -Number.MAX_VALUE],
                values: {}
            };
        }

        if (!metric) return undefined;

        const { values, residueRange, valueRange } = metric;
        const residueKey: AtomicIndex.ResidueLabelKey = {
            label_entity_id: '',
            label_asym_id: '',
            label_seq_id: 0,
            pdbx_PDB_ins_code: undefined,
        };

        for (let i = 0, il = ma_qa_metric_local_pairwise._rowCount; i < il; i++) {
            if (model_id.value(i) !== model.modelNum || metric_id.value(i) !== metricId) continue;

            let labelAsymId = label_asym_id_1.value(i);
            let entityIndex = index.findEntity(labelAsymId);
            residueKey.label_entity_id = model.entities.data.id.value(entityIndex);
            residueKey.label_asym_id = labelAsymId;
            residueKey.label_seq_id = label_seq_id_1.value(i);

            const rI_1 = index.findResidueLabel(residueKey);
            if (rI_1 < 0) continue;

            labelAsymId = label_asym_id_2.value(i);
            entityIndex = index.findEntity(labelAsymId);
            residueKey.label_entity_id = model.entities.data.id.value(entityIndex);
            residueKey.label_asym_id = labelAsymId;
            residueKey.label_seq_id = label_seq_id_2.value(i);

            const rI_2 = index.findResidueLabel(residueKey);
            if (rI_1 < 0) continue;

            let r1 = values[rI_1];
            if (!r1) {
                r1 = {};
                values[rI_1] = r1;
            }
            const value = metric_value.value(i);
            r1[rI_2] = value;

            if (rI_1 < residueRange[0]) residueRange[0] = rI_1;
            if (rI_2 < residueRange[0]) residueRange[0] = rI_2;
            if (rI_1 > residueRange[1]) residueRange[1] = rI_1;
            if (rI_2 > residueRange[1]) residueRange[1] = rI_2;
            if (value < valueRange[0]) valueRange[0] = value;
            if (value > valueRange[1]) valueRange[1] = value;
        }

        return metric;
    }

    export const symbols = {
        pLDDT: QuerySymbolRuntime.Dynamic(CustomPropSymbol('ma', 'quality-assessment.pLDDT', Type.Num),
            ctx => {
                const { unit, element } = ctx.element;
                if (!Unit.isAtomic(unit)) return -1;
                const qualityAssessment = QualityAssessmentProvider.get(unit.model).value;
                return qualityAssessment?.pLDDT?.get(unit.model.atomicHierarchy.residueAtomSegments.index[element]) ?? -1;
            }
        ),
        qmean: QuerySymbolRuntime.Dynamic(CustomPropSymbol('ma', 'quality-assessment.qmean', Type.Num),
            ctx => {
                const { unit, element } = ctx.element;
                if (!Unit.isAtomic(unit)) return -1;
                const qualityAssessment = QualityAssessmentProvider.get(unit.model).value;
                return qualityAssessment?.qmean?.get(unit.model.atomicHierarchy.residueAtomSegments.index[element]) ?? -1;
            }
        ),
    };
}

export const QualityAssessmentParams = { };
export type QualityAssessmentParams = typeof QualityAssessmentParams
export type QualityAssessmentProps = PD.Values<QualityAssessmentParams>

export const QualityAssessmentProvider: CustomModelProperty.Provider<QualityAssessmentParams, QualityAssessment> = CustomModelProperty.createProvider({
    label: 'QualityAssessment',
    descriptor: CustomPropertyDescriptor({
        name: 'ma_quality_assessment',
        symbols: QualityAssessment.symbols
    }),
    type: 'static',
    defaultParams: QualityAssessmentParams,
    getParams: (data: Model) => QualityAssessmentParams,
    isApplicable: (data: Model) => QualityAssessment.isApplicable(data),
    obtain: async (ctx: CustomProperty.Context, data: Model, props: Partial<QualityAssessmentProps>) => {
        const p = { ...PD.getDefaultValues(QualityAssessmentParams), ...props };
        return await QualityAssessment.obtain(ctx, data, p);
    }
});