/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Unit } from '../../../mol-model/structure';
import { CustomProperty } from '../../../mol-model-props/common/custom-property';
import { CustomModelProperty } from '../../../mol-model-props/common/custom-model-property';
import { Model, ResidueIndex } from '../../../mol-model/structure/model';
import { QuerySymbolRuntime } from '../../../mol-script/runtime/query/compiler';
import { CustomPropSymbol } from '../../../mol-script/language/symbol';
import { Type } from '../../../mol-script/language/type';
import { CustomPropertyDescriptor } from '../../../mol-model/custom-property';
import { MmcifFormat } from '../../../mol-model-formats/structure/mmcif';
import { AtomicIndex } from '../../../mol-model/structure/model/properties/atomic';

export { QualityAssessment };

interface LocalPairwiseMetricInfo {
    minResidueIndex: ResidueIndex
    maxResidueIndex: ResidueIndex
    minMetric: number
    maxMetric: number
}

interface QualityAssessment {
    localMetrics: Map<string, Map<ResidueIndex, number>>
    // TODO: use arrays instead of maps?
    localPairwiseMetrics: Map<string, Map<ResidueIndex, Map<ResidueIndex, number>>>
    localPairwiseMetricInfo: Map<string, LocalPairwiseMetricInfo>
    pLDDT?: Map<ResidueIndex, number>
    qmean?: Map<ResidueIndex, number>
}

namespace QualityAssessment {
    const Empty = {
        value: {
            localMetrics: new Map(),
            localPairwiseMetrics: new Map(),
            localPairwiseMetricInfo: new Map(),
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
            for (let i = 0, il = db.ma_qa_metric._rowCount; i < il; i++) {
                if (db.ma_qa_metric.mode.value(i) !== 'local') continue;
                if (localMetricName === db.ma_qa_metric.name.value(i)) return true;
            }
            return false;
        } else {
            return hasLocalMetric;
        }
    }

    export async function obtain(ctx: CustomProperty.Context, model: Model, props: QualityAssessmentProps): Promise<CustomProperty.Data<QualityAssessment>> {
        if (!model || !MmcifFormat.is(model.sourceData)) return Empty;
        const { ma_qa_metric, ma_qa_metric_local, ma_qa_metric_local_pairwise } = model.sourceData.data.db;
        const { model_id, label_asym_id, label_seq_id, metric_id, metric_value } = ma_qa_metric_local;
        const { model_id: model_id_pairwise, label_asym_id_1, label_seq_id_1, label_asym_id_2, label_seq_id_2, metric_id: metric_id_pairwise, metric_value: metric_value_pairwise } = ma_qa_metric_local_pairwise;
        const { index } = model.atomicHierarchy;

        // for simplicity we assume names in ma_qa_metric for mode 'local' are unique
        const localMetrics = new Map<string, Map<ResidueIndex, number>>();
        const localNames = new Map<number, string>();

        const localPairwiseMetrics = new Map<string, Map<ResidueIndex, Map<ResidueIndex, number>>>();
        const localPairwiseMetricInfo = new Map<string, LocalPairwiseMetricInfo>();
        const pairwiseNames = new Map<number, string>();

        for (let i = 0, il = ma_qa_metric._rowCount; i < il; i++) {
            if (ma_qa_metric.mode.value(i) === 'local') {
                const name = ma_qa_metric.name.value(i);
                if (localMetrics.has(name)) {
                    console.warn(`local ma_qa_metric with name '${name}' already added`);
                    continue;
                }

                localMetrics.set(name, new Map());
                localNames.set(ma_qa_metric.id.value(i), name);
            } else if (ma_qa_metric.mode.value(i) === 'local-pairwise') {
                const name = ma_qa_metric.name.value(i);
                if (localPairwiseMetrics.has(name)) {
                    console.warn(`local pairwise ma_qa_metric with name '${name}' already added`);
                    continue;
                }

                localPairwiseMetrics.set(name, new Map());
                pairwiseNames.set(ma_qa_metric.id.value(i), name);
                localPairwiseMetricInfo.set(name, {
                    minResidueIndex: Number.MAX_SAFE_INTEGER as ResidueIndex,
                    maxResidueIndex: Number.MIN_SAFE_INTEGER as ResidueIndex,
                    minMetric: Number.MAX_VALUE,
                    maxMetric: -Number.MAX_VALUE,
                });
            }
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
                const name = localNames.get(metric_id.value(i))!;
                localMetrics.get(name)!.set(rI, metric_value.value(i));
            }
        }

        for (let i = 0, il = ma_qa_metric_local_pairwise._rowCount; i < il; i++) {
            if (model_id_pairwise.value(i) !== model.modelNum) continue;

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

            const name = pairwiseNames.get(metric_id_pairwise.value(i))!;
            const metrics = localPairwiseMetrics.get(name)!;

            let r1 = metrics.get(rI_1);
            if (!r1) {
                r1 = new Map();
                metrics.set(rI_1, r1);
            }
            const value = metric_value_pairwise.value(i);
            r1.set(rI_2, value);

            const info = localPairwiseMetricInfo.get(name)!;
            if (rI_1 > info.maxResidueIndex) info.maxResidueIndex = rI_1;
            if (rI_1 < info.minResidueIndex) info.minResidueIndex = rI_1;
            if (rI_2 > info.maxResidueIndex) info.maxResidueIndex = rI_2;
            if (rI_2 < info.minResidueIndex) info.minResidueIndex = rI_2;
            if (value > info.maxMetric) info.maxMetric = value;
            if (value < info.minMetric) info.minMetric = value;
        }

        return {
            value: {
                localMetrics,
                localPairwiseMetrics,
                localPairwiseMetricInfo,
                pLDDT: localMetrics.get('pLDDT'),
                qmean: localMetrics.get('qmean'),
            }
        };
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