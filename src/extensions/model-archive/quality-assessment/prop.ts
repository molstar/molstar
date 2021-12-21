/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
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

export { QualityAssessment };

interface QualityAssessment {
    localMetrics: Map<string, Map<ResidueIndex, number>>
    pLDDT?: Map<ResidueIndex, number>
    qmean?: Map<ResidueIndex, number>
}

namespace QualityAssessment {
    const Empty = {
        value: {
            localMetrics: new Map()
        }
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
        const { ma_qa_metric, ma_qa_metric_local } = model.sourceData.data.db;
        const { model_id, label_asym_id, label_seq_id, metric_id, metric_value } = ma_qa_metric_local;
        const { index } = model.atomicHierarchy;

        // for simplicity we assume names in ma_qa_metric for mode 'local' are unique
        const localMetrics = new Map<string, Map<ResidueIndex, number>>();
        const localNames = new Map<number, string>();

        for (let i = 0, il = ma_qa_metric._rowCount; i < il; i++) {
            if (ma_qa_metric.mode.value(i) !== 'local') continue;

            const name = ma_qa_metric.name.value(i);
            if (localMetrics.has(name)) {
                console.warn(`local ma_qa_metric with name '${name}' already added`);
                continue;
            }

            localMetrics.set(name, new Map());
            localNames.set(ma_qa_metric.id.value(i), name);
        }

        for (let i = 0, il = ma_qa_metric_local._rowCount; i < il; i++) {
            if (model_id.value(i) !== model.modelNum) continue;

            const labelAsymId = label_asym_id.value(i);
            const entityIndex = index.findEntity(labelAsymId);
            const rI = index.findResidue(model.entities.data.id.value(entityIndex), labelAsymId, label_seq_id.value(i));
            const name = localNames.get(metric_id.value(i))!;
            localMetrics.get(name)!.set(rI, metric_value.value(i));
        }

        return {
            value: {
                localMetrics,
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