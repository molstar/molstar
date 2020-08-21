/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column, Table } from '../../mol-data/db';
import { createModels } from '../../mol-model-formats/structure/basic/parser';
import { BasicSchema, createBasic } from '../../mol-model-formats/structure/basic/schema';
import { EntityBuilder } from '../../mol-model-formats/structure/common/entity';
import { CustomModelProperty } from '../../mol-model-props/common/custom-model-property';
import { Model, Trajectory, Unit } from '../../mol-model/structure';
import { MoleculeType } from '../../mol-model/structure/model/types';
import { CustomPropSymbol } from '../../mol-script/language/symbol';
import Type from '../../mol-script/language/type';
import { QuerySymbolRuntime } from '../../mol-script/runtime/query/base';
import { Task, RuntimeContext } from '../../mol-task';
import { objectForEach } from '../../mol-util/object';
import { ParamDefinition } from '../../mol-util/param-definition';
import { MolScriptBuilder as MS } from '../../mol-script/language/builder';
import { G3dDataBlock } from './data';
import { Loci } from '../../mol-model/loci';
import { LociLabelProvider } from '../../mol-plugin-state/manager/loci-label';
import { OrderedSet } from '../../mol-data/int';

interface NormalizedData {
    entity_id: string[],
    chromosome: string[],
    seq_id_begin: Int32Array,
    seq_id_end: Int32Array,
    start: Int32Array,
    x: Float32Array,
    y: Float32Array,
    z: Float32Array,
    haplotype: string[]
}

function getColumns(block: G3dDataBlock) {
    const { data } = block;
    let size = 0;

    objectForEach(data, h => objectForEach(h, g => size += g.start.length));

    const normalized: NormalizedData = {
        entity_id: new Array(size),
        chromosome: new Array(size),
        seq_id_begin: new Int32Array(size),
        seq_id_end: new Int32Array(size),
        start: new Int32Array(size),
        x: new Float32Array(size),
        y: new Float32Array(size),
        z: new Float32Array(size),
        haplotype: new Array(size)
    };

    let o = 0;
    objectForEach(data, (hs, h) => {
        objectForEach(hs, (chs, ch) => {
            const entity_id = `${ch}-${h}`;
            for (let i = 0, _i = chs.start.length; i < _i; i++) {
                normalized.entity_id[o] = entity_id;
                normalized.chromosome[o] = ch;
                normalized.start[o] = chs.start[i];
                normalized.seq_id_begin[o] = o;
                normalized.seq_id_end[o] = o;
                normalized.x[o] = 10 * chs.x[i];
                normalized.y[o] = 10 * chs.y[i];
                normalized.z[o] = 10 * chs.z[i];
                normalized.haplotype[o] = h;
                o++;
            }
        });
    });

    return normalized;
}

async function getTraj(ctx: RuntimeContext, data: G3dDataBlock) {
    const normalized = getColumns(data);

    const rowCount = normalized.seq_id_begin.length;
    const entityIds = new Array<string>(rowCount);
    const entityBuilder = new EntityBuilder();

    const stride = normalized.seq_id_begin[1] - normalized.seq_id_begin[0];
    const objectRadius = stride / 3500;

    const eName = { customName: '' };
    for (let i = 0; i < rowCount; ++i) {
        const e = normalized.entity_id[i];
        eName.customName = e;
        const entityId = entityBuilder.getEntityId(e, MoleculeType.DNA, e, eName);
        entityIds[i] = entityId;
    }

    const ihm_sphere_obj_site = Table.ofPartialColumns(BasicSchema.ihm_sphere_obj_site, {
        id: Column.range(0, rowCount),
        entity_id: Column.ofStringArray(entityIds),
        seq_id_begin: Column.ofIntArray(normalized.seq_id_begin),
        seq_id_end: Column.ofIntArray(normalized.seq_id_end),
        asym_id: Column.ofStringArray(normalized.chromosome),

        Cartn_x: Column.ofFloatArray(normalized.x),
        Cartn_y: Column.ofFloatArray(normalized.y),
        Cartn_z: Column.ofFloatArray(normalized.z),

        object_radius: Column.ofConst(objectRadius, rowCount, Column.Schema.float),
        rmsf: Column.ofConst(0, rowCount, Column.Schema.float),
        model_id: Column.ofConst(1, rowCount, Column.Schema.int),
    }, rowCount);

    const basic = createBasic({
        entity: entityBuilder.getEntityTable(),
        ihm_model_list: Table.ofPartialColumns(BasicSchema.ihm_model_list, {
            model_id: Column.ofIntArray([1]),
            model_name: Column.ofStringArray(['3DG Model']),
        }, 1),
        ihm_sphere_obj_site
    });

    const models = await createModels(basic, { kind: 'g3d', name: 'G3D', data }, ctx);

    await G3dInfoPropertyProvider.attach({ runtime: ctx, assetManager: void 0 as any }, models.representative, {
        resolution: data.resolution,
        haplotype: normalized.haplotype,
        start: normalized.start
    });

    return models;
}

export function trajectoryFromG3D(data: G3dDataBlock): Task<Trajectory> {
    return Task.create('Parse G3D', async ctx => {
        return getTraj(ctx, data);
    });
}

export const G3dSymbols = {
    haplotype: QuerySymbolRuntime.Dynamic(CustomPropSymbol('g3d', 'haplotype', Type.Str),
        ctx => {
            if (Unit.isAtomic(ctx.element.unit)) return '';
            const info = G3dInfoPropertyProvider.get(ctx.element.unit.model);
            const seqId = ctx.element.unit.model.coarseHierarchy.spheres.seq_id_begin.value(ctx.element.element);
            return info.value?.haplotype[seqId] || '';
        }
    )
};

export function g3dHaplotypeQuery(haplotype: 'maternal' | 'paternal') {
    return MS.struct.generator.atomGroups({
        'chain-test': MS.core.rel.eq([G3dSymbols.haplotype.symbol(), haplotype]),
    });
}

export const G3dInfoPropertyParams = {
    haplotype: ParamDefinition.Value<string[]>([]),
    start: ParamDefinition.Value<Int32Array>(new Int32Array(0)),
    resolution: ParamDefinition.Numeric(0)
};

export type G3dInfoPropertyParams = typeof G3dInfoPropertyParams
export type G3dInfoPropertyParamsProps = ParamDefinition.Values<G3dInfoPropertyParams>

export const G3dInfoPropertyProvider: CustomModelProperty.Provider<G3dInfoPropertyParams, G3dInfoPropertyParamsProps> = CustomModelProperty.createProvider({
    label: 'G3d Info',
    type: 'static',
    isHidden: true,
    defaultParams: G3dInfoPropertyParams,
    descriptor: { name: 'g3d_info', symbols: G3dSymbols },
    getParams: (data: Model) => G3dInfoPropertyParams,
    isApplicable: (data: Model) => true,
    obtain: async (ctx, data, props) => {
        return { value: props as any };
    }
});

export const G3dLabelProvider: LociLabelProvider = {
    label: (e: Loci): string | undefined => {
        if (e.kind !== 'element-loci' || Loci.isEmpty(e)) return;

        const first = e.elements[0];
        if (e.elements.length !== 1 || Unit.isAtomic(first.unit)) return;
        const info = G3dInfoPropertyProvider.get(first.unit.model);
        if (!info) return;

        const eI = first.unit.elements[OrderedSet.getAt(first.indices, 0)];
        const seqId = first.unit.model.coarseHierarchy.spheres.seq_id_begin.value(eI);
        return `<b>Start:</b> ${info.value?.start[seqId]} <small>| resolution ${info.value?.resolution}<small>`;
    }
};