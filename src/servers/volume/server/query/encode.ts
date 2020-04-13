/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CifWriter } from '../../../../mol-io/writer/cif';
import * as Data from './data-model';
import * as Coords from '../algebra/coordinate';
import { VOLUME_SERVER_VERSION as VERSION } from '../version';
import * as DataFormat from '../../common/data-format';
import { Column } from '../../../../mol-data/db';
import { ArrayEncoding, ArrayEncoder } from '../../../../mol-io/common/binary-cif';
import { TypedArrayValueType, TypedArrayValueArray } from '../../../../mol-io/common/typed-array';

export default function encode(query: Data.QueryContext, output: Data.QueryOutputStream) {
    let w = CifWriter.createEncoder({ binary: query.params.asBinary, encoderName: `VolumeServer ${VERSION}` });
    write(w, query);
    w.encode();
    w.writeTo(output);
}

interface ResultContext {
    query: Data.QueryContext.Data,
    channelIndex: number
}

function string<T>(name: string, str: (data: T) => string, isSpecified?: (data: T) => boolean): CifWriter.Field<number, T> {
    if (isSpecified) {
        return CifWriter.Field.str(name, (i, d) => str(d), { valueKind: (i, d) => isSpecified(d) ? Column.ValueKind.Present : Column.ValueKind.NotPresent });
    }
    return CifWriter.Field.str(name, (i, d) => str(d));
}

function int32<T>(name: string, value: (data: T) => number): CifWriter.Field<number, T> {
    return CifWriter.Field.int(name, (i, d) => value(d));
}

function float64<T>(name: string, value: (data: T) => number, digitCount: number = 6): CifWriter.Field<number, T> {
    return CifWriter.Field.float(name, (i, d) => value(d), { digitCount: digitCount, typedArray: Float64Array });
}

interface _vd3d_Ctx {
    header: DataFormat.Header,
    channelIndex: number,
    grid: Coords.GridDomain<'Query'>,
    sampleRate: number,
    globalValuesInfo: DataFormat.ValuesInfo,
    sampledValuesInfo: DataFormat.ValuesInfo,
}

const _volume_data_3d_info_fields = [
    string<_vd3d_Ctx>('name', ctx => ctx.header.channels[ctx.channelIndex]),

    int32<_vd3d_Ctx>('axis_order[0]', ctx => ctx.header.axisOrder[0]),
    int32<_vd3d_Ctx>('axis_order[1]', ctx => ctx.header.axisOrder[1]),
    int32<_vd3d_Ctx>('axis_order[2]', ctx => ctx.header.axisOrder[2]),

    float64<_vd3d_Ctx>('origin[0]', ctx => ctx.grid.origin[0]),
    float64<_vd3d_Ctx>('origin[1]', ctx => ctx.grid.origin[1]),
    float64<_vd3d_Ctx>('origin[2]', ctx => ctx.grid.origin[2]),

    float64<_vd3d_Ctx>('dimensions[0]', ctx => ctx.grid.dimensions[0]),
    float64<_vd3d_Ctx>('dimensions[1]', ctx => ctx.grid.dimensions[1]),
    float64<_vd3d_Ctx>('dimensions[2]', ctx => ctx.grid.dimensions[2]),

    int32<_vd3d_Ctx>('sample_rate', ctx => ctx.sampleRate),
    int32<_vd3d_Ctx>('sample_count[0]', ctx => ctx.grid.sampleCount[0]),
    int32<_vd3d_Ctx>('sample_count[1]', ctx => ctx.grid.sampleCount[1]),
    int32<_vd3d_Ctx>('sample_count[2]', ctx => ctx.grid.sampleCount[2]),

    int32<_vd3d_Ctx>('spacegroup_number', ctx => ctx.header.spacegroup.number),

    float64<_vd3d_Ctx>('spacegroup_cell_size[0]', ctx => ctx.header.spacegroup.size[0], 3),
    float64<_vd3d_Ctx>('spacegroup_cell_size[1]', ctx => ctx.header.spacegroup.size[1], 3),
    float64<_vd3d_Ctx>('spacegroup_cell_size[2]', ctx => ctx.header.spacegroup.size[2], 3),

    float64<_vd3d_Ctx>('spacegroup_cell_angles[0]', ctx => ctx.header.spacegroup.angles[0], 3),
    float64<_vd3d_Ctx>('spacegroup_cell_angles[1]', ctx => ctx.header.spacegroup.angles[1], 3),
    float64<_vd3d_Ctx>('spacegroup_cell_angles[2]', ctx => ctx.header.spacegroup.angles[2], 3),

    float64<_vd3d_Ctx>('mean_source', ctx => ctx.globalValuesInfo.mean),
    float64<_vd3d_Ctx>('mean_sampled', ctx => ctx.sampledValuesInfo.mean),
    float64<_vd3d_Ctx>('sigma_source', ctx => ctx.globalValuesInfo.sigma),
    float64<_vd3d_Ctx>('sigma_sampled', ctx => ctx.sampledValuesInfo.sigma),
    float64<_vd3d_Ctx>('min_source', ctx => ctx.globalValuesInfo.min),
    float64<_vd3d_Ctx>('min_sampled', ctx => ctx.sampledValuesInfo.min),
    float64<_vd3d_Ctx>('max_source', ctx => ctx.globalValuesInfo.max),
    float64<_vd3d_Ctx>('max_sampled', ctx => ctx.sampledValuesInfo.max)
];

const _volume_data_3d_info: CifWriter.Category<ResultContext> = {
    name: 'volume_data_3d_info',
    instance(result) {
        const ctx: _vd3d_Ctx = {
            header: result.query.data.header,
            channelIndex: result.channelIndex,
            grid: result.query.samplingInfo.gridDomain,
            sampleRate: result.query.samplingInfo.sampling.rate,
            globalValuesInfo: result.query.data.header.sampling[0].valuesInfo[result.channelIndex],
            sampledValuesInfo: result.query.data.header.sampling[result.query.samplingInfo.sampling.index].valuesInfo[result.channelIndex]
        };

        return { fields: _volume_data_3d_info_fields, source: [{ data: ctx, rowCount: 1 }] };
    }
};

function _volume_data_3d_number(i: number, ctx: TypedArrayValueArray): number {
    return ctx[i];
}

const _volume_data_3d: CifWriter.Category<ResultContext> = {
    name: 'volume_data_3d',
    instance(ctx) {
        const data = ctx.query.values[ctx.channelIndex];

        const E = ArrayEncoding;
        let encoder: ArrayEncoder;
        let typedArray: any;
        if (ctx.query.data.header.valueType === TypedArrayValueType.Float32 || ctx.query.data.header.valueType === TypedArrayValueType.Int16) {
            let min: number, max: number;
            min = data[0], max = data[0];
            for (let i = 0, n = data.length; i < n; i++) {
                let v = data[i];
                if (v < min) min = v;
                else if (v > max) max = v;
            }
            typedArray = Float32Array;
            // encode into 255 steps and store each value in 1 byte.
            encoder = E.by(E.intervalQuantizaiton(min, max, 255, Uint8Array)).and(E.byteArray);
        } else {
            typedArray = Int8Array;
            // just encode the bytes
            encoder = E.by(E.byteArray);
        }

        const fields = [CifWriter.Field.float('values', _volume_data_3d_number, { encoder, typedArray, digitCount: 6 })];
        return CifWriter.categoryInstance(fields, { data, rowCount: data.length });
    }
};

function pickQueryBoxDimension(ctx: Data.QueryContext, e: 'a' | 'b', d: number) {
    const box = ctx.params.box;
    switch (box.kind) {
        case 'Cartesian':
        case 'Fractional':
            return `${Math.round(1000000 * box[e][d]) / 1000000}`;
        default: return '';
    }
}

function queryBoxDimension(e: 'a' | 'b', d: number) {
    return string<Data.QueryContext>(`query_box_${e}[${d}]`, ctx => pickQueryBoxDimension(ctx, e, d), ctx => ctx.params.box.kind !== 'Cell');
}

const _density_server_result_fields = [
    string<Data.QueryContext>('server_version', ctx => VERSION),
    string<Data.QueryContext>('datetime_utc', ctx => new Date().toISOString().replace(/T/, ' ').replace(/\..+/, '')),
    string<Data.QueryContext>('guid', ctx => ctx.guid),
    string<Data.QueryContext>('is_empty', ctx => ctx.kind === 'Empty' || ctx.kind === 'Error' ? 'yes' : 'no'),
    string<Data.QueryContext>('has_error', ctx => ctx.kind === 'Error' ? 'yes' : 'no'),
    string<Data.QueryContext>('error', ctx => ctx.kind === 'Error' ? ctx.message : '', (ctx) => ctx.kind === 'Error'),
    string<Data.QueryContext>('query_source_id', ctx => ctx.params.sourceId),
    string<Data.QueryContext>('query_type', ctx => 'box'),
    string<Data.QueryContext>('query_box_type', ctx => ctx.params.box.kind.toLowerCase()),
    queryBoxDimension('a', 0),
    queryBoxDimension('a', 1),
    queryBoxDimension('a', 2),
    queryBoxDimension('b', 0),
    queryBoxDimension('b', 1),
    queryBoxDimension('b', 2)
];

const _density_server_result: CifWriter.Category<Data.QueryContext> = {
    name: 'density_server_result',
    instance: ctx => CifWriter.categoryInstance(_density_server_result_fields, { data: ctx, rowCount: 1 })
};

function write(encoder: CifWriter.Encoder, query: Data.QueryContext) {
    encoder.startDataBlock('SERVER');
    encoder.writeCategory(_density_server_result, query);

    switch (query.kind) {
        case 'Data':
    }

    if (query.kind === 'Data') {
        const header = query.data.header;
        for (let i = 0; i < header.channels.length; i++) {
            encoder.startDataBlock(header.channels[i]);
            const ctx: ResultContext = { query, channelIndex: i };

            encoder.writeCategory(_volume_data_3d_info, ctx);
            encoder.writeCategory(_volume_data_3d, ctx);
        }
    }
}