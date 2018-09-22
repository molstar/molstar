/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util'
import { Mat4 } from 'mol-math/linear-algebra'
import { transformPositionArray/* , transformDirectionArray, getNormalMatrix */ } from '../../util';
import { Geometry } from '../geometry';
import { RuntimeContext } from 'mol-task';
import { createColors } from '../color-data';
import { createMarkers } from '../marker-data';
import { createSizes } from '../size-data';
import { TransformData } from '../transform-data';
import { LocationIterator } from '../../util/location-iterator';
import { SizeThemeProps } from 'mol-view/theme/size';
import { LinesValues } from 'mol-gl/renderable/lines';

/** Wide line */
export interface Lines {
    readonly kind: 'lines',
    /** Number of lines */
    lineCount: number,
    /** Mapping buffer as array of xy values wrapped in a value cell */
    readonly mappingBuffer: ValueCell<Float32Array>,
    /** Index buffer as array of vertex index triplets wrapped in a value cell */
    readonly indexBuffer: ValueCell<Uint32Array>,
    /** Group buffer as array of group ids for each vertex wrapped in a value cell */
    readonly groupBuffer: ValueCell<Float32Array>,
    /** Line start buffer as array of xyz values wrapped in a value cell */
    readonly startBuffer: ValueCell<Float32Array>,
    /** Line end buffer as array of xyz values wrapped in a value cell */
    readonly endBuffer: ValueCell<Float32Array>,
}

export namespace Lines {
    export function createEmpty(lines?: Lines): Lines {
        const mb = lines ? lines.mappingBuffer.ref.value : new Float32Array(0)
        const ib = lines ? lines.indexBuffer.ref.value : new Uint32Array(0)
        const gb = lines ? lines.groupBuffer.ref.value : new Float32Array(0)
        const sb = lines ? lines.startBuffer.ref.value : new Float32Array(0)
        const eb = lines ? lines.endBuffer.ref.value : new Float32Array(0)
        return {
            kind: 'lines',
            lineCount: 0,
            mappingBuffer: lines ? ValueCell.update(lines.mappingBuffer, mb) : ValueCell.create(mb),
            indexBuffer: lines ? ValueCell.update(lines.indexBuffer, ib) : ValueCell.create(ib),
            groupBuffer: lines ? ValueCell.update(lines.groupBuffer, gb) : ValueCell.create(gb),
            startBuffer: lines ? ValueCell.update(lines.startBuffer, sb) : ValueCell.create(sb),
            endBuffer: lines ? ValueCell.update(lines.endBuffer, eb) : ValueCell.create(eb),
        }
    }

    export function transformImmediate(line: Lines, t: Mat4) {
        transformRangeImmediate(line, t, 0, line.lineCount)
    }

    export function transformRangeImmediate(lines: Lines, t: Mat4, offset: number, count: number) {
        const start = lines.startBuffer.ref.value
        transformPositionArray(t, start, offset, count * 4)
        ValueCell.update(lines.startBuffer, start);
        const end = lines.endBuffer.ref.value
        transformPositionArray(t, end, offset, count * 4)
        ValueCell.update(lines.endBuffer, end);
    }

    //

    export const DefaultProps = {
        ...Geometry.DefaultProps,
        lineSizeAttenuation: false,
        sizeTheme: { name: 'uniform', value: 1 } as SizeThemeProps,
    }
    export type Props = typeof DefaultProps

    export async function createValues(ctx: RuntimeContext, lines: Lines, transform: TransformData, locationIt: LocationIterator, props: Props): Promise<LinesValues> {
        const { instanceCount, groupCount } = locationIt
        const color = await createColors(ctx, locationIt, props.colorTheme)
        const size = await createSizes(ctx, locationIt, props.sizeTheme)
        const marker = createMarkers(instanceCount * groupCount)

        const counts = { drawCount: lines.lineCount * 2 * 3, groupCount, instanceCount }

        return {
            aMapping: lines.mappingBuffer,
            aGroup: lines.groupBuffer,
            aStart: lines.startBuffer,
            aEnd: lines.endBuffer,
            elements: lines.indexBuffer,
            ...color,
            ...size,
            ...marker,
            ...transform,

            ...Geometry.createValues(props, counts),
            dLineSizeAttenuation: ValueCell.create(props.lineSizeAttenuation),
            dDoubleSided: ValueCell.create(true),
            dFlipSided: ValueCell.create(false),
        }
    }

    export function updateValues(values: LinesValues, props: Props) {
        Geometry.updateValues(values, props)
        ValueCell.updateIfChanged(values.dLineSizeAttenuation, props.lineSizeAttenuation)
    }
}