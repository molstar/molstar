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
import { LinesValues } from 'mol-gl/renderable/lines';
import { Mesh } from '../mesh/mesh';
import { LinesBuilder } from './lines-builder';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { calculateBoundingSphere } from 'mol-gl/renderable/util';
import { Sphere3D } from 'mol-math/geometry';
import { Theme } from 'mol-theme/theme';

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

    export function fromMesh(mesh: Mesh, lines?: Lines) {
        const vb = mesh.vertexBuffer.ref.value
        const ib = mesh.indexBuffer.ref.value
        const gb = mesh.groupBuffer.ref.value

        const builder = LinesBuilder.create(mesh.triangleCount * 3, mesh.triangleCount / 10, lines)

        // TODO avoid duplicate lines
        for (let i = 0, il = mesh.triangleCount * 3; i < il; i += 3) {
            const i0 = ib[i], i1 = ib[i + 1], i2 = ib[i + 2];
            const x0 = vb[i0 * 3], y0 = vb[i0 * 3 + 1], z0 = vb[i0 * 3 + 2];
            const x1 = vb[i1 * 3], y1 = vb[i1 * 3 + 1], z1 = vb[i1 * 3 + 2];
            const x2 = vb[i2 * 3], y2 = vb[i2 * 3 + 1], z2 = vb[i2 * 3 + 2];
            builder.add(x0, y0, z0, x1, y1, z1, gb[i0])
            builder.add(x0, y0, z0, x2, y2, z2, gb[i0])
            builder.add(x1, y1, z1, x2, y2, z2, gb[i1])
        }

        return builder.getLines();
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

    export const Params = {
        ...Geometry.Params,
        lineSizeAttenuation: PD.Boolean(false),
    }
    export type Params = typeof Params

    export async function createValues(ctx: RuntimeContext, lines: Lines, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): Promise<LinesValues> {
        const { instanceCount, groupCount } = locationIt
        const color = await createColors(ctx, locationIt, theme.color)
        const size = await createSizes(ctx, locationIt, theme.size)
        const marker = createMarkers(instanceCount * groupCount)

        const counts = { drawCount: lines.lineCount * 2 * 3, groupCount, instanceCount }

        const boundingSphere = Sphere3D.addSphere(
            calculateBoundingSphere(
                lines.startBuffer.ref.value, lines.lineCount,
                transform.aTransform.ref.value, transform.instanceCount.ref.value
            ),
            calculateBoundingSphere(
                lines.startBuffer.ref.value, lines.lineCount,
                transform.aTransform.ref.value, transform.instanceCount.ref.value
            )
        )

        return {
            aMapping: lines.mappingBuffer,
            aGroup: lines.groupBuffer,
            aStart: lines.startBuffer,
            aEnd: lines.endBuffer,
            elements: lines.indexBuffer,
            boundingSphere: ValueCell.create(boundingSphere),
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

    export function updateValues(values: LinesValues, props: PD.Values<Params>) {
        Geometry.updateValues(values, props)
        ValueCell.updateIfChanged(values.dLineSizeAttenuation, props.lineSizeAttenuation)
    }

    export function updateBoundingSphere(values: LinesValues, lines: Lines) {
        const boundingSphere = Sphere3D.addSphere(
            calculateBoundingSphere(
                values.aStart.ref.value, lines.lineCount,
                values.aTransform.ref.value, values.instanceCount.ref.value
            ),
            calculateBoundingSphere(
                values.aEnd.ref.value, lines.lineCount,
                values.aTransform.ref.value, values.instanceCount.ref.value
            ),
        )
        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere)
        }
    }
}