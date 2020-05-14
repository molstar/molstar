/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../../mol-util';
import { Mat4, Vec4 } from '../../../mol-math/linear-algebra';
import { transformPositionArray, GroupMapping, createGroupMapping} from '../../util';
import { GeometryUtils } from '../geometry';
import { createColors } from '../color-data';
import { createMarkers } from '../marker-data';
import { createSizes } from '../size-data';
import { TransformData } from '../transform-data';
import { LocationIterator } from '../../util/location-iterator';
import { LinesValues } from '../../../mol-gl/renderable/lines';
import { Mesh } from '../mesh/mesh';
import { LinesBuilder } from './lines-builder';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { calculateInvariantBoundingSphere, calculateTransformBoundingSphere } from '../../../mol-gl/renderable/util';
import { Sphere3D } from '../../../mol-math/geometry';
import { Theme } from '../../../mol-theme/theme';
import { Color } from '../../../mol-util/color';
import { BaseGeometry } from '../base';
import { createEmptyOverpaint } from '../overpaint-data';
import { createEmptyTransparency } from '../transparency-data';
import { hashFnv32a } from '../../../mol-data/util';
import { createEmptyClipping } from '../clipping-data';

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

    /** Bounding sphere of the lines */
    readonly boundingSphere: Sphere3D
    /** Maps group ids to line indices */
    readonly groupMapping: GroupMapping

    setBoundingSphere(boundingSphere: Sphere3D): void
}

export namespace Lines {
    export function create(mappings: Float32Array, indices: Uint32Array, groups: Float32Array, starts: Float32Array, ends: Float32Array, lineCount: number, lines?: Lines): Lines {
        return lines ?
            update(mappings, indices, groups, starts, ends, lineCount, lines) :
            fromArrays(mappings, indices, groups, starts, ends, lineCount);
    }

    export function createEmpty(lines?: Lines): Lines {
        const mb = lines ? lines.mappingBuffer.ref.value : new Float32Array(0);
        const ib = lines ? lines.indexBuffer.ref.value : new Uint32Array(0);
        const gb = lines ? lines.groupBuffer.ref.value : new Float32Array(0);
        const sb = lines ? lines.startBuffer.ref.value : new Float32Array(0);
        const eb = lines ? lines.endBuffer.ref.value : new Float32Array(0);
        return create(mb, ib, gb, sb, eb, 0, lines);
    }

    export function fromMesh(mesh: Mesh, lines?: Lines) {
        const vb = mesh.vertexBuffer.ref.value;
        const ib = mesh.indexBuffer.ref.value;
        const gb = mesh.groupBuffer.ref.value;

        const builder = LinesBuilder.create(mesh.triangleCount * 3, mesh.triangleCount / 10, lines);

        // TODO avoid duplicate lines
        for (let i = 0, il = mesh.triangleCount * 3; i < il; i += 3) {
            const i0 = ib[i], i1 = ib[i + 1], i2 = ib[i + 2];
            const x0 = vb[i0 * 3], y0 = vb[i0 * 3 + 1], z0 = vb[i0 * 3 + 2];
            const x1 = vb[i1 * 3], y1 = vb[i1 * 3 + 1], z1 = vb[i1 * 3 + 2];
            const x2 = vb[i2 * 3], y2 = vb[i2 * 3 + 1], z2 = vb[i2 * 3 + 2];
            builder.add(x0, y0, z0, x1, y1, z1, gb[i0]);
            builder.add(x0, y0, z0, x2, y2, z2, gb[i0]);
            builder.add(x1, y1, z1, x2, y2, z2, gb[i1]);
        }

        return builder.getLines();
    }

    function hashCode(lines: Lines) {
        return hashFnv32a([
            lines.lineCount, lines.mappingBuffer.ref.version, lines.indexBuffer.ref.version,
            lines.groupBuffer.ref.version, lines.startBuffer.ref.version, lines.startBuffer.ref.version
        ]);
    }

    function fromArrays(mappings: Float32Array, indices: Uint32Array, groups: Float32Array, starts: Float32Array, ends: Float32Array, lineCount: number): Lines {

        const boundingSphere = Sphere3D();
        let groupMapping: GroupMapping;

        let currentHash = -1;
        let currentGroup = -1;

        const lines = {
            kind: 'lines' as const,
            lineCount,
            mappingBuffer: ValueCell.create(mappings),
            indexBuffer: ValueCell.create(indices),
            groupBuffer: ValueCell.create(groups),
            startBuffer: ValueCell.create(starts),
            endBuffer: ValueCell.create(ends),
            get boundingSphere() {
                const newHash = hashCode(lines);
                if (newHash !== currentHash) {
                    const s = calculateInvariantBoundingSphere(lines.startBuffer.ref.value, lines.lineCount * 4, 4);
                    const e = calculateInvariantBoundingSphere(lines.endBuffer.ref.value, lines.lineCount * 4, 4);

                    Sphere3D.expandBySphere(boundingSphere, s, e);
                    currentHash = newHash;
                }
                return boundingSphere;
            },
            get groupMapping() {
                if (lines.groupBuffer.ref.version !== currentGroup) {
                    groupMapping = createGroupMapping(lines.groupBuffer.ref.value, lines.lineCount, 4);
                    currentGroup = lines.groupBuffer.ref.version;
                }
                return groupMapping;
            },
            setBoundingSphere(sphere: Sphere3D) {
                Sphere3D.copy(boundingSphere, sphere);
                currentHash = hashCode(lines);
            }
        };
        return lines;
    }

    function update(mappings: Float32Array, indices: Uint32Array, groups: Float32Array, starts: Float32Array, ends: Float32Array, lineCount: number, lines: Lines) {
        lines.lineCount = lineCount;
        ValueCell.update(lines.mappingBuffer, mappings);
        ValueCell.update(lines.indexBuffer, indices);
        ValueCell.update(lines.groupBuffer, groups);
        ValueCell.update(lines.startBuffer, starts);
        ValueCell.update(lines.endBuffer, ends);
        return lines;
    }

    export function transform(lines: Lines, t: Mat4) {
        const start = lines.startBuffer.ref.value;
        transformPositionArray(t, start, 0, lines.lineCount * 4);
        ValueCell.update(lines.startBuffer, start);
        const end = lines.endBuffer.ref.value;
        transformPositionArray(t, end, 0, lines.lineCount * 4);
        ValueCell.update(lines.endBuffer, end);
    }

    //

    export const Params = {
        ...BaseGeometry.Params,
        sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
        lineSizeAttenuation: PD.Boolean(false),
    };
    export type Params = typeof Params

    export const Utils: GeometryUtils<Lines, Params> = {
        Params,
        createEmpty,
        createValues,
        createValuesSimple,
        updateValues,
        updateBoundingSphere,
        createRenderableState: BaseGeometry.createRenderableState,
        updateRenderableState: BaseGeometry.updateRenderableState
    };

    function createValues(lines: Lines, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): LinesValues {
        const { instanceCount, groupCount } = locationIt;
        const color = createColors(locationIt, theme.color);
        const size = createSizes(locationIt, theme.size);
        const marker = createMarkers(instanceCount * groupCount);
        const overpaint = createEmptyOverpaint();
        const transparency = createEmptyTransparency();
        const clipping = createEmptyClipping();

        const counts = { drawCount: lines.lineCount * 2 * 3, groupCount, instanceCount };

        const invariantBoundingSphere = Sphere3D.clone(lines.boundingSphere);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, transform.aTransform.ref.value, instanceCount);

        return {
            aMapping: lines.mappingBuffer,
            aGroup: lines.groupBuffer,
            aStart: lines.startBuffer,
            aEnd: lines.endBuffer,
            elements: lines.indexBuffer,
            boundingSphere: ValueCell.create(boundingSphere),
            invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
            uInvariantBoundingSphere: ValueCell.create(Vec4.ofSphere(invariantBoundingSphere)),
            ...color,
            ...size,
            ...marker,
            ...overpaint,
            ...transparency,
            ...clipping,
            ...transform,

            ...BaseGeometry.createValues(props, counts),
            uSizeFactor: ValueCell.create(props.sizeFactor),
            dLineSizeAttenuation: ValueCell.create(props.lineSizeAttenuation),
            dDoubleSided: ValueCell.create(true),
            dFlipSided: ValueCell.create(false),
        };
    }

    function createValuesSimple(lines: Lines, props: Partial<PD.Values<Params>>, colorValue: Color, sizeValue: number, transform?: TransformData) {
        const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
        const p = { ...PD.getDefaultValues(Params), ...props };
        return createValues(lines, s.transform, s.locationIterator, s.theme, p);
    }

    function updateValues(values: LinesValues, props: PD.Values<Params>) {
        BaseGeometry.updateValues(values, props);
        ValueCell.updateIfChanged(values.uSizeFactor, props.sizeFactor);
        ValueCell.updateIfChanged(values.dLineSizeAttenuation, props.lineSizeAttenuation);
    }

    function updateBoundingSphere(values: LinesValues, lines: Lines) {
        const invariantBoundingSphere = Sphere3D.clone(lines.boundingSphere);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, values.aTransform.ref.value, values.instanceCount.ref.value);

        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere);
        }
        if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
            ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere);
            ValueCell.update(values.uInvariantBoundingSphere, Vec4.fromSphere(values.uInvariantBoundingSphere.ref.value, invariantBoundingSphere));
        }
    }
}