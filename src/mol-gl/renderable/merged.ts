/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../mol-util';
import { Vec2 } from '../../mol-math/linear-algebra/3d/vec2';
import { Vec4 } from '../../mol-math/linear-algebra/3d/vec4';
import { Sphere3D } from '../../mol-math/geometry/primitives/sphere3d';
import { Frustum3D } from '../../mol-math/geometry/primitives/frustum3d';
import { Plane3D } from '../../mol-math/geometry/primitives/plane3d';
import { createEmptyInstanceGrid } from '../../mol-math/geometry/instance-grid';
import { BoundaryHelper } from '../../mol-math/geometry/boundary-helper';
import { clamp } from '../../mol-math/interpolate';
import { createTextureImage, TextureImage } from './util';
import { GlobalUniformSchema, GlobalTextureSchema, GlobalDefineSchema, InternalSchema, InternalValues, SegmentSchema, SegmentValues, GlobalDefines, GlobalDefineValues, RenderableSchema, RenderableValues, BaseValues } from './schema';
import { SpheresSchema } from './spheres';
import { MeshSchema } from './mesh';
import { CylindersSchema } from './cylinders';
import { LinesSchema } from './lines';
import { PointsSchema } from './points';
import { SpheresShaderCode, MeshShaderCode, CylindersShaderCode, LinesShaderCode, PointsShaderCode, ShaderCode } from '../shader-code';
import { createGraphicsRenderItem, Transparency, GraphicsRenderVariant, DrawMode } from '../webgl/render-item';
import { WebGLContext, WebGLStats } from '../webgl/context';
import { Renderable, RenderableState, createSegmentedMdbList, CullSegment, LodLevelsValue, CullValues, createCullCache } from '../renderable';

export type MergeableValues = RenderableValues & BaseValues

export const MergedTypes = ['spheres', 'mesh', 'cylinders', 'lines', 'points'] as const;
export type MergedType = typeof MergedTypes[number]

export function isMergeableType(type: string): type is MergedType {
    return MergedTypes.includes(type as MergedType);
}

//

type Granularity = 'instance' | 'group' | 'groupInstance' | 'vertex' | 'none'

type MergedTexSpec = {
    tex: string
    dim: string
    itemSize: number
    float?: boolean
    granularity: (v: MergeableValues) => Granularity
    /** texel count multiplier, e.g. for dual colors */
    factor?: (v: MergeableValues) => number
}

function pickGranularity(type: string): Granularity {
    switch (type) {
        case 'instance': return 'instance';
        case 'group': return 'group';
        case 'groupInstance': return 'groupInstance';
        case 'vertex': return 'vertex';
        default: return 'none';
    }
}

const BaseTexSpecs: MergedTexSpec[] = [
    { tex: 'tColor', dim: 'uColorTexDim', itemSize: 3, granularity: v => pickGranularity(v.dColorType.ref.value), factor: v => v.dDualColor?.ref.value ? 2 : 1 },
    { tex: 'tSize', dim: 'uSizeTexDim', itemSize: 3, granularity: v => pickGranularity(v.dSizeType.ref.value) },
    { tex: 'tMarker', dim: 'uMarkerTexDim', itemSize: 1, granularity: v => pickGranularity(v.dMarkerType.ref.value) },
    { tex: 'tOverpaint', dim: 'uOverpaintTexDim', itemSize: 4, granularity: v => v.dOverpaint.ref.value ? pickGranularity(v.dOverpaintType.ref.value) : 'none' },
    { tex: 'tTransparency', dim: 'uTransparencyTexDim', itemSize: 1, granularity: v => v.dTransparency.ref.value ? pickGranularity(v.dTransparencyType.ref.value) : 'none' },
    { tex: 'tEmissive', dim: 'uEmissiveTexDim', itemSize: 1, granularity: v => v.dEmissive.ref.value ? pickGranularity(v.dEmissiveType.ref.value) : 'none' },
    { tex: 'tSubstance', dim: 'uSubstanceTexDim', itemSize: 4, granularity: v => v.dSubstance.ref.value ? pickGranularity(v.dSubstanceType.ref.value) : 'none' },
    { tex: 'tClipping', dim: 'uClippingTexDim', itemSize: 1, granularity: v => v.dClipping.ref.value ? pickGranularity(v.dClippingType.ref.value) : 'none' },
    { tex: 'tWiggle', dim: 'uWiggleTexDim', itemSize: 1, granularity: v => v.dWiggle.ref.value ? pickGranularity(v.dWiggleType.ref.value) : 'none' },
];

type VertexAttribute = { key: string, itemSize: number }

function getVertexAttributes(schema: RenderableSchema): VertexAttribute[] {
    const ret: VertexAttribute[] = [];
    for (const k of Object.keys(schema)) {
        const spec = schema[k];
        if (spec.type === 'attribute' && spec.divisor === 0) {
            ret.push({ key: k, itemSize: spec.itemSize });
        }
    }
    return ret;
}

type MergedTypeDescriptor = {
    readonly drawMode: DrawMode
    readonly shaderCode: ShaderCode
    readonly schema: RenderableSchema
    readonly hasElements: boolean
    /** texel count of 'vertex'-granularity textures */
    vertexTexels: (v: MergeableValues) => number
    /** kind-specific textures to merge */
    readonly extraTextures: MergedTexSpec[]
    readonly vertexAttributes: VertexAttribute[]
}

const MergedTypeDescriptors: { [k in MergedType]: MergedTypeDescriptor } = {
    spheres: {
        drawMode: 'triangles',
        shaderCode: SpheresShaderCode,
        schema: SpheresSchema,
        hasElements: false,
        // six vertices per sphere impostor share one position texel
        vertexTexels: v => v.uVertexCount.ref.value / 6,
        extraTextures: [{ tex: 'tPositionGroup', dim: 'uTexDim', itemSize: 4, float: true, granularity: () => 'vertex' }],
        vertexAttributes: getVertexAttributes(SpheresSchema),
    },
    mesh: {
        drawMode: 'triangles',
        shaderCode: MeshShaderCode,
        schema: MeshSchema,
        hasElements: true,
        vertexTexels: v => v.uVertexCount.ref.value,
        extraTextures: [],
        vertexAttributes: getVertexAttributes(MeshSchema),
    },
    cylinders: {
        drawMode: 'triangles',
        shaderCode: CylindersShaderCode,
        schema: CylindersSchema,
        hasElements: true,
        vertexTexels: v => v.uVertexCount.ref.value,
        extraTextures: [],
        vertexAttributes: getVertexAttributes(CylindersSchema),
    },
    lines: {
        drawMode: 'triangles',
        shaderCode: LinesShaderCode,
        schema: LinesSchema,
        hasElements: true,
        vertexTexels: v => v.uVertexCount.ref.value,
        extraTextures: [],
        vertexAttributes: getVertexAttributes(LinesSchema),
    },
    points: {
        drawMode: 'points',
        shaderCode: PointsShaderCode,
        schema: PointsSchema,
        hasElements: false,
        vertexTexels: v => v.uVertexCount.ref.value,
        extraTextures: [],
        vertexAttributes: getVertexAttributes(PointsSchema),
    },
};

//

function isDefineKey(k: string) {
    return k.length > 1 && k[0] === 'd' && k[1] >= 'A' && k[1] <= 'Z';
}

/**
 * Signature of everything that must be equal across members sharing a merged
 * render item: define values, granularity resolution, and lod configuration.
 *
 * Note: `lodLevels[j]`'s `scale` (index 4) is intentionally excluded - it is
 * data-dependent per member (see `Spheres.getLodLevelsValue`, based on that
 * member's own sphere count), not part of the user-configured lod preset
 * (only min/max/overlap are). Per-member scales are instead carried in the
 * `aSegmentLod` per-instance attribute (see `createMergedValues`) and
 * selected in the shader via `uLodLevel`, so they never need to match here.
 */
export function getMergeKey(values: MergeableValues): string {
    const parts: string[] = [];
    const keys = Object.keys(values).filter(isDefineKey).sort();
    for (const k of keys) {
        parts.push(`${k}:${values[k].ref.value}`);
    }
    parts.push(`ig:${values.instanceGranularity.ref.value}`);
    parts.push(`uLod:${(values.uLod.ref.value as Vec4).join(',')}`);
    const lodLevels = values.lodLevels?.ref.value as LodLevelsValue | undefined;
    parts.push('lod:' + (lodLevels ? lodLevels.map(l => `${l[0]},${l[1]},${l[2]}`).join(';') : ''));
    return parts.join('|');
}

const SupportedVertexTypes = ['uniform', 'instance', 'group', 'groupInstance', 'vertex'];
const SupportedIgTypes = ['instance', 'groupInstance'];

export function canMergeValues(type: string, members: readonly MergeableValues[]): boolean {
    if (!isMergeableType(type)) return false;
    if (members.length < 2) return false;

    const first = members[0];
    if (first.dUsePalette?.ref.value) return false;
    if (!SupportedVertexTypes.includes(first.dColorType.ref.value)) return false;
    if (first.dSizeType && !SupportedVertexTypes.includes(first.dSizeType.ref.value)) return false;
    if (first.dOverpaint.ref.value && !SupportedIgTypes.includes(first.dOverpaintType.ref.value)) return false;
    if (first.dTransparency.ref.value && !SupportedIgTypes.includes(first.dTransparencyType.ref.value)) return false;
    if (first.dEmissive.ref.value && !SupportedIgTypes.includes(first.dEmissiveType.ref.value)) return false;
    if (first.dSubstance.ref.value && !SupportedIgTypes.includes(first.dSubstanceType.ref.value)) return false;
    if (first.dClipping.ref.value && !SupportedIgTypes.includes(first.dClippingType.ref.value)) return false;
    if (first.dWiggle.ref.value && !SupportedIgTypes.includes(first.dWiggleType.ref.value)) return false;

    const key = getMergeKey(first);
    for (let i = 1, il = members.length; i < il; ++i) {
        if (getMergeKey(members[i]) !== key) return false;
    }
    return true;
}

export type MergedSegments = {
    count: number
    /** first instance (into the merged instanced attributes) per segment */
    instanceBases: number[]
    /** instance count per segment */
    instanceCounts: number[]
    /** first vertex (into the merged geometry) per segment */
    vertexBases: number[]
    /** byte offset (into the merged element buffer) per segment */
    elementOffsets: number[]
}

export type Merged = {
    readonly type: MergedType
    readonly values: MergeableValues & SegmentValues
    readonly segments: MergedSegments
    /** sync merged values with member values, returns true if the layout changed */
    sync: () => boolean
}

function ensureArray<T extends Uint8Array | Uint32Array | Float32Array>(array: T, length: number, ctor: new (length: number) => T): T {
    return array.length >= length ? array : new ctor(length);
}

const boundaryHelper = new BoundaryHelper('98');

export function createMergedValues(type: MergedType, members: readonly MergeableValues[]): Merged {
    const descriptor = MergedTypeDescriptors[type];
    const { hasElements, vertexAttributes } = descriptor;
    const first = members[0];
    const n = members.length;

    const texSpecs = [...descriptor.extraTextures, ...BaseTexSpecs].filter(spec => first[spec.tex] !== undefined);

    const segments: MergedSegments = {
        count: n,
        instanceBases: new Array(n).fill(0),
        instanceCounts: new Array(n).fill(0),
        vertexBases: new Array(n).fill(0),
        elementOffsets: new Array(n).fill(0),
    };

    // own cells, identity must stay stable for the lifetime of the merged render item
    const own: { [k: string]: ValueCell<any> } = {
        aTransform: ValueCell.create(new Float32Array(0)),
        aInstance: ValueCell.create(new Float32Array(0)),
        aSegment: ValueCell.create(new Float32Array(0)),
        aSegmentSphere: ValueCell.create(new Float32Array(0)),
        aSegmentLod: ValueCell.create(new Float32Array(0)),
        drawCount: ValueCell.create(0),
        instanceCount: ValueCell.create(0),
        uVertexCount: ValueCell.create(0),
        uInstanceCount: ValueCell.create(0),
        uGroupCount: ValueCell.create(0),
        boundingSphere: ValueCell.create(Sphere3D()),
        invariantBoundingSphere: ValueCell.create(Sphere3D()),
        uInvariantBoundingSphere: ValueCell.create(Vec4()),
        uMarker: ValueCell.create(0),
        markerAverage: ValueCell.create(0),
        markerStatus: ValueCell.create(-1),
        transparencyAverage: ValueCell.create(0),
        transparencyMin: ValueCell.create(0),
        emissiveAverage: ValueCell.create(0),
        wiggleAverage: ValueCell.create(0),
        hasReflection: ValueCell.create(false),
        instanceGrid: ValueCell.create({ ...createEmptyInstanceGrid(), cellSize: 2 }),
        dSegmented: ValueCell.create(true),
    };
    for (const { key } of vertexAttributes) {
        own[key] = ValueCell.create(new Float32Array(0));
    }
    if (hasElements) {
        own.elements = ValueCell.create(new Uint32Array(0));
    }
    for (const spec of texSpecs) {
        const ctor = (spec.float ? Float32Array : Uint8Array) as new (length: number) => any;
        own[spec.tex] = ValueCell.create(createTextureImage(1, spec.itemSize, ctor));
        own[spec.dim] = ValueCell.create(Vec2.create(1, 1));
    }

    const values = {
        ...first,
        ...own,
    } as unknown as MergeableValues & SegmentValues;

    /**
     * Extent of a member's 'group'-granularity data. `uGroupCount` comes from the
     * visual's location iterator and can understate it - e.g. particle target
     * structures iterate per particle (groupCount 1) but carry per-element sizes.
     */
    function getGroupSpan(i: number): number {
        const v = members[i];
        let span = v.uGroupCount.ref.value;
        for (const spec of texSpecs) {
            if (spec.granularity(first) !== 'group') continue;
            const image = v[spec.tex].ref.value as TextureImage<any>;
            const texels = image.width * image.height;
            if (texels > span) span = texels;
        }
        return span;
    }

    function getTexelCount(i: number, spec: MergedTexSpec): number {
        const v = members[i];
        const factor = spec.factor ? spec.factor(v) : 1;
        switch (spec.granularity(first)) {
            case 'instance': return v.instanceCount.ref.value * factor;
            case 'group': return getGroupSpan(i) * factor;
            case 'groupInstance': return v.instanceCount.ref.value * v.uGroupCount.ref.value * factor;
            case 'vertex': return descriptor.vertexTexels(v) * factor;
            case 'none': return 0;
        }
    }

    /** texel base per member for each merged texture */
    const texBases = new Map<string, number[]>();
    for (const spec of texSpecs) {
        texBases.set(spec.tex, new Array(n).fill(0));
    }

    const memberVersions: { [key: string]: number }[] = members.map(() => ({}));
    let layoutKey = '';

    function getLayoutKey() {
        let key = '';
        for (let i = 0; i < n; ++i) {
            const v = members[i];
            key += `${v.instanceCount.ref.value},${v.uGroupCount.ref.value},${getGroupSpan(i)},${v.uVertexCount.ref.value},${v.drawCount.ref.value};`;
        }
        return key + '#' + getMergeKey(first);
    }

    function copyTexBlock(image: TextureImage<any>, spec: MergedTexSpec, i: number, base: number) {
        const texels = getTexelCount(i, spec);
        const src = members[i][spec.tex].ref.value as TextureImage<any>;
        image.array.set(src.array.subarray(0, texels * spec.itemSize), base * spec.itemSize);
    }

    function mergeTextureFull(spec: MergedTexSpec) {
        const cell = own[spec.tex];
        const dimCell = own[spec.dim];
        const bases = texBases.get(spec.tex)!;
        const ctor = (spec.float ? Float32Array : Uint8Array) as new (length: number) => any;

        if (spec.granularity(first) === 'none') {
            // no per-segment data used, mirror the content of the first member
            const src = first[spec.tex].ref.value as TextureImage<any>;
            const array = ensureArray(cell.ref.value.array, src.array.length, ctor);
            array.set(src.array);
            bases.fill(0);
            ValueCell.update(cell, { array, width: src.width, height: src.height });
            ValueCell.update(dimCell, Vec2.set(dimCell.ref.value, src.width, src.height));
            return;
        }

        let total = 0;
        for (let i = 0; i < n; ++i) {
            bases[i] = total;
            total += getTexelCount(i, spec);
        }
        const image = createTextureImage(Math.max(1, total), spec.itemSize, ctor, cell.ref.value.array);
        for (let i = 0; i < n; ++i) {
            copyTexBlock(image, spec, i, bases[i]);
        }
        ValueCell.update(cell, image);
        ValueCell.update(dimCell, Vec2.set(dimCell.ref.value, image.width, image.height));
    }

    function unionSpheres(out: Sphere3D, key: 'boundingSphere' | 'invariantBoundingSphere') {
        boundaryHelper.reset();
        for (let i = 0; i < n; ++i) {
            boundaryHelper.includeSphere(members[i][key].ref.value as Sphere3D);
        }
        boundaryHelper.finishedIncludeStep();
        for (let i = 0; i < n; ++i) {
            boundaryHelper.radiusSphere(members[i][key].ref.value as Sphere3D);
        }
        return boundaryHelper.getSphere(out);
    }

    function updateBoundingSpheres() {
        ValueCell.update(own.boundingSphere, unionSpheres(own.boundingSphere.ref.value, 'boundingSphere'));
        const ibs = unionSpheres(own.invariantBoundingSphere.ref.value, 'invariantBoundingSphere');
        ValueCell.update(own.invariantBoundingSphere, ibs);
        ValueCell.update(own.uInvariantBoundingSphere, Vec4.fromSphere(own.uInvariantBoundingSphere.ref.value, ibs));
    }

    function updateMarkerAggregates() {
        let average = 0;
        let weight = 0;
        let status = members[0].markerStatus.ref.value;
        for (let i = 0; i < n; ++i) {
            const v = members[i];
            const w = Math.max(1, v.dMarkerType.ref.value === 'instance'
                ? v.instanceCount.ref.value
                : v.instanceCount.ref.value * v.uGroupCount.ref.value);
            average += v.markerAverage.ref.value * w;
            weight += w;
            if (v.markerStatus.ref.value !== status) status = -1;
        }
        ValueCell.updateIfChanged(own.markerAverage, weight > 0 ? average / weight : 0);
        ValueCell.updateIfChanged(own.markerStatus, status);
        ValueCell.updateIfChanged(own.uMarker, status);
    }

    function updateAggregates() {
        let transparencyAverage = 0;
        let transparencyMin = 1;
        let emissiveAverage = 0;
        let wiggleAverage = 0;
        let hasReflection = false;
        let weight = 0;
        for (let i = 0; i < n; ++i) {
            const v = members[i];
            const w = Math.max(1, v.instanceCount.ref.value * v.uGroupCount.ref.value);
            transparencyAverage += v.transparencyAverage.ref.value * w;
            emissiveAverage += v.emissiveAverage.ref.value * w;
            wiggleAverage += v.wiggleAverage.ref.value * w;
            transparencyMin = Math.min(transparencyMin, v.transparencyMin.ref.value);
            hasReflection = hasReflection || v.hasReflection.ref.value;
            weight += w;
        }
        ValueCell.updateIfChanged(own.transparencyAverage, weight > 0 ? transparencyAverage / weight : 0);
        ValueCell.updateIfChanged(own.transparencyMin, transparencyMin);
        ValueCell.updateIfChanged(own.emissiveAverage, weight > 0 ? emissiveAverage / weight : 0);
        ValueCell.updateIfChanged(own.wiggleAverage, weight > 0 ? wiggleAverage / weight : 0);
        ValueCell.updateIfChanged(own.hasReflection, hasReflection);
    }

    function copyVertexBlock(attr: VertexAttribute, i: number) {
        const vertexCount = members[i].uVertexCount.ref.value;
        const array = own[attr.key].ref.value as Float32Array;
        array.set((members[i][attr.key].ref.value as Float32Array).subarray(0, vertexCount * attr.itemSize), segments.vertexBases[i] * attr.itemSize);
    }

    function copyElementsBlock(i: number) {
        const drawCount = members[i].drawCount.ref.value;
        const vertexBase = segments.vertexBases[i];
        const elementBase = segments.elementOffsets[i] / 4;
        const src = members[i].elements.ref.value as Uint32Array;
        const dst = own.elements.ref.value as Uint32Array;
        for (let k = 0; k < drawCount; ++k) {
            dst[elementBase + k] = src[k] + vertexBase;
        }
    }

    /** per-segment values for the aSegment attribute */
    const segmentGroupCounts = new Array(n).fill(0);
    const segmentGroupInstanceOffsets = new Array(n).fill(0);
    const segmentGroupBases = new Array(n).fill(0);

    /** fill member i's instance block of aSegmentLod with its own (sphere-count-dependent) per-level scales */
    function copySegmentLodBlock(i: number, array: Float32Array) {
        const instanceCount = segments.instanceCounts[i];
        const instanceBase = segments.instanceBases[i];
        const member = members[i];
        // gridless members are drawn at full (undecimated) detail, see the
        // `cellSize <= 1` branch in renderable.ts's cullSegment - no per-level
        // scale compensation applies then, else spheres would render too large
        const gridless = member.instanceGrid.ref.value.cellSize <= 1;
        const lodLevels = !gridless ? member.lodLevels?.ref.value as LodLevelsValue | undefined : undefined;
        const s0 = lodLevels?.[0]?.[4] ?? 1;
        const s1 = lodLevels?.[1]?.[4] ?? 1;
        const s2 = lodLevels?.[2]?.[4] ?? 1;
        const s3 = lodLevels?.[3]?.[4] ?? 1;
        for (let k = 0; k < instanceCount; ++k) {
            const o = (instanceBase + k) * 4;
            array[o + 0] = s0;
            array[o + 1] = s1;
            array[o + 2] = s2;
            array[o + 3] = s3;
        }
    }

    function computeSegments() {
        let instanceTotal = 0;
        let vertexTotal = 0;
        let drawTotal = 0;
        let maxGroupCount = 0;
        let groupInstanceBase = 0;
        let groupBase = 0;
        for (let i = 0; i < n; ++i) {
            const v = members[i];
            const instanceCount = v.instanceCount.ref.value;
            const groupCount = v.uGroupCount.ref.value;
            segments.instanceBases[i] = instanceTotal;
            segments.instanceCounts[i] = instanceCount;
            segments.vertexBases[i] = vertexTotal;
            segments.elementOffsets[i] = drawTotal * 4;
            segmentGroupCounts[i] = groupCount;
            segmentGroupInstanceOffsets[i] = groupInstanceBase - instanceTotal * groupCount;
            segmentGroupBases[i] = groupBase;
            instanceTotal += instanceCount;
            vertexTotal += v.uVertexCount.ref.value;
            drawTotal += v.drawCount.ref.value;
            groupInstanceBase += instanceCount * groupCount;
            groupBase += getGroupSpan(i);
            if (groupCount > maxGroupCount) maxGroupCount = groupCount;
        }
        segments.count = n;
        return { instanceTotal, vertexTotal, drawTotal, maxGroupCount };
    }

    function rebuildInstances(instanceTotal: number) {
        const segmentArray = ensureArray(own.aSegment.ref.value, instanceTotal * 3, Float32Array);
        const segmentSphereArray = ensureArray(own.aSegmentSphere.ref.value, instanceTotal * 4, Float32Array);
        const segmentLodArray = ensureArray(own.aSegmentLod.ref.value, instanceTotal * 4, Float32Array);

        // instanced attributes in per-segment blocks
        const transformArray = ensureArray(own.aTransform.ref.value, instanceTotal * 16, Float32Array);
        const instanceArray = ensureArray(own.aInstance.ref.value, instanceTotal, Float32Array);
        for (let i = 0; i < n; ++i) {
            const v = members[i];
            const instanceCount = segments.instanceCounts[i];
            const instanceBase = segments.instanceBases[i];
            transformArray.set((v.aTransform.ref.value as Float32Array).subarray(0, instanceCount * 16), instanceBase * 16);
            const memberInstance = v.aInstance.ref.value as Float32Array;
            const ibs = v.invariantBoundingSphere.ref.value as Sphere3D;
            copySegmentLodBlock(i, segmentLodArray);
            for (let k = 0; k < instanceCount; ++k) {
                instanceArray[instanceBase + k] = memberInstance[k] + instanceBase;
                segmentArray[(instanceBase + k) * 3 + 0] = segmentGroupCounts[i];
                segmentArray[(instanceBase + k) * 3 + 1] = segmentGroupInstanceOffsets[i];
                segmentArray[(instanceBase + k) * 3 + 2] = segmentGroupBases[i];
                segmentSphereArray[(instanceBase + k) * 4 + 0] = ibs.center[0];
                segmentSphereArray[(instanceBase + k) * 4 + 1] = ibs.center[1];
                segmentSphereArray[(instanceBase + k) * 4 + 2] = ibs.center[2];
                segmentSphereArray[(instanceBase + k) * 4 + 3] = ibs.radius;
            }
        }
        ValueCell.update(own.aTransform, transformArray);
        ValueCell.update(own.aInstance, instanceArray);
        // cellSize > 1 keeps the renderer on the cull path, each member's OWN instanceGrid is used instead
        ValueCell.update(own.instanceGrid, { ...createEmptyInstanceGrid(), cellSize: 2 });
        ValueCell.update(own.aSegment, segmentArray);
        ValueCell.update(own.aSegmentSphere, segmentSphereArray);
        ValueCell.update(own.aSegmentLod, segmentLodArray);
    }

    function rebuild() {
        const { instanceTotal, vertexTotal, drawTotal, maxGroupCount } = computeSegments();

        rebuildInstances(instanceTotal);

        for (const attr of vertexAttributes) {
            ValueCell.update(own[attr.key], ensureArray(own[attr.key].ref.value, vertexTotal * attr.itemSize, Float32Array));
            for (let i = 0; i < n; ++i) copyVertexBlock(attr, i);
            ValueCell.update(own[attr.key], own[attr.key].ref.value);
        }
        if (hasElements) {
            ValueCell.update(own.elements, ensureArray(own.elements.ref.value, drawTotal, Uint32Array));
            for (let i = 0; i < n; ++i) copyElementsBlock(i);
            ValueCell.update(own.elements, own.elements.ref.value);
        }

        for (const spec of texSpecs) {
            mergeTextureFull(spec);
        }

        ValueCell.updateIfChanged(own.drawCount, drawTotal);
        ValueCell.updateIfChanged(own.instanceCount, instanceTotal);
        ValueCell.updateIfChanged(own.uVertexCount, vertexTotal);
        ValueCell.updateIfChanged(own.uInstanceCount, instanceTotal);
        ValueCell.updateIfChanged(own.uGroupCount, maxGroupCount);

        updateBoundingSpheres();
        updateMarkerAggregates();
        updateAggregates();
    }

    function snapshotAll() {
        for (let i = 0; i < n; ++i) {
            const v = members[i];
            const mv = memberVersions[i];
            mv.aTransform = v.aTransform.ref.version;
            mv.aInstance = v.aInstance.ref.version;
            mv.boundingSphere = v.boundingSphere.ref.version;
            mv.invariantBoundingSphere = v.invariantBoundingSphere.ref.version;
            mv.markerStatus = v.markerStatus.ref.version;
            mv.markerAverage = v.markerAverage.ref.version;
            mv.transparencyAverage = v.transparencyAverage.ref.version;
            mv.transparencyMin = v.transparencyMin.ref.version;
            mv.emissiveAverage = v.emissiveAverage.ref.version;
            mv.wiggleAverage = v.wiggleAverage.ref.version;
            mv.hasReflection = v.hasReflection.ref.version;
            mv.lodLevels = v.lodLevels ? v.lodLevels.ref.version : -1;
            for (const attr of vertexAttributes) {
                mv[attr.key] = v[attr.key].ref.version;
            }
            if (hasElements) {
                mv.elements = v.elements.ref.version;
            }
            for (const spec of texSpecs) {
                mv[spec.tex] = v[spec.tex].ref.version;
            }
        }
    }

    function sync() {
        const key = getLayoutKey();
        if (key !== layoutKey) {
            layoutKey = key;
            rebuild();
            snapshotAll();
            return true;
        }

        let instancesChanged = false;
        let boundsChanged = false;
        let markersChanged = false;
        let aggregatesChanged = false;
        let elementsChanged = false;
        let lodChanged = false;
        const attributesChanged = new Set<string>();

        for (let i = 0; i < n; ++i) {
            const v = members[i];
            const mv = memberVersions[i];

            // the shared grid placement depends on transforms and invariant spheres
            if (mv.aTransform !== v.aTransform.ref.version ||
                mv.aInstance !== v.aInstance.ref.version ||
                mv.invariantBoundingSphere !== v.invariantBoundingSphere.ref.version
            ) {
                instancesChanged = true;
            }
            if (mv.boundingSphere !== v.boundingSphere.ref.version || mv.invariantBoundingSphere !== v.invariantBoundingSphere.ref.version) {
                boundsChanged = true;
            }
            if (mv.markerStatus !== v.markerStatus.ref.version || mv.markerAverage !== v.markerAverage.ref.version) {
                markersChanged = true;
            }
            if (mv.transparencyAverage !== v.transparencyAverage.ref.version ||
                mv.transparencyMin !== v.transparencyMin.ref.version ||
                mv.emissiveAverage !== v.emissiveAverage.ref.version ||
                mv.wiggleAverage !== v.wiggleAverage.ref.version ||
                mv.hasReflection !== v.hasReflection.ref.version
            ) {
                aggregatesChanged = true;
            }
            for (const attr of vertexAttributes) {
                if (mv[attr.key] !== v[attr.key].ref.version) {
                    copyVertexBlock(attr, i);
                    attributesChanged.add(attr.key);
                }
            }
            if (hasElements && mv.elements !== v.elements.ref.version) {
                copyElementsBlock(i);
                elementsChanged = true;
            }

            const lodVersion = v.lodLevels ? v.lodLevels.ref.version : -1;
            if (mv.lodLevels !== lodVersion) {
                copySegmentLodBlock(i, own.aSegmentLod.ref.value as Float32Array);
                lodChanged = true;
            }
        }

        for (const spec of texSpecs) {
            const bases = texBases.get(spec.tex)!;
            let texChanged = false;
            if (spec.granularity(first) === 'none') {
                if (memberVersions[0][spec.tex] !== first[spec.tex].ref.version) {
                    mergeTextureFull(spec);
                }
                continue;
            }
            const image = own[spec.tex].ref.value as TextureImage<any>;
            for (let i = 0; i < n; ++i) {
                const v = members[i];
                const changed = memberVersions[i][spec.tex] !== v[spec.tex].ref.version ||
                    // uniform-status marking updates the array without bumping the texture cell
                    (spec.tex === 'tMarker' && memberVersions[i].markerStatus !== v.markerStatus.ref.version);
                if (changed) {
                    copyTexBlock(image, spec, i, bases[i]);
                    texChanged = true;
                }
            }
            if (texChanged) {
                ValueCell.update(own[spec.tex], image);
            }
        }

        if (instancesChanged) rebuildInstances(own.instanceCount.ref.value as number);
        for (const key of attributesChanged) ValueCell.update(own[key], own[key].ref.value);
        if (elementsChanged) ValueCell.update(own.elements, own.elements.ref.value);
        if (!instancesChanged && lodChanged) ValueCell.update(own.aSegmentLod, own.aSegmentLod.ref.value);
        if (boundsChanged) updateBoundingSpheres();
        if (markersChanged) updateMarkerAggregates();
        if (aggregatesChanged) updateAggregates();

        snapshotAll();
        return false;
    }

    layoutKey = getLayoutKey();
    rebuild();
    snapshotAll();

    return { type, values, segments, sync };
}

//

export function MergedRenderable(ctx: WebGLContext, id: number, merged: Merged, members: readonly MergeableValues[], state: RenderableState, materialId: number, transparency: Transparency, globals: GlobalDefines): Renderable<MergeableValues & SegmentValues> {
    const descriptor = MergedTypeDescriptors[merged.type];
    const schema = { ...GlobalUniformSchema, ...GlobalTextureSchema, ...GlobalDefineSchema, ...InternalSchema, ...descriptor.schema, ...SegmentSchema };
    const renderValues: MergeableValues & SegmentValues & InternalValues & GlobalDefineValues = {
        ...merged.values,
        uObjectId: ValueCell.create(id),
        dLightCount: ValueCell.create(globals.dLightCount),
        dColorMarker: ValueCell.create(globals.dColorMarker),
    };
    const renderItem = createGraphicsRenderItem(ctx, descriptor.drawMode, descriptor.shaderCode, schema, renderValues, materialId, transparency);

    const mdb = createSegmentedMdbList();
    let mode: 'none' | 'full' | 'cull' = 'none';
    const cullCache = createCullCache();

    const segment: CullSegment = { first: 0, offset: 0, instanceBase: 0 };
    function setSegment(i: number) {
        const { segments } = merged;
        segment.first = descriptor.hasElements ? 0 : segments.vertexBases[i];
        segment.offset = segments.elementOffsets[i];
        segment.instanceBase = segments.instanceBases[i];
        return segment;
    }

    function getCapacity() {
        let capacity = 0;
        for (let i = 0, il = members.length; i < il; ++i) {
            capacity += Math.max(1, (members[i] as CullValues).instanceGrid.ref.value.cellCount) + 1;
        }
        return capacity;
    }

    function getMemberLod(i: number, hasLodLevels: boolean): LodLevelsValue | undefined {
        // per-member lod counts, the level ranges are equal across members
        return hasLodLevels ? members[i].lodLevels!.ref.value as LodLevelsValue : undefined;
    }

    function buildFullList() {
        const hasLodLevels = !!mdb.prepare(members[0].lodLevels, getCapacity(), true);
        for (let i = 0, il = members.length; i < il; ++i) {
            mdb.fullSegment(members[i] as CullValues, setSegment(i), getMemberLod(i, hasLodLevels));
        }
        mode = 'full';
    }

    function cull(cameraPlane: Plane3D, frustum: Frustum3D, isOccluded: ((s: Sphere3D) => boolean) | null, stats: WebGLStats) {
        // skip recomputation if nothing relevant to culling has actually changed
        if (cullCache.unchanged(merged.values, cameraPlane, frustum, isOccluded)) return;

        const lodLevels = mdb.prepare(members[0].lodLevels, getCapacity(), true);
        mode = 'cull';

        for (let i = 0, il = members.length; i < il; ++i) {
            mdb.cullSegment(members[i] as CullValues, setSegment(i), getMemberLod(i, !!lodLevels), cameraPlane, frustum, isOccluded, stats, true);
        }

        cullCache.commit(merged.values, cameraPlane, frustum, isOccluded);
    }

    function cullSimple(d: number, radius: number, scale: number) {
        const lodLevels = mdb.prepare(members[0].lodLevels, getCapacity(), true);
        if (!lodLevels) {
            buildFullList();
            return;
        }
        mode = 'cull';

        for (let i = 0, il = members.length; i < il; ++i) {
            mdb.cullSimpleSegment(members[i] as CullValues, setSegment(i), getMemberLod(i, true)!, d, radius, scale);
        }
    }

    buildFullList();

    return {
        id,
        materialId,
        values: renderValues,
        state,

        cull,
        uncull: () => {
            cullCache.invalidate();
            if (mode !== 'full') buildFullList();
        },
        cullSimple,
        render: (variant: GraphicsRenderVariant, sharedTexturesCount: number) => {
            if (renderValues.uAlpha && renderValues.alpha) {
                ValueCell.updateIfChanged(renderValues.uAlpha, clamp(renderValues.alpha.ref.value * state.alphaFactor, 0, 1));
            }
            renderItem.render(variant, sharedTexturesCount, mdb.list);
        },
        getByteCount: () => renderItem.getByteCount(),
        getProgram: (variant: GraphicsRenderVariant) => renderItem.getProgram(variant),
        setTransparency: (value: Transparency) => renderItem.setTransparency(value),
        update: () => {
            merged.sync();
            mode = 'none';
            cullCache.invalidate();
            renderItem.update();
        },
        dispose: () => renderItem.destroy(),
    };
}
