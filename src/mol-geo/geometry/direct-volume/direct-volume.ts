/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { hashFnv32a } from '../../../mol-data/util';
import { LocationIterator, PositionLocation } from '../../../mol-geo/util/location-iterator';
import { RenderableState } from '../../../mol-gl/renderable';
import { DirectVolumeValues } from '../../../mol-gl/renderable/direct-volume';
import { calculateTransformBoundingSphere } from '../../../mol-gl/renderable/util';
import { createNullTexture, Texture } from '../../../mol-gl/webgl/texture';
import { Box3D, Sphere3D } from '../../../mol-math/geometry';
import { Mat4, Vec2, Vec3, Vec4 } from '../../../mol-math/linear-algebra';
import { Theme } from '../../../mol-theme/theme';
import { ValueCell } from '../../../mol-util';
import { Color } from '../../../mol-util/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { Box } from '../../primitive/box';
import { BaseGeometry } from '../base';
import { createColors } from '../color-data';
import { GeometryUtils } from '../geometry';
import { createMarkers } from '../marker-data';
import { createEmptyOverpaint } from '../overpaint-data';
import { TransformData } from '../transform-data';
import { createEmptyTransparency } from '../transparency-data';
import { createTransferFunctionTexture, getControlPointsFromVec2Array } from './transfer-function';
import { createEmptyClipping } from '../clipping-data';
import { Grid } from '../../../mol-model/volume';
import { createEmptySubstance } from '../substance-data';
import { createEmptyEmissive } from '../emissive-data';

const VolumeBox = Box();

export interface DirectVolume {
    readonly kind: 'direct-volume',

    readonly gridTexture: ValueCell<Texture>
    readonly gridTextureDim: ValueCell<Vec3>
    readonly gridDimension: ValueCell<Vec3>
    readonly gridStats: ValueCell<Vec4> // [min, max, mean, sigma]
    readonly bboxSize: ValueCell<Vec3>
    readonly bboxMin: ValueCell<Vec3>
    readonly bboxMax: ValueCell<Vec3>
    readonly transform: ValueCell<Mat4>

    readonly cellDim: ValueCell<Vec3>
    readonly unitToCartn: ValueCell<Mat4>
    readonly cartnToUnit: ValueCell<Mat4>
    readonly packedGroup: ValueCell<boolean>
    readonly axisOrder: ValueCell<Vec3>

    /** Bounding sphere of the volume */
    readonly boundingSphere: Sphere3D

    setBoundingSphere(boundingSphere: Sphere3D): void
}

export namespace DirectVolume {
    export function create(bbox: Box3D, gridDimension: Vec3, transform: Mat4, unitToCartn: Mat4, cellDim: Vec3, texture: Texture, stats: Grid['stats'], packedGroup: boolean, axisOrder: Vec3, directVolume?: DirectVolume): DirectVolume {
        return directVolume ?
            update(bbox, gridDimension, transform, unitToCartn, cellDim, texture, stats, packedGroup, axisOrder, directVolume) :
            fromData(bbox, gridDimension, transform, unitToCartn, cellDim, texture, stats, packedGroup, axisOrder);
    }

    function hashCode(directVolume: DirectVolume) {
        return hashFnv32a([
            directVolume.bboxSize.ref.version, directVolume.gridDimension.ref.version,
            directVolume.gridTexture.ref.version, directVolume.transform.ref.version,
            directVolume.gridStats.ref.version
        ]);
    }

    function fromData(bbox: Box3D, gridDimension: Vec3, transform: Mat4, unitToCartn: Mat4, cellDim: Vec3, texture: Texture, stats: Grid['stats'], packedGroup: boolean, axisOrder: Vec3): DirectVolume {
        const boundingSphere = Sphere3D();
        let currentHash = -1;

        const width = texture.getWidth();
        const height = texture.getHeight();
        const depth = texture.getDepth();

        const directVolume = {
            kind: 'direct-volume' as const,
            gridDimension: ValueCell.create(gridDimension),
            gridTexture: ValueCell.create(texture),
            gridTextureDim: ValueCell.create(Vec3.create(width, height, depth)),
            gridStats: ValueCell.create(Vec4.create(stats.min, stats.max, stats.mean, stats.sigma)),
            bboxMin: ValueCell.create(bbox.min),
            bboxMax: ValueCell.create(bbox.max),
            bboxSize: ValueCell.create(Vec3.sub(Vec3(), bbox.max, bbox.min)),
            transform: ValueCell.create(transform),
            cellDim: ValueCell.create(cellDim),
            unitToCartn: ValueCell.create(unitToCartn),
            cartnToUnit: ValueCell.create(Mat4.invert(Mat4(), unitToCartn)),
            get boundingSphere() {
                const newHash = hashCode(directVolume);
                if (newHash !== currentHash) {
                    const b = getBoundingSphere(directVolume.gridDimension.ref.value, directVolume.transform.ref.value);
                    Sphere3D.copy(boundingSphere, b);
                    currentHash = newHash;
                }
                return boundingSphere;
            },
            packedGroup: ValueCell.create(packedGroup),
            axisOrder: ValueCell.create(axisOrder),
            setBoundingSphere(sphere: Sphere3D) {
                Sphere3D.copy(boundingSphere, sphere);
                currentHash = hashCode(directVolume);
            }
        };
        return directVolume;
    }

    function update(bbox: Box3D, gridDimension: Vec3, transform: Mat4, unitToCartn: Mat4, cellDim: Vec3, texture: Texture, stats: Grid['stats'], packedGroup: boolean, axisOrder: Vec3, directVolume: DirectVolume): DirectVolume {
        const width = texture.getWidth();
        const height = texture.getHeight();
        const depth = texture.getDepth();

        ValueCell.update(directVolume.gridDimension, gridDimension);
        ValueCell.update(directVolume.gridTexture, texture);
        ValueCell.update(directVolume.gridTextureDim, Vec3.set(directVolume.gridTextureDim.ref.value, width, height, depth));
        ValueCell.update(directVolume.gridStats, Vec4.set(directVolume.gridStats.ref.value, stats.min, stats.max, stats.mean, stats.sigma));
        ValueCell.update(directVolume.bboxMin, bbox.min);
        ValueCell.update(directVolume.bboxMax, bbox.max);
        ValueCell.update(directVolume.bboxSize, Vec3.sub(directVolume.bboxSize.ref.value, bbox.max, bbox.min));
        ValueCell.update(directVolume.transform, transform);
        ValueCell.update(directVolume.cellDim, cellDim);
        ValueCell.update(directVolume.unitToCartn, unitToCartn);
        ValueCell.update(directVolume.cartnToUnit, Mat4.invert(Mat4(), unitToCartn));
        ValueCell.updateIfChanged(directVolume.packedGroup, packedGroup);
        ValueCell.updateIfChanged(directVolume.axisOrder, Vec3.fromArray(directVolume.axisOrder.ref.value, axisOrder, 0));
        return directVolume;
    }

    export function createEmpty(directVolume?: DirectVolume): DirectVolume {
        const bbox = Box3D();
        const gridDimension = Vec3();
        const transform = Mat4.identity();
        const unitToCartn = Mat4.identity();
        const cellDim = Vec3();
        const texture = createNullTexture();
        const stats = Grid.One.stats;
        const packedGroup = false;
        const axisOrder = Vec3.create(0, 1, 2);
        return create(bbox, gridDimension, transform, unitToCartn, cellDim, texture, stats, packedGroup, axisOrder, directVolume);
    }

    export const Params = {
        ...BaseGeometry.Params,
        ignoreLight: PD.Boolean(false, BaseGeometry.ShadingCategory),
        celShaded: PD.Boolean(false, BaseGeometry.ShadingCategory),
        xrayShaded: PD.Select<boolean | 'inverted'>(false, [[false, 'Off'], [true, 'On'], ['inverted', 'Inverted']], BaseGeometry.ShadingCategory),
        controlPoints: PD.LineGraph([
            Vec2.create(0.19, 0.0), Vec2.create(0.2, 0.05), Vec2.create(0.25, 0.05), Vec2.create(0.26, 0.0),
            Vec2.create(0.79, 0.0), Vec2.create(0.8, 0.05), Vec2.create(0.85, 0.05), Vec2.create(0.86, 0.0),
        ], { isEssential: true }),
        stepsPerCell: PD.Numeric(3, { min: 1, max: 10, step: 1 }),
        jumpLength: PD.Numeric(0, { min: 0, max: 20, step: 0.1 }),
    };
    export type Params = typeof Params

    export const Utils: GeometryUtils<DirectVolume, Params> = {
        Params,
        createEmpty,
        createValues,
        createValuesSimple,
        updateValues,
        updateBoundingSphere,
        createRenderableState,
        updateRenderableState,
        createPositionIterator
    };

    function createPositionIterator(directVolume: DirectVolume, transform: TransformData): LocationIterator {
        const t = directVolume.transform.ref.value;
        const [x, y, z] = directVolume.gridDimension.ref.value;
        const groupCount = x * y * z;
        const instanceCount = transform.instanceCount.ref.value;
        const location = PositionLocation();
        const p = location.position;
        const m = transform.aTransform.ref.value;
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            const k = Math.floor(groupIndex / z);
            p[0] = Math.floor(k / y);
            p[1] = k % y;
            p[2] = groupIndex % z;
            Vec3.transformMat4(p, p, t);
            if (instanceIndex >= 0) {
                Vec3.transformMat4Offset(p, p, m, 0, 0, instanceIndex * 16);
            }
            return location;
        };
        return LocationIterator(groupCount, instanceCount, 1, getLocation);
    }

    function getMaxSteps(gridDim: Vec3, stepsPerCell: number) {
        return Math.ceil(Vec3.magnitude(gridDim) * stepsPerCell);
    }

    function getStepScale(cellDim: Vec3, stepsPerCell: number) {
        return Math.min(...cellDim) * (1 / stepsPerCell);
    }

    function getTransferScale(stepsPerCell: number) {
        return (1 / stepsPerCell);
    }

    function createValues(directVolume: DirectVolume, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): DirectVolumeValues {
        const { gridTexture, gridTextureDim, gridStats } = directVolume;
        const { bboxSize, bboxMin, bboxMax, gridDimension, transform: gridTransform } = directVolume;

        const { instanceCount, groupCount } = locationIt;
        const positionIt = Utils.createPositionIterator(directVolume, transform);

        const color = createColors(locationIt, positionIt, theme.color);
        const marker = props.instanceGranularity
            ? createMarkers(instanceCount, 'instance')
            : createMarkers(instanceCount * groupCount, 'groupInstance');
        const overpaint = createEmptyOverpaint();
        const transparency = createEmptyTransparency();
        const emissive = createEmptyEmissive();
        const material = createEmptySubstance();
        const clipping = createEmptyClipping();

        const [x, y, z] = gridDimension.ref.value;
        const counts = { drawCount: VolumeBox.indices.length, vertexCount: x * y * z, groupCount, instanceCount };

        const invariantBoundingSphere = Sphere3D.clone(directVolume.boundingSphere);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, transform.aTransform.ref.value, instanceCount, 0);

        const controlPoints = getControlPointsFromVec2Array(props.controlPoints);
        const transferTex = createTransferFunctionTexture(controlPoints);

        return {
            dGeometryType: ValueCell.create('directVolume'),

            ...color,
            ...marker,
            ...overpaint,
            ...transparency,
            ...emissive,
            ...material,
            ...clipping,
            ...transform,
            ...BaseGeometry.createValues(props, counts),

            aPosition: ValueCell.create(VolumeBox.vertices as Float32Array),
            elements: ValueCell.create(VolumeBox.indices as Uint32Array),
            boundingSphere: ValueCell.create(boundingSphere),
            invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
            uInvariantBoundingSphere: ValueCell.create(Vec4.ofSphere(invariantBoundingSphere)),

            uBboxMin: bboxMin,
            uBboxMax: bboxMax,
            uBboxSize: bboxSize,
            uMaxSteps: ValueCell.create(getMaxSteps(gridDimension.ref.value, props.stepsPerCell)),
            uStepScale: ValueCell.create(getStepScale(directVolume.cellDim.ref.value, props.stepsPerCell)),
            uJumpLength: ValueCell.create(props.jumpLength),
            uTransform: gridTransform,
            uGridDim: gridDimension,
            tTransferTex: transferTex,
            uTransferScale: ValueCell.create(getTransferScale(props.stepsPerCell)),

            dGridTexType: ValueCell.create(gridTexture.ref.value.getDepth() > 0 ? '3d' : '2d'),
            uGridTexDim: gridTextureDim,
            tGridTex: gridTexture,
            uGridStats: gridStats,

            uCellDim: directVolume.cellDim,
            uCartnToUnit: directVolume.cartnToUnit,
            uUnitToCartn: directVolume.unitToCartn,
            dPackedGroup: directVolume.packedGroup,
            dAxisOrder: ValueCell.create(directVolume.axisOrder.ref.value.join('')),

            dIgnoreLight: ValueCell.create(props.ignoreLight),
            dCelShaded: ValueCell.create(props.celShaded),
            dXrayShaded: ValueCell.create(props.xrayShaded === 'inverted' ? 'inverted' : props.xrayShaded === true ? 'on' : 'off'),
        };
    }

    function createValuesSimple(directVolume: DirectVolume, props: Partial<PD.Values<Params>>, colorValue: Color, sizeValue: number, transform?: TransformData) {
        const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
        const p = { ...PD.getDefaultValues(Params), ...props };
        return createValues(directVolume, s.transform, s.locationIterator, s.theme, p);
    }

    function updateValues(values: DirectVolumeValues, props: PD.Values<Params>) {
        BaseGeometry.updateValues(values, props);
        ValueCell.updateIfChanged(values.dIgnoreLight, props.ignoreLight);
        ValueCell.updateIfChanged(values.dCelShaded, props.celShaded);
        ValueCell.updateIfChanged(values.dXrayShaded, props.xrayShaded === 'inverted' ? 'inverted' : props.xrayShaded === true ? 'on' : 'off');

        const controlPoints = getControlPointsFromVec2Array(props.controlPoints);
        createTransferFunctionTexture(controlPoints, values.tTransferTex);

        ValueCell.updateIfChanged(values.uMaxSteps, getMaxSteps(values.uGridDim.ref.value, props.stepsPerCell));
        ValueCell.updateIfChanged(values.uStepScale, getStepScale(values.uCellDim.ref.value, props.stepsPerCell));
        ValueCell.updateIfChanged(values.uTransferScale, getTransferScale(props.stepsPerCell));
        ValueCell.updateIfChanged(values.uJumpLength, props.jumpLength);
    }

    function updateBoundingSphere(values: DirectVolumeValues, directVolume: DirectVolume) {
        const invariantBoundingSphere = Sphere3D.clone(directVolume.boundingSphere);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, values.aTransform.ref.value, values.instanceCount.ref.value, 0);

        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere);
        }
        if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
            ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere);
            ValueCell.update(values.uInvariantBoundingSphere, Vec4.fromSphere(values.uInvariantBoundingSphere.ref.value, invariantBoundingSphere));
        }
    }

    function createRenderableState(props: PD.Values<Params>): RenderableState {
        const state = BaseGeometry.createRenderableState(props);
        state.opaque = false;
        state.writeDepth = false;
        return state;
    }

    function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
        BaseGeometry.updateRenderableState(state, props);
        state.opaque = false;
        state.writeDepth = false;
    }
}

//

function getBoundingSphere(gridDimension: Vec3, gridTransform: Mat4) {
    return Sphere3D.fromDimensionsAndTransform(Sphere3D(), gridDimension, gridTransform);
}