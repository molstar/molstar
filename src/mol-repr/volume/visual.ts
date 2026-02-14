/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Visual, VisualContext } from '../visual';
import { Grid, Volume } from '../../mol-model/volume';
import { Geometry, GeometryUtils } from '../../mol-geo/geometry/geometry';
import { LocationIterator } from '../../mol-geo/util/location-iterator';
import { Theme } from '../../mol-theme/theme';
import { createIdentityTransform, createTransform, TransformData } from '../../mol-geo/geometry/transform-data';
import { createRenderObject, GraphicsRenderObject } from '../../mol-gl/render-object';
import { PickingId } from '../../mol-geo/geometry/picking';
import { Loci, isEveryLoci, EmptyLoci } from '../../mol-model/loci';
import { Interval, OrderedSet } from '../../mol-data/int';
import { LocationCallback, VisualUpdateState } from '../util';
import { ColorTheme } from '../../mol-theme/color';
import { ValueCell } from '../../mol-util';
import { createSizes } from '../../mol-geo/geometry/size-data';
import { createColors } from '../../mol-geo/geometry/color-data';
import { MarkerAction } from '../../mol-util/marker-action';
import { Mat4 } from '../../mol-math/linear-algebra';
import { Overpaint } from '../../mol-theme/overpaint';
import { Transparency } from '../../mol-theme/transparency';
import { SizeValues } from '../../mol-gl/renderable/schema';
import { Clipping } from '../../mol-theme/clipping';
import { isPromiseLike } from '../../mol-util/type-helpers';
import { Substance } from '../../mol-theme/substance';
import { createMarkers } from '../../mol-geo/geometry/marker-data';
import { Emissive } from '../../mol-theme/emissive';
import { SizeTheme } from '../../mol-theme/size';
import { Sphere3D } from '../../mol-math/geometry/primitives/sphere3d';
import { BaseGeometry } from '../../mol-geo/geometry/base';

export const VolumeParams = {
    ...BaseGeometry.Params,
};
export type VolumeParams = typeof VolumeParams

export type VolumeKey = { volume: Volume, key: number }
export interface VolumeVisual<P extends VolumeParams> extends Visual<VolumeKey, P> { }

function createVolumeInstancesTransform(volume: Volume, invariantBoundingSphere: Sphere3D, cellSize: number, batchSize: number, transformData?: TransformData) {
    const instanceCount = volume.instances.length;
    const transformArray = new Float32Array(instanceCount * 16);
    for (let i = 0; i < instanceCount; ++i) {
        Mat4.toArray(volume.instances[i].transform, transformArray, i * 16);
    }
    return createTransform(transformArray, instanceCount, invariantBoundingSphere, cellSize, batchSize, transformData);
}

function createVolumeRenderObject<G extends Geometry>(volume: Volume, geometry: G, locationIt: LocationIterator, theme: Theme, props: PD.Values<Geometry.Params<G>>, materialId: number) {
    const { createValues, createRenderableState } = Geometry.getUtils(geometry);
    const transform = locationIt.nonInstanceable
        ? createIdentityTransform()
        : createVolumeInstancesTransform(volume, geometry.boundingSphere, props.cellSize, props.batchSize);
    const values = createValues(geometry, transform, locationIt, theme, props);
    const state = createRenderableState(props);
    return createRenderObject(geometry.kind, values, state, materialId);
}

interface VolumeVisualBuilder<P extends VolumeParams, G extends Geometry> {
    defaultProps: PD.Values<P>
    createGeometry(ctx: VisualContext, volume: Volume, key: number, theme: Theme, props: PD.Values<P>, geometry?: G): Promise<G> | G
    createLocationIterator(volume: Volume, key: number): LocationIterator
    getLoci(pickingId: PickingId, volume: Volume, key: number, props: PD.Values<P>, id: number): Loci
    eachLocation(loci: Loci, volume: Volume, key: number, props: PD.Values<P>, apply: (interval: Interval) => boolean): boolean
    setUpdateState(state: VisualUpdateState, volume: Volume, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme): void
    initUpdateState?: (state: VisualUpdateState, volume: Volume, newProps: PD.Values<P>, newTheme: Theme) => void
    mustRecreate?: (volumeKey: VolumeKey, props: PD.Values<P>) => boolean
    dispose?: (geometry: G) => void
}

interface VolumeVisualGeometryBuilder<P extends VolumeParams, G extends Geometry> extends VolumeVisualBuilder<P, G> {
    geometryUtils: GeometryUtils<G>
}

export function VolumeVisual<G extends Geometry, P extends VolumeParams & Geometry.Params<G>>(builder: VolumeVisualGeometryBuilder<P, G>, materialId: number): VolumeVisual<P> {
    const { defaultProps, createGeometry, createLocationIterator, getLoci, eachLocation, setUpdateState, initUpdateState, mustRecreate, dispose } = builder;
    const { updateValues, updateBoundingSphere, updateRenderableState, createPositionIterator } = builder.geometryUtils;
    const updateState = VisualUpdateState.create();

    let renderObject: GraphicsRenderObject<G['kind']> | undefined;

    let newProps: PD.Values<P>;
    let newTheme: Theme;
    let newVolume: Volume;
    let newKey: number;

    let currentProps: PD.Values<P> = Object.assign({}, defaultProps);
    let currentTheme: Theme = Theme.createEmpty();
    let currentVolume: Volume;
    let currentKey: number;

    let geometry: G;
    let geometryVersion = -1;
    let locationIt: LocationIterator;
    let positionIt: LocationIterator;

    function prepareUpdate(theme: Theme, props: Partial<PD.Values<P>>, volume: Volume, key: number) {
        if (!volume && !currentVolume) {
            throw new Error('missing volume');
        }

        newProps = Object.assign({}, currentProps, props);
        newTheme = theme;
        newVolume = volume;
        newKey = key;

        VisualUpdateState.reset(updateState);

        if (!renderObject) {
            updateState.createNew = true;
        } else if (Grid.areEquivalent(newVolume.grid, currentVolume.grid) && !Volume.areInstanceTransformsEqual(newVolume, currentVolume)) {
            updateState.updateTransform = true;
        } else if (!Volume.areEquivalent(newVolume, currentVolume) || newKey !== currentKey) {
            updateState.createNew = true;
        }

        if (updateState.createNew) {
            initUpdateState?.(updateState, newVolume, newProps, newTheme);
            updateState.createGeometry = true;
            return;
        }

        setUpdateState(updateState, volume, newProps, currentProps, newTheme, currentTheme);

        if (!ColorTheme.areEqual(theme.color, currentTheme.color)) {
            updateState.updateColor = true;
        }

        if (!SizeTheme.areEqual(theme.size, currentTheme.size)) {
            updateState.updateSize = true;
        }

        if (locationIt.nonInstanceable) {
            if (newProps.instanceGranularity !== currentProps.instanceGranularity) {
                updateState.updateTransform = true;
            }
        } else {
            if (newProps.instanceGranularity !== currentProps.instanceGranularity || newProps.cellSize !== currentProps.cellSize || newProps.batchSize !== currentProps.batchSize) {
                updateState.updateTransform = true;
            }

            if (updateState.updateTransform) {
                updateState.updateMatrix = true;
            }
        }

        if (updateState.updateSize && !('uSize' in renderObject!.values)) {
            updateState.createGeometry = true;
        }

        if (updateState.createGeometry) {
            updateState.updateColor = true;
            updateState.updateSize = true;
        }
    }

    function update(newGeometry?: G) {
        if (updateState.createNew) {
            locationIt = createLocationIterator(newVolume, newKey);
            if (newGeometry) {
                renderObject = createVolumeRenderObject(newVolume, newGeometry, locationIt, newTheme, newProps, materialId);
                positionIt = createPositionIterator(newGeometry, renderObject.values);
            } else {
                throw new Error('expected geometry to be given');
            }
        } else {
            if (!renderObject) {
                throw new Error('expected renderObject to be available');
            }

            if (updateState.updateColor || updateState.updateSize || updateState.updateTransform) {
                // console.log('update locationIterator');
                locationIt = createLocationIterator(newVolume, newKey);
            }

            if (updateState.updateTransform) {
                // console.log('update transform');
                const { instanceCount, groupCount } = locationIt;
                if (newProps.instanceGranularity) {
                    createMarkers(instanceCount, 'instance', renderObject.values);
                } else {
                    createMarkers(instanceCount * groupCount, 'groupInstance', renderObject.values);
                }
            }

            if (updateState.updateMatrix) {
                // console.log('update matrix');
                createVolumeInstancesTransform(newVolume, geometry.boundingSphere, newProps.cellSize, newProps.batchSize, renderObject.values);
            }

            if (updateState.createGeometry) {
                if (newGeometry) {
                    ValueCell.updateIfChanged(renderObject.values.drawCount, Geometry.getDrawCount(newGeometry));
                    ValueCell.updateIfChanged(renderObject.values.uVertexCount, Geometry.getVertexCount(newGeometry));
                    ValueCell.updateIfChanged(renderObject.values.uGroupCount, locationIt.groupCount);
                } else {
                    throw new Error('expected geometry to be given');
                }
            }

            if (updateState.updateTransform || updateState.createGeometry) {
                updateBoundingSphere(renderObject.values, newGeometry || geometry);
                positionIt = createPositionIterator(newGeometry || geometry, renderObject.values);
            }

            if (updateState.updateSize) {
                // not all geometries have size data, so check here
                if ('uSize' in renderObject.values) {
                    // console.log('update size');
                    createSizes(locationIt, positionIt, newTheme.size, renderObject.values as SizeValues);
                }
            }

            if (updateState.updateColor) {
                // console.log('update color');
                createColors(locationIt, positionIt, newTheme.color, renderObject.values);
            }

            updateValues(renderObject.values, newProps);
            updateRenderableState(renderObject.state, newProps);
        }

        currentProps = newProps;
        currentTheme = newTheme;
        currentVolume = newVolume;
        currentKey = newKey;
        if (newGeometry) {
            geometry = newGeometry;
            geometryVersion += 1;
        }
    }

    function eachInstance(loci: Loci, volume: Volume, key: number, apply: (interval: Interval) => boolean) {
        let changed = false;
        if (Volume.Cell.isLoci(loci)) {
            if (Volume.Cell.isLociEmpty(loci)) return false;
            if (!Volume.areEquivalent(loci.volume, volume)) return false;
            for (const { instances } of loci.elements) {
                OrderedSet.forEach(instances, j => {
                    if (apply(Interval.ofSingleton(j))) changed = true;
                });
            }
        } else if (Volume.Segment.isLoci(loci)) {
            if (Volume.Segment.isLociEmpty(loci)) return false;
            if (!Volume.areEquivalent(loci.volume, volume)) return false;
            for (const { segments, instances } of loci.elements) {
                if (OrderedSet.has(segments, key)) {
                    OrderedSet.forEach(instances, j => {
                        if (apply(Interval.ofSingleton(j))) changed = true;
                    });
                }
            }
        } else if (Volume.Isosurface.isLoci(loci)) {
            if (Volume.Isosurface.isLociEmpty(loci)) return false;
            if (Interval.is(loci.instances)) {
                if (apply(loci.instances)) changed = true;
            } else {
                for (let i = 0, il = loci.instances.length; i < il; ++i) {
                    if (apply(Interval.ofSingleton(i))) changed = true;
                }
            }
        } else if (Volume.isLoci(loci)) {
            if (Volume.isLociEmpty(loci)) return false;
            if (Interval.is(loci.instances)) {
                if (apply(loci.instances)) changed = true;
            } else {
                for (let i = 0, il = loci.instances.length; i < il; ++i) {
                    if (apply(Interval.ofSingleton(i))) changed = true;
                }
            }
        }
        return changed;
    }

    function lociApply(loci: Loci, apply: (interval: Interval) => boolean) {
        if (isEveryLoci(loci)) {
            if (currentProps.instanceGranularity) {
                return apply(Interval.ofBounds(0, locationIt.instanceCount));
            } else {
                return apply(Interval.ofBounds(0, locationIt.groupCount * locationIt.instanceCount));
            }
        } else {
            if (currentProps.instanceGranularity) {
                return eachInstance(loci, currentVolume, currentKey, apply);
            } else {
                return eachLocation(loci, currentVolume, currentKey, currentProps, apply);
            }
        }
    }

    return {
        get groupCount() { return locationIt ? locationIt.count : 0; },
        get renderObject() { return renderObject; },
        get geometryVersion() { return geometryVersion; },
        async createOrUpdate(ctx: VisualContext, theme: Theme, props: Partial<PD.Values<P>> = {}, volumeKey?: VolumeKey) {
            prepareUpdate(theme, props, volumeKey?.volume || currentVolume, volumeKey?.key || currentKey);
            if (updateState.createGeometry) {
                const newGeometry = createGeometry(ctx, newVolume, newKey, newTheme, newProps, geometry);
                return isPromiseLike(newGeometry) ? newGeometry.then(update) : update(newGeometry);
            } else {
                update();
            }
        },
        getLoci(pickingId: PickingId) {
            return renderObject ? getLoci(pickingId, currentVolume, currentKey, currentProps, renderObject.id) : EmptyLoci;
        },
        eachLocation(cb: LocationCallback) {
            locationIt.reset();
            while (locationIt.hasNext) {
                const { location, isSecondary } = locationIt.move();
                cb(location, isSecondary);
            }
        },
        mark(loci: Loci, action: MarkerAction) {
            return Visual.mark(renderObject, loci, action, lociApply);
        },
        setVisibility(visible: boolean) {
            Visual.setVisibility(renderObject, visible);
        },
        setAlphaFactor(alphaFactor: number) {
            Visual.setAlphaFactor(renderObject, alphaFactor);
        },
        setPickable(pickable: boolean) {
            Visual.setPickable(renderObject, pickable);
        },
        setColorOnly(colorOnly: boolean) {
            Visual.setColorOnly(renderObject, colorOnly);
        },
        setTransform(matrix?: Mat4, instanceMatrices?: Float32Array | null) {
            Visual.setTransform(renderObject, matrix, instanceMatrices);
        },
        setOverpaint(overpaint: Overpaint) {
            return Visual.setOverpaint(renderObject, overpaint, lociApply, true);
        },
        setTransparency(transparency: Transparency) {
            return Visual.setTransparency(renderObject, transparency, lociApply, true);
        },
        setEmissive(emissive: Emissive) {
            return Visual.setEmissive(renderObject, emissive, lociApply, true);
        },
        setSubstance(substance: Substance) {
            return Visual.setSubstance(renderObject, substance, lociApply, true);
        },
        setClipping(clipping: Clipping) {
            return Visual.setClipping(renderObject, clipping, lociApply, true);
        },
        setThemeStrength(strength: { overpaint: number, transparency: number, emissive: number, substance: number }) {
            Visual.setThemeStrength(renderObject, strength);
        },
        destroy() {
            dispose?.(geometry);
            if (renderObject) {
                renderObject.state.disposed = true;
                renderObject = undefined;
            }
        },
        mustRecreate
    };
}
