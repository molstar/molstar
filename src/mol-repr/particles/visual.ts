/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Visual, VisualContext } from '../visual';
import { ParticleList, getParticleTransforms, Particle } from '../../mol-model/particles/particle-list';
import { Geometry, GeometryUtils } from '../../mol-geo/geometry/geometry';
import { LocationIterator } from '../../mol-geo/util/location-iterator';
import { Theme } from '../../mol-theme/theme';
import { createTransform, createIdentityTransform, TransformData } from '../../mol-geo/geometry/transform-data';
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
import { Wiggle } from '../../mol-theme/wiggle';
import { SizeTheme } from '../../mol-theme/size';
import { Sphere3D } from '../../mol-math/geometry/primitives/sphere3d';
import { BaseGeometry, resolveInstanceGranularity } from '../../mol-geo/geometry/base';
import { WebGLContext } from '../../mol-gl/webgl/context';

export const ParticleParams = {
    ...BaseGeometry.Params,
};
export type ParticleParams = typeof ParticleParams;

export type ParticleKey = { particles: ParticleList }
export interface ParticleVisual<P extends ParticleParams> extends Visual<ParticleKey, P> { }

export function createParticleTransform(particles: ParticleList, invariantBoundingSphere: Sphere3D, cellSize: number, batchSize: number, transformData?: TransformData) {
    const transformArray = getParticleTransforms(particles);
    return createTransform(transformArray, particles.count, invariantBoundingSphere, cellSize, batchSize, transformData);
}

function createParticleRenderObject<G extends Geometry>(particles: ParticleList, geometry: G, locationIt: LocationIterator, theme: Theme, props: PD.Values<Geometry.Params<G>>, materialId: number) {
    const { createValues, createRenderableState } = Geometry.getUtils(geometry);
    const transform = locationIt.nonInstanceable
        ? createIdentityTransform()
        : createParticleTransform(particles, geometry.boundingSphere, props.cellSize, props.batchSize);
    const values = createValues(geometry, transform, locationIt, theme, props);
    const state = createRenderableState(props);
    ValueCell.update(values.boundingSphere, Particle.getBoundary(particles).sphere);
    return createRenderObject(geometry.kind, values, state, materialId);
}

function eachParticleLoci(loci: Loci, particles: ParticleList, apply: (interval: Interval) => boolean): boolean {
    if (!Particle.isLoci(loci)) return false;
    if (Particle.isLociEmpty(loci)) return false;
    if (loci.particles !== particles) return false;
    let changed = false;
    OrderedSet.forEach(loci.indices, idx => {
        if (apply(Interval.ofSingleton(idx))) changed = true;
    });
    return changed;
}

interface ParticleVisualBuilder<P extends ParticleParams, G extends Geometry> {
    defaultProps: PD.Values<P>
    createGeometry(ctx: VisualContext, particles: ParticleList, theme: Theme, props: PD.Values<P>, geometry?: G): Promise<G> | G
    createLocationIterator(particles: ParticleList, geometry: G): LocationIterator
    getLoci(pickingId: PickingId, particles: ParticleList, props: PD.Values<P>, id: number, geometry: G): Loci
    eachLocation(loci: Loci, particles: ParticleList, props: PD.Values<P>, apply: (interval: Interval) => boolean, geometry: G): boolean
    /**
     * Note: a new `ParticleList` reference does not automatically trigger `createGeometry`.
     * If `createGeometry` builds geometry from per-particle data directly (e.g. baked
     * positions), set `state.createGeometry = true` here when `newParticles !== currentParticles`.
     * Otherwise (e.g. geometry is an instanced template), leave it to the base
     * implementation to only refresh the transform/color/size.
     */
    setUpdateState(state: VisualUpdateState, newParticles: ParticleList, currentParticles: ParticleList, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme): void
    /** Optional hook to override the theme before geometry creation and color/size updates. */
    overrideTheme?: (theme: Theme, props: PD.Values<P>) => Theme
    mustRecreate?: (particleKey: ParticleKey, props: PD.Values<P>, webgl?: WebGLContext) => boolean
    dispose?: (geometry: G) => void
}

interface ParticleVisualGeometryBuilder<P extends ParticleParams, G extends Geometry> extends ParticleVisualBuilder<P, G> {
    geometryUtils: GeometryUtils<G>
}

export function ParticleVisual<G extends Geometry, P extends ParticleParams & Geometry.Params<G>>(builder: ParticleVisualGeometryBuilder<P, G>, materialId: number): ParticleVisual<P> {
    const { defaultProps, createGeometry, createLocationIterator, getLoci, eachLocation, setUpdateState, overrideTheme, mustRecreate, dispose } = builder;
    const { updateValues, updateBoundingSphere, updateRenderableState, createPositionIterator } = builder.geometryUtils;
    const updateState = VisualUpdateState.create();

    let renderObject: GraphicsRenderObject<G['kind']> | undefined;

    let newProps: PD.Values<P>;
    let newTheme: Theme;
    let newParticles: ParticleList;

    let currentProps: PD.Values<P> = Object.assign({}, defaultProps);
    let currentTheme: Theme = Theme.createEmpty();
    let currentParticles: ParticleList;

    let geometry: G;
    let geometryVersion = -1;
    let locationIt: LocationIterator;
    let positionIt: LocationIterator;

    function prepareUpdate(theme: Theme, props: Partial<PD.Values<P>>, particles: ParticleList) {
        if (!particles && !currentParticles) {
            throw new Error('missing particles');
        }

        newProps = Object.assign({}, currentProps, props);
        newTheme = theme;
        newParticles = particles;

        VisualUpdateState.reset(updateState);

        if (!renderObject) {
            updateState.createNew = true;
            updateState.createGeometry = true;
            return;
        }

        setUpdateState(updateState, newParticles, currentParticles, newProps, currentProps, newTheme, currentTheme);

        if (newParticles !== currentParticles) {
            // A new `ParticleList` does not necessarily mean the geometry must be
            // rebuilt from scratch (e.g. a new trajectory frame with the same topology
            // but different coordinates). Refresh transform/color/size unconditionally,
            // but leave the decision of whether the geometry itself needs to be
            // recreated (e.g. because it bakes per-particle positions, like fiber
            // curves) to `setUpdateState` above - mirroring how `ComplexVisual` only
            // recreates geometry when the structure is not equivalent or its
            // hierarchy changed, instead of on every new `Structure` reference.
            updateState.updateTransform = true;
            updateState.updateColor = true;
            updateState.updateSize = true;
        }

        if (!ColorTheme.areEqual(newTheme.color, currentTheme.color)) {
            updateState.updateColor = true;
        }

        if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) {
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
        const effectiveTheme = overrideTheme ? overrideTheme(newTheme, newProps) : newTheme;

        if (updateState.createNew) {
            if (newGeometry) {
                locationIt = createLocationIterator(newParticles, newGeometry);
                renderObject = createParticleRenderObject(newParticles, newGeometry, locationIt, effectiveTheme, newProps, materialId);
                positionIt = createPositionIterator(newGeometry, renderObject.values);
            } else {
                throw new Error('expected geometry to be given');
            }
        } else {
            if (!renderObject) {
                throw new Error('expected renderObject to be available');
            }

            if (updateState.updateColor || updateState.updateSize || updateState.updateTransform || updateState.updateLocation) {
                locationIt = createLocationIterator(newParticles, newGeometry || geometry);
            }

            if (updateState.updateTransform || updateState.updateLocation) {
                const { instanceCount, groupCount } = locationIt;
                if (resolveInstanceGranularity(newProps.instanceGranularity, groupCount, instanceCount)) {
                    createMarkers(instanceCount, 'instance', renderObject.values);
                } else {
                    createMarkers(instanceCount * groupCount, 'groupInstance', renderObject.values);
                }
            }

            if (updateState.updateMatrix) {
                createParticleTransform(newParticles, geometry.boundingSphere, newProps.cellSize, newProps.batchSize, renderObject.values);
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
                    createSizes(locationIt, positionIt, effectiveTheme.size, renderObject.values as SizeValues);
                }
            }

            if (updateState.updateColor) {
                createColors(locationIt, positionIt, effectiveTheme.color, renderObject.values);
            }

            updateValues(renderObject.values, newProps);
            updateRenderableState(renderObject.state, newProps);
        }

        currentProps = newProps;
        currentTheme = newTheme;
        currentParticles = newParticles;
        if (newGeometry) {
            geometry = newGeometry;
            geometryVersion += 1;
        }
    }

    function lociApply(loci: Loci, apply: (interval: Interval) => boolean) {
        const instanceGranularity = resolveInstanceGranularity(currentProps.instanceGranularity, locationIt.groupCount, locationIt.instanceCount);
        if (isEveryLoci(loci)) {
            if (instanceGranularity) {
                return apply(Interval.ofBounds(0, locationIt.instanceCount));
            } else {
                return apply(Interval.ofBounds(0, locationIt.groupCount * locationIt.instanceCount));
            }
        } else {
            if (instanceGranularity) {
                return eachParticleLoci(loci, currentParticles, apply);
            } else {
                return eachLocation(loci, currentParticles, currentProps, apply, geometry);
            }
        }
    }

    return {
        get groupCount() { return locationIt ? locationIt.count : 0; },
        get renderObject() { return renderObject; },
        get geometryVersion() { return geometryVersion; },
        async createOrUpdate(ctx: VisualContext, theme: Theme, props: Partial<PD.Values<P>> = {}, particleKey?: ParticleKey) {
            prepareUpdate(theme, props, particleKey?.particles || currentParticles);
            if (updateState.createGeometry) {
                const newGeometry = createGeometry(ctx, newParticles, newTheme, newProps, geometry);
                return isPromiseLike(newGeometry) ? newGeometry.then(update) : update(newGeometry);
            } else {
                update();
            }
        },
        getLoci(pickingId: PickingId) {
            if (!renderObject) return EmptyLoci;
            if (resolveInstanceGranularity(currentProps.instanceGranularity, locationIt.groupCount, locationIt.instanceCount)) {
                pickingId = { ...pickingId, groupId: PickingId.Null };
            }
            return getLoci(pickingId, currentParticles, currentProps, renderObject.id, geometry);
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
        setWiggle(wiggle: Wiggle) {
            return Visual.setWiggle(renderObject, wiggle, lociApply, true);
        },
        setThemeStrength(strength: { overpaint: number, transparency: number, emissive: number, substance: number, wiggle: number }) {
            Visual.setThemeStrength(renderObject, strength);
        },
        destroy() {
            dispose?.(geometry);
            if (renderObject) {
                renderObject.state.disposed = true;
                renderObject = undefined;
            }
        },
        mustRecreate,
    };
}
