/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Visual, VisualContext } from '../visual';
import { Bond, Structure, StructureElement } from '../../mol-model/structure';
import { Geometry, GeometryUtils } from '../../mol-geo/geometry/geometry';
import { LocationIterator } from '../../mol-geo/util/location-iterator';
import { Theme } from '../../mol-theme/theme';
import { createIdentityTransform } from '../../mol-geo/geometry/transform-data';
import { createRenderObject, GraphicsRenderObject, RenderObjectValues } from '../../mol-gl/render-object';
import { PickingId } from '../../mol-geo/geometry/picking';
import { Loci, isEveryLoci, EmptyLoci } from '../../mol-model/loci';
import { Interval } from '../../mol-data/int';
import { LocationCallback, VisualUpdateState } from '../util';
import { ColorTheme } from '../../mol-theme/color';
import { ValueCell, deepEqual } from '../../mol-util';
import { createSizes, SizeData } from '../../mol-geo/geometry/size-data';
import { createColors } from '../../mol-geo/geometry/color-data';
import { MarkerAction } from '../../mol-util/marker-action';
import { Mat4 } from '../../mol-math/linear-algebra';
import { Overpaint } from '../../mol-theme/overpaint';
import { Transparency } from '../../mol-theme/transparency';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { Cylinders } from '../../mol-geo/geometry/cylinders/cylinders';
import { Lines } from '../../mol-geo/geometry/lines/lines';
import { Text } from '../../mol-geo/geometry/text/text';
import { SizeTheme } from '../../mol-theme/size';
import { DirectVolume } from '../../mol-geo/geometry/direct-volume/direct-volume';
import { createMarkers } from '../../mol-geo/geometry/marker-data';
import { StructureParams, StructureMeshParams, StructureTextParams, StructureDirectVolumeParams, StructureLinesParams, StructureCylindersParams, StructureTextureMeshParams, StructureSpheresParams, StructurePointsParams, StructureImageParams } from './params';
import { Clipping } from '../../mol-theme/clipping';
import { TextureMesh } from '../../mol-geo/geometry/texture-mesh/texture-mesh';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { isPromiseLike } from '../../mol-util/type-helpers';
import { Substance } from '../../mol-theme/substance';
import { Spheres } from '../../mol-geo/geometry/spheres/spheres';
import { Emissive } from '../../mol-theme/emissive';
import { Points } from '../../mol-geo/geometry/points/points';
import { Image } from '../../mol-geo/geometry/image/image';

export interface ComplexVisual<P extends StructureParams> extends Visual<Structure, P> { }

function createComplexRenderObject<G extends Geometry>(structure: Structure, geometry: G, locationIt: LocationIterator, theme: Theme, props: PD.Values<Geometry.Params<G>>, materialId: number) {
    const { createValues, createRenderableState } = Geometry.getUtils(geometry);
    const transform = createIdentityTransform();
    const values = createValues(geometry, transform, locationIt, theme, props);
    const state = createRenderableState(props);
    return createRenderObject(geometry.kind, values, state, materialId);
}

interface ComplexVisualBuilder<P extends StructureParams, G extends Geometry> {
    defaultProps: PD.Values<P>
    createGeometry(ctx: VisualContext, structure: Structure, theme: Theme, props: PD.Values<P>, geometry?: G): Promise<G> | G
    createLocationIterator(structure: Structure, props: PD.Values<P>): LocationIterator
    getLoci(pickingId: PickingId, structure: Structure, id: number): Loci
    eachLocation(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean, isMarking: boolean): boolean,
    setUpdateState(state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure): void
    mustRecreate?: (structure: Structure, props: PD.Values<P>) => boolean
    processValues?: (values: RenderObjectValues<G['kind']>, geometry: G, props: PD.Values<P>, theme: Theme, webgl?: WebGLContext) => void
    dispose?: (geometry: G) => void
}

interface ComplexVisualGeometryBuilder<P extends StructureParams, G extends Geometry> extends ComplexVisualBuilder<P, G> {
    geometryUtils: GeometryUtils<G>
}

export function ComplexVisual<G extends Geometry, P extends StructureParams & Geometry.Params<G>>(builder: ComplexVisualGeometryBuilder<P, G>, materialId: number): ComplexVisual<P> {
    const { defaultProps, createGeometry, createLocationIterator, getLoci, eachLocation, setUpdateState, mustRecreate, processValues, dispose } = builder;
    const { updateValues, updateBoundingSphere, updateRenderableState, createPositionIterator } = builder.geometryUtils;
    const updateState = VisualUpdateState.create();
    const previousMark: Visual.PreviousMark = { loci: EmptyLoci, action: MarkerAction.None, status: -1 };

    let renderObject: GraphicsRenderObject<G['kind']> | undefined;

    let newProps: PD.Values<P>;
    let newTheme: Theme;
    let newStructure: Structure;

    let currentProps: PD.Values<P> = Object.assign({}, defaultProps);
    let currentTheme: Theme = Theme.createEmpty();
    let currentStructure: Structure;

    let geometry: G;
    let geometryVersion = -1;
    let locationIt: LocationIterator;
    let positionIt: LocationIterator;

    function prepareUpdate(theme: Theme, props: Partial<PD.Values<P>>, structure: Structure) {
        if (!structure && !currentStructure) {
            throw new Error('missing structure');
        }

        newProps = Object.assign({}, currentProps, props);
        newTheme = theme;
        newStructure = structure;

        VisualUpdateState.reset(updateState);

        if (!renderObject || !currentStructure) {
            updateState.createNew = true;
            updateState.createGeometry = true;
            return;
        }

        setUpdateState(updateState, newProps, currentProps, newTheme, currentTheme, newStructure, currentStructure);

        if (!Structure.areEquivalent(newStructure, currentStructure)) {
            updateState.createGeometry = true;
        }

        if (!Structure.areHierarchiesEqual(newStructure, currentStructure)) {
            updateState.updateTransform = true;
            updateState.createGeometry = true;
        }

        if (!ColorTheme.areEqual(theme.color, currentTheme.color)) {
            updateState.updateColor = true;
        }

        if (!SizeTheme.areEqual(theme.size, currentTheme.size)) {
            updateState.updateSize = true;
        }

        if (!deepEqual(newProps.unitKinds, currentProps.unitKinds)) {
            updateState.createGeometry = true;
        }

        if (currentStructure.child !== newStructure.child) {
            // console.log('new child');
            updateState.createGeometry = true;
            updateState.updateTransform = true;
        }

        if (newProps.instanceGranularity !== currentProps.instanceGranularity) {
            updateState.updateTransform = true;
        }

        if (updateState.updateSize && !('uSize' in renderObject.values)) {
            updateState.createGeometry = true;
        }

        if (updateState.createGeometry) {
            updateState.updateColor = true;
            updateState.updateSize = true;
        }
    }

    function update(newGeometry?: G) {
        if (updateState.createNew) {
            locationIt = createLocationIterator(newStructure, newProps);
            if (newGeometry) {
                renderObject = createComplexRenderObject(newStructure, newGeometry, locationIt, newTheme, newProps, materialId);
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
                locationIt = createLocationIterator(newStructure, newProps);
            }

            if (updateState.updateTransform) {
                // console.log('update transform')
                const { instanceCount, groupCount } = locationIt;
                if (newProps.instanceGranularity) {
                    createMarkers(instanceCount, 'instance', renderObject.values);
                } else {
                    createMarkers(instanceCount * groupCount, 'groupInstance', renderObject.values);
                }
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
                positionIt = createPositionIterator(geometry, renderObject.values);
            }

            if (updateState.updateSize) {
                // not all geometries have size data, so check here
                if ('uSize' in renderObject.values) {
                    createSizes(locationIt, newTheme.size, renderObject.values as SizeData);
                }
            }

            if (updateState.updateColor) {
                createColors(locationIt, positionIt, newTheme.color, renderObject.values);
            }

            updateValues(renderObject.values, newProps);
            updateRenderableState(renderObject.state, newProps);
        }

        currentProps = newProps;
        currentTheme = newTheme;
        currentStructure = newStructure;
        if (newGeometry) {
            geometry = newGeometry;
            geometryVersion += 1;
        }
    }

    function lociIsSuperset(loci: Loci) {
        if (isEveryLoci(loci)) return true;
        if (Structure.isLoci(loci) && Structure.areRootsEquivalent(loci.structure, currentStructure)) return true;
        if (StructureElement.Loci.is(loci) && Structure.areRootsEquivalent(loci.structure, currentStructure)) {
            if (StructureElement.Loci.isWholeStructure(loci)) return true;
        }
        return false;
    }

    function eachInstance(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean) {
        let changed = false;
        if (!StructureElement.Loci.is(loci) && !Bond.isLoci(loci)) return false;
        if (!Structure.areEquivalent(loci.structure, structure)) return false;
        if (apply(Interval.ofSingleton(0))) changed = true;
        return changed;
    }

    function lociApply(loci: Loci, apply: (interval: Interval) => boolean, isMarking: boolean) {
        if (lociIsSuperset(loci)) {
            if (currentProps.instanceGranularity) {
                return apply(Interval.ofBounds(0, locationIt.instanceCount));
            } else {
                return apply(Interval.ofBounds(0, locationIt.groupCount * locationIt.instanceCount));
            }
        } else {
            if (currentProps.instanceGranularity) {
                return eachInstance(loci, currentStructure, apply);
            } else {
                return eachLocation(loci, currentStructure, apply, isMarking);
            }
        }
    }

    function finalize(ctx: VisualContext) {
        if (renderObject) {
            processValues?.(renderObject.values, geometry, currentProps, currentTheme, ctx.webgl);
        }
    }

    return {
        get groupCount() { return locationIt ? locationIt.count : 0; },
        get renderObject() { return locationIt && locationIt.count ? renderObject : undefined; },
        get geometryVersion() { return geometryVersion; },
        createOrUpdate(ctx: VisualContext, theme: Theme, props: Partial<PD.Values<P>> = {}, structure?: Structure) {
            prepareUpdate(theme, props, structure || currentStructure);
            if (updateState.createGeometry) {
                const newGeometry = createGeometry(ctx, newStructure, newTheme, newProps, geometry);
                if (isPromiseLike(newGeometry)) {
                    return newGeometry.then(g => {
                        update(g);
                        finalize(ctx);
                    });
                }
                update(newGeometry);
            } else {
                update();
            }
            finalize(ctx);
        },
        getLoci(pickingId: PickingId) {
            return renderObject ? getLoci(pickingId, currentStructure, renderObject.id) : EmptyLoci;
        },
        eachLocation(cb: LocationCallback) {
            locationIt.reset();
            while (locationIt.hasNext) {
                const { location, isSecondary } = locationIt.move();
                cb(location, isSecondary);
            }
        },
        mark(loci: Loci, action: MarkerAction) {
            return Visual.mark(renderObject, loci, action, lociApply, previousMark);
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
        setOverpaint(overpaint: Overpaint, webgl?: WebGLContext) {
            const smoothing = { geometry, props: currentProps, webgl };
            Visual.setOverpaint(renderObject, overpaint, lociApply, true, smoothing);
        },
        setTransparency(transparency: Transparency, webgl?: WebGLContext) {
            const smoothing = { geometry, props: currentProps, webgl };
            Visual.setTransparency(renderObject, transparency, lociApply, true, smoothing);
        },
        setEmissive(emissive: Emissive, webgl?: WebGLContext) {
            const smoothing = { geometry, props: currentProps, webgl };
            Visual.setEmissive(renderObject, emissive, lociApply, true, smoothing);
        },
        setSubstance(substance: Substance, webgl?: WebGLContext) {
            const smoothing = { geometry, props: currentProps, webgl };
            Visual.setSubstance(renderObject, substance, lociApply, true, smoothing);
        },
        setClipping(clipping: Clipping) {
            Visual.setClipping(renderObject, clipping, lociApply, true);
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

// mesh

export const ComplexMeshParams = { ...StructureMeshParams, ...StructureParams };
export type ComplexMeshParams = typeof ComplexMeshParams

export interface ComplexMeshVisualBuilder<P extends ComplexMeshParams> extends ComplexVisualBuilder<P, Mesh> { }

export function ComplexMeshVisual<P extends ComplexMeshParams>(builder: ComplexMeshVisualBuilder<P>, materialId: number): ComplexVisual<P> {
    return ComplexVisual<Mesh, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructure, currentStructure);
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true;
        },
        geometryUtils: Mesh.Utils
    }, materialId);
}

// spheres

export const ComplexSpheresParams = { ...StructureSpheresParams, ...StructureParams };
export type ComplexSpheresParams = typeof ComplexSpheresParams

export interface ComplexSpheresVisualBuilder<P extends ComplexSpheresParams> extends ComplexVisualBuilder<P, Spheres> { }

export function ComplexSpheresVisual<P extends ComplexSpheresParams>(builder: ComplexSpheresVisualBuilder<P>, materialId: number): ComplexVisual<P> {
    return ComplexVisual<Spheres, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructure, currentStructure);
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true;
        },
        geometryUtils: Spheres.Utils
    }, materialId);
}

// cylinders

export const ComplexCylindersParams = { ...StructureCylindersParams, ...StructureParams };
export type ComplexCylindersParams = typeof ComplexCylindersParams

export interface ComplexCylindersVisualBuilder<P extends ComplexCylindersParams> extends ComplexVisualBuilder<P, Cylinders> { }

export function ComplexCylindersVisual<P extends ComplexCylindersParams>(builder: ComplexCylindersVisualBuilder<P>, materialId: number): ComplexVisual<P> {
    return ComplexVisual<Cylinders, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructure, currentStructure);
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true;
        },
        geometryUtils: Cylinders.Utils
    }, materialId);
}

// points

export const ComplexPointsParams = { ...StructurePointsParams, ...StructureParams };
export type ComplexPointsParams = typeof ComplexPointsParams

export interface ComplexPointsVisualBuilder<P extends ComplexPointsParams> extends ComplexVisualBuilder<P, Points> { }

export function ComplexPointsVisual<P extends ComplexPointsParams>(builder: ComplexPointsVisualBuilder<P>, materialId: number): ComplexVisual<P> {
    return ComplexVisual<Points, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructure, currentStructure);
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true;
        },
        geometryUtils: Points.Utils
    }, materialId);
}

// lines

export const ComplexLinesParams = { ...StructureLinesParams, ...StructureParams };
export type ComplexLinesParams = typeof ComplexLinesParams

export interface ComplexLinesVisualBuilder<P extends ComplexLinesParams> extends ComplexVisualBuilder<P, Lines> { }

export function ComplexLinesVisual<P extends ComplexLinesParams>(builder: ComplexLinesVisualBuilder<P>, materialId: number): ComplexVisual<P> {
    return ComplexVisual<Lines, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructure, currentStructure);
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true;
        },
        geometryUtils: Lines.Utils
    }, materialId);
}

// text

export const ComplexTextParams = { ...StructureTextParams, ...StructureParams };
export type ComplexTextParams = typeof ComplexTextParams

export interface ComplexTextVisualBuilder<P extends ComplexTextParams> extends ComplexVisualBuilder<P, Text> { }

export function ComplexTextVisual<P extends ComplexTextParams>(builder: ComplexTextVisualBuilder<P>, materialId: number): ComplexVisual<P> {
    return ComplexVisual<Text, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructure, currentStructure);
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true;
            if (newProps.background !== currentProps.background) state.createGeometry = true;
            if (newProps.backgroundMargin !== currentProps.backgroundMargin) state.createGeometry = true;
            if (newProps.tether !== currentProps.tether) state.createGeometry = true;
            if (newProps.tetherLength !== currentProps.tetherLength) state.createGeometry = true;
            if (newProps.tetherBaseWidth !== currentProps.tetherBaseWidth) state.createGeometry = true;
            if (newProps.attachment !== currentProps.attachment) state.createGeometry = true;

            if (newProps.fontFamily !== currentProps.fontFamily) state.createGeometry = true;
            if (newProps.fontQuality !== currentProps.fontQuality) state.createGeometry = true;
            if (newProps.fontStyle !== currentProps.fontStyle) state.createGeometry = true;
            if (newProps.fontVariant !== currentProps.fontVariant) state.createGeometry = true;
            if (newProps.fontWeight !== currentProps.fontWeight) state.createGeometry = true;
        },
        geometryUtils: Text.Utils
    }, materialId);
}

// direct-volume

export const ComplexDirectVolumeParams = { ...StructureDirectVolumeParams, ...StructureParams };
export type ComplexDirectVolumeParams = typeof ComplexDirectVolumeParams

export interface ComplexDirectVolumeVisualBuilder<P extends ComplexDirectVolumeParams> extends ComplexVisualBuilder<P, DirectVolume> { }

export function ComplexDirectVolumeVisual<P extends ComplexDirectVolumeParams>(builder: ComplexDirectVolumeVisualBuilder<P>, materialId: number): ComplexVisual<P> {
    return ComplexVisual<DirectVolume, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructure, currentStructure);
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true;
        },
        geometryUtils: DirectVolume.Utils
    }, materialId);
}

// texture-mesh

export const ComplexTextureMeshParams = { ...StructureTextureMeshParams, ...StructureParams };
export type ComplexTextureMeshParams = typeof ComplexTextureMeshParams

export interface ComplexTextureMeshVisualBuilder<P extends ComplexTextureMeshParams> extends ComplexVisualBuilder<P, TextureMesh> { }

export function ComplexTextureMeshVisual<P extends ComplexTextureMeshParams>(builder: ComplexTextureMeshVisualBuilder<P>, materialId: number): ComplexVisual<P> {
    return ComplexVisual<TextureMesh, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructure, currentStructure);
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true;
        },
        geometryUtils: TextureMesh.Utils
    }, materialId);
}

// image

export const ComplexImageParams = { ...StructureImageParams, ...StructureParams };
export type ComplexImageParams = typeof ComplexImageParams

export interface ComplexImageVisualBuilder<P extends ComplexImageParams> extends ComplexVisualBuilder<P, Image> { }

export function ComplexImageVisual<P extends ComplexImageParams>(builder: ComplexImageVisualBuilder<P>, materialId: number): ComplexVisual<P> {
    return ComplexVisual<Image, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructure, currentStructure);
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true;
        },
        geometryUtils: Image.Utils
    }, materialId);
}
