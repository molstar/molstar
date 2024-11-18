/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Structure, Unit, StructureElement, Bond } from '../../mol-model/structure';
import { RepresentationProps } from '../representation';
import { Visual, VisualContext } from '../visual';
import { Geometry, GeometryUtils } from '../../mol-geo/geometry/geometry';
import { LocationIterator } from '../../mol-geo/util/location-iterator';
import { Theme } from '../../mol-theme/theme';
import { createUnitsTransform, includesUnitKind, StructureGroup } from './visual/util/common';
import { createRenderObject, GraphicsRenderObject, RenderObjectValues } from '../../mol-gl/render-object';
import { PickingId } from '../../mol-geo/geometry/picking';
import { Loci, isEveryLoci, EmptyLoci } from '../../mol-model/loci';
import { Interval } from '../../mol-data/int';
import { LocationCallback, VisualUpdateState } from '../util';
import { ColorTheme } from '../../mol-theme/color';
import { createMarkers } from '../../mol-geo/geometry/marker-data';
import { MarkerAction } from '../../mol-util/marker-action';
import { ValueCell, deepEqual } from '../../mol-util';
import { createSizes } from '../../mol-geo/geometry/size-data';
import { createColors } from '../../mol-geo/geometry/color-data';
import { Mat4 } from '../../mol-math/linear-algebra';
import { Overpaint } from '../../mol-theme/overpaint';
import { Transparency } from '../../mol-theme/transparency';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { SizeTheme } from '../../mol-theme/size';
import { Spheres } from '../../mol-geo/geometry/spheres/spheres';
import { Cylinders } from '../../mol-geo/geometry/cylinders/cylinders';
import { Points } from '../../mol-geo/geometry/points/points';
import { Lines } from '../../mol-geo/geometry/lines/lines';
import { Text } from '../../mol-geo/geometry/text/text';
import { DirectVolume } from '../../mol-geo/geometry/direct-volume/direct-volume';
import { TextureMesh } from '../../mol-geo/geometry/texture-mesh/texture-mesh';
import { SizeValues } from '../../mol-gl/renderable/schema';
import { StructureParams, StructureMeshParams, StructureSpheresParams, StructurePointsParams, StructureLinesParams, StructureTextParams, StructureDirectVolumeParams, StructureTextureMeshParams, StructureCylindersParams } from './params';
import { Clipping } from '../../mol-theme/clipping';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { isPromiseLike } from '../../mol-util/type-helpers';
import { Substance } from '../../mol-theme/substance';
import { Emissive } from '../../mol-theme/emissive';

export interface UnitsVisual<P extends RepresentationProps = {}> extends Visual<StructureGroup, P> { }

function createUnitsRenderObject<G extends Geometry>(structureGroup: StructureGroup, geometry: G, locationIt: LocationIterator, theme: Theme, props: PD.Values<StructureParams & Geometry.Params<G>>, materialId: number) {
    const { createValues, createRenderableState } = Geometry.getUtils(geometry);
    const transform = createUnitsTransform(structureGroup, props.includeParent, geometry.boundingSphere, props.cellSize, props.batchSize);
    const values = createValues(geometry, transform, locationIt, theme, props);
    const state = createRenderableState(props);
    return createRenderObject(geometry.kind, values, state, materialId);
}

interface UnitsVisualBuilder<P extends StructureParams, G extends Geometry> {
    defaultProps: PD.Values<P>
    createGeometry(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<P>, geometry?: G): Promise<G> | G
    createLocationIterator(structureGroup: StructureGroup, props: PD.Values<P>): LocationIterator
    getLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number): Loci
    eachLocation(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean, isMarking: boolean): boolean
    setUpdateState(state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructureGroup: StructureGroup, currentStructureGroup: StructureGroup): void
    initUpdateState?: (state: VisualUpdateState, newProps: PD.Values<P>, newTheme: Theme, newStructureGroup: StructureGroup) => void
    mustRecreate?: (structureGroup: StructureGroup, props: PD.Values<P>) => boolean
    processValues?: (values: RenderObjectValues<G['kind']>, geometry: G, props: PD.Values<P>, theme: Theme, webgl?: WebGLContext) => void
    dispose?: (geometry: G) => void
}

interface UnitsVisualGeometryBuilder<P extends StructureParams, G extends Geometry> extends UnitsVisualBuilder<P, G> {
    geometryUtils: GeometryUtils<G>
}

export function UnitsVisual<G extends Geometry, P extends StructureParams & Geometry.Params<G>>(builder: UnitsVisualGeometryBuilder<P, G>, materialId: number): UnitsVisual<P> {
    const { defaultProps, createGeometry, createLocationIterator, getLoci, eachLocation, setUpdateState, initUpdateState, mustRecreate, processValues, dispose } = builder;
    const { createEmpty: createEmptyGeometry, updateValues, updateBoundingSphere, updateRenderableState, createPositionIterator } = builder.geometryUtils;
    const updateState = VisualUpdateState.create();
    const previousMark: Visual.PreviousMark = { loci: EmptyLoci, action: MarkerAction.None, status: -1 };

    let renderObject: GraphicsRenderObject<G['kind']> | undefined;

    let newProps: PD.Values<P> = Object.assign({}, defaultProps);
    let newTheme: Theme = Theme.createEmpty();
    let newStructureGroup: StructureGroup;

    let currentProps: PD.Values<P>;
    let currentTheme: Theme;
    let currentStructureGroup: StructureGroup;

    let geometry: G;
    let geometryVersion = -1;
    let locationIt: LocationIterator;
    let positionIt: LocationIterator;

    function prepareUpdate(theme: Theme, props: PD.Values<P>, structureGroup: StructureGroup) {
        if (!structureGroup && !currentStructureGroup) {
            throw new Error('missing structureGroup');
        }

        newProps = props;
        newTheme = theme;
        newStructureGroup = structureGroup;

        VisualUpdateState.reset(updateState);

        if (!renderObject || !currentStructureGroup) {
            initUpdateState?.(updateState, newProps, newTheme, newStructureGroup);
            // console.log('create new');
            updateState.createNew = true;
            updateState.createGeometry = true;
            return;
        }

        setUpdateState(updateState, newProps, currentProps, newTheme, currentTheme, newStructureGroup, currentStructureGroup);

        if (!Structure.areHierarchiesEqual(currentStructureGroup.structure, newStructureGroup.structure)) {
            // console.log('new hierarchy');
            updateState.updateTransform = true;
            updateState.updateColor = true;
            updateState.updateSize = true;
        }

        if (!ColorTheme.areEqual(newTheme.color, currentTheme.color)) {
            // console.log('new colorTheme');
            updateState.updateColor = true;
        }

        if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) {
            // console.log('new sizeTheme');
            updateState.updateSize = true;
        }

        if (currentStructureGroup.structure.child !== newStructureGroup.structure.child) {
            // console.log('new child');
            updateState.createGeometry = true;
            updateState.updateTransform = true;
        }

        if (newProps.instanceGranularity !== currentProps.instanceGranularity || newProps.cellSize !== currentProps.cellSize || newProps.batchSize !== currentProps.batchSize) {
            updateState.updateTransform = true;
        }

        if (!deepEqual(newProps.unitKinds, currentProps.unitKinds)) {
            // console.log('new unitKinds');
            updateState.createGeometry = true;
        }

        if (newStructureGroup.group.transformHash !== currentStructureGroup.group.transformHash) {
            // console.log('new transformHash');
            if (newStructureGroup.group.units.length !== currentStructureGroup.group.units.length || updateState.updateColor) {
                updateState.updateTransform = true;
            } else {
                updateState.updateMatrix = true;
            }
        }

        // check if the operator or conformation of unit has changed
        const newUnit = newStructureGroup.group.units[0];
        const currentUnit = currentStructureGroup.group.units[0];
        if (!Unit.areOperatorsEqual(newUnit, currentUnit)) {
            // console.log('new operators');
            updateState.updateTransform = true;
        }
        if (!Unit.areConformationsEqual(newUnit, currentUnit)) {
            // console.log('new conformation');
            updateState.createGeometry = true;
        }

        if (updateState.updateTransform) {
            updateState.updateMatrix = true;
        }

        if (updateState.updateSize && !('uSize' in renderObject.values)) {
            updateState.createGeometry = true;
        }

        if (updateState.createGeometry || updateState.updateTransform) {
            if (currentStructureGroup.structure.hashCode !== newStructureGroup.structure.hashCode) {
                // console.log('new hashCode');
                updateState.updateColor = true;
                updateState.updateSize = true;
            }
            if (newTheme.color.granularity.startsWith('vertex') ||
                renderObject.values.dColorType.ref.value.startsWith('vertex') ||
                newTheme.color.granularity.startsWith('volume') ||
                renderObject.values.dColorType.ref.value.startsWith('volume')
            ) {
                updateState.updateColor = true;
            }
        }
    }

    function update(newGeometry?: G) {
        if (updateState.createNew) {
            locationIt = createLocationIterator(newStructureGroup, newProps);
            if (newGeometry) {
                renderObject = createUnitsRenderObject(newStructureGroup, newGeometry, locationIt, newTheme, newProps, materialId);
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
                locationIt = createLocationIterator(newStructureGroup, newProps);
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
                createUnitsTransform(newStructureGroup, newProps.includeParent, renderObject.values.invariantBoundingSphere.ref.value, newProps.cellSize, newProps.batchSize, renderObject.values);
                if ('lodLevels' in renderObject.values) {
                    // to trigger `uLod` update in `renderable.cull`
                    ValueCell.update(renderObject.values.lodLevels, renderObject.values.lodLevels.ref.value);
                }
            }

            if (updateState.createGeometry) {
                // console.log('update geometry');
                if (newGeometry) {
                    ValueCell.updateIfChanged(renderObject.values.drawCount, Geometry.getDrawCount(newGeometry));
                    ValueCell.updateIfChanged(renderObject.values.uVertexCount, Geometry.getVertexCount(newGeometry));
                    ValueCell.updateIfChanged(renderObject.values.uGroupCount, Geometry.getGroupCount(newGeometry));
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
                    createSizes(locationIt, newTheme.size, renderObject.values as SizeValues);
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
        currentStructureGroup = newStructureGroup;
        if (newGeometry) {
            geometry = newGeometry;
            geometryVersion += 1;
        }
    }

    function _createGeometry(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: PD.Values<P>, geometry?: G) {
        return includesUnitKind(props.unitKinds, unit)
            ? createGeometry(ctx, unit, structure, theme, props, geometry)
            : createEmptyGeometry(geometry);
    }

    function lociIsSuperset(loci: Loci) {
        if (isEveryLoci(loci)) return true;
        if (Structure.isLoci(loci) && Structure.areRootsEquivalent(loci.structure, currentStructureGroup.structure)) return true;
        if (StructureElement.Loci.is(loci) && Structure.areRootsEquivalent(loci.structure, currentStructureGroup.structure)) {
            if (StructureElement.Loci.isWholeStructure(loci)) return true;
        }
        return false;
    }

    function eachInstance(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean) {
        let changed = false;
        if (Bond.isLoci(loci)) {
            const { structure, group } = structureGroup;
            if (!Structure.areEquivalent(loci.structure, structure)) return false;
            for (const b of loci.bonds) {
                if (b.aUnit !== b.bUnit) continue;
                const unitIdx = group.unitIndexMap.get(b.aUnit.id);
                if (unitIdx !== undefined) {
                    if (apply(Interval.ofSingleton(unitIdx))) changed = true;
                }
            }
        } else if (StructureElement.Loci.is(loci)) {
            const { structure, group } = structureGroup;
            if (!Structure.areEquivalent(loci.structure, structure)) return false;
            for (const e of loci.elements) {
                const unitIdx = group.unitIndexMap.get(e.unit.id);
                if (unitIdx !== undefined) {
                    if (apply(Interval.ofSingleton(unitIdx))) changed = true;
                }
            }
        }
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
                return eachInstance(loci, currentStructureGroup, apply);
            } else {
                return eachLocation(loci, currentStructureGroup, apply, isMarking);
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
        createOrUpdate(ctx: VisualContext, theme: Theme, props: PD.Values<P>, structureGroup?: StructureGroup) {
            prepareUpdate(theme, props, structureGroup || currentStructureGroup);
            if (updateState.createGeometry) {
                const newGeometry = _createGeometry(ctx, newStructureGroup.group.units[0], newStructureGroup.structure, newTheme, newProps, geometry);
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
            return renderObject ? getLoci(pickingId, currentStructureGroup, renderObject.id) : EmptyLoci;
        },
        eachLocation(cb: LocationCallback) {
            locationIt.reset();
            while (locationIt.hasNext) {
                const { location, isSecondary } = locationIt.move();
                cb(location, isSecondary);
            }
        },
        mark(loci: Loci, action: MarkerAction) {
            let hasInvariantId = true;
            if (StructureElement.Loci.is(loci)) {
                hasInvariantId = false;
                const { invariantId } = currentStructureGroup.group.units[0];
                for (const e of loci.elements) {
                    if (e.unit.invariantId === invariantId) {
                        hasInvariantId = true;
                        break;
                    }
                }
            }
            return hasInvariantId ? Visual.mark(renderObject, loci, action, lociApply, previousMark) : false;
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

export const UnitsMeshParams = { ...StructureMeshParams, ...StructureParams };
export type UnitsMeshParams = typeof UnitsMeshParams
export interface UnitsMeshVisualBuilder<P extends UnitsMeshParams> extends UnitsVisualBuilder<P, Mesh> { }

export function UnitsMeshVisual<P extends UnitsMeshParams>(builder: UnitsMeshVisualBuilder<P>, materialId: number): UnitsVisual<P> {
    return UnitsVisual<Mesh, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructureGroup: StructureGroup, currentStructureGroup: StructureGroup) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructureGroup, currentStructureGroup);
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true;
        },
        geometryUtils: Mesh.Utils
    }, materialId);
}

// spheres

export const UnitsSpheresParams = { ...StructureSpheresParams, ...StructureParams };
export type UnitsSpheresParams = typeof UnitsSpheresParams
export interface UnitsSpheresVisualBuilder<P extends UnitsSpheresParams> extends UnitsVisualBuilder<P, Spheres> { }

export function UnitsSpheresVisual<P extends UnitsSpheresParams>(builder: UnitsSpheresVisualBuilder<P>, materialId: number): UnitsVisual<P> {
    return UnitsVisual<Spheres, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructureGroup: StructureGroup, currentStructureGroup: StructureGroup) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructureGroup, currentStructureGroup);
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true;
        },
        geometryUtils: Spheres.Utils
    }, materialId);
}

// cylinders

export const UnitsCylindersParams = { ...StructureCylindersParams, ...StructureParams };
export type UnitsCylindersParams = typeof UnitsCylindersParams
export interface UnitsCylindersVisualBuilder<P extends UnitsCylindersParams> extends UnitsVisualBuilder<P, Cylinders> { }

export function UnitsCylindersVisual<P extends UnitsCylindersParams>(builder: UnitsCylindersVisualBuilder<P>, materialId: number): UnitsVisual<P> {
    return UnitsVisual<Cylinders, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructureGroup: StructureGroup, currentStructureGroup: StructureGroup) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructureGroup, currentStructureGroup);
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true;
        },
        geometryUtils: Cylinders.Utils
    }, materialId);
}

// points

export const UnitsPointsParams = { ...StructurePointsParams, ...StructureParams };
export type UnitsPointsParams = typeof UnitsPointsParams
export interface UnitsPointVisualBuilder<P extends UnitsPointsParams> extends UnitsVisualBuilder<P, Points> { }

export function UnitsPointsVisual<P extends UnitsPointsParams>(builder: UnitsPointVisualBuilder<P>, materialId: number): UnitsVisual<P> {
    return UnitsVisual<Points, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructureGroup: StructureGroup, currentStructureGroup: StructureGroup) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructureGroup, currentStructureGroup);
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true;
        },
        geometryUtils: Points.Utils
    }, materialId);
}

// lines

export const UnitsLinesParams = { ...StructureLinesParams, ...StructureParams };
export type UnitsLinesParams = typeof UnitsLinesParams
export interface UnitsLinesVisualBuilder<P extends UnitsLinesParams> extends UnitsVisualBuilder<P, Lines> { }

export function UnitsLinesVisual<P extends UnitsLinesParams>(builder: UnitsLinesVisualBuilder<P>, materialId: number): UnitsVisual<P> {
    return UnitsVisual<Lines, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructureGroup: StructureGroup, currentStructureGroup: StructureGroup) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructureGroup, currentStructureGroup);
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true;
        },
        geometryUtils: Lines.Utils
    }, materialId);
}

// text

export const UnitsTextParams = { ...StructureTextParams, ...StructureParams };
export type UnitsTextParams = typeof UnitsTextParams
export interface UnitsTextVisualBuilder<P extends UnitsTextParams> extends UnitsVisualBuilder<P, Text> { }

export function UnitsTextVisual<P extends UnitsTextParams>(builder: UnitsTextVisualBuilder<P>, materialId: number): UnitsVisual<P> {
    return UnitsVisual<Text, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructureGroup: StructureGroup, currentStructureGroup: StructureGroup) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructureGroup, currentStructureGroup);
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

export const UnitsDirectVolumeParams = { ...StructureDirectVolumeParams, ...StructureParams };
export type UnitsDirectVolumeParams = typeof UnitsDirectVolumeParams
export interface UnitsDirectVolumeVisualBuilder<P extends UnitsDirectVolumeParams> extends UnitsVisualBuilder<P, DirectVolume> { }

export function UnitsDirectVolumeVisual<P extends UnitsDirectVolumeParams>(builder: UnitsDirectVolumeVisualBuilder<P>, materialId: number): UnitsVisual<P> {
    return UnitsVisual<DirectVolume, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructureGroup: StructureGroup, currentStructureGroup: StructureGroup) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructureGroup, currentStructureGroup);
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true;
        },
        geometryUtils: DirectVolume.Utils
    }, materialId);
}

// texture-mesh

export const UnitsTextureMeshParams = { ...StructureTextureMeshParams, ...StructureParams };
export type UnitsTextureMeshParams = typeof UnitsTextureMeshParams
export interface UnitsTextureMeshVisualBuilder<P extends UnitsTextureMeshParams> extends UnitsVisualBuilder<P, TextureMesh> { }

export function UnitsTextureMeshVisual<P extends UnitsTextureMeshParams>(builder: UnitsTextureMeshVisualBuilder<P>, materialId: number): UnitsVisual<P> {
    return UnitsVisual<TextureMesh, P>({
        ...builder,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructureGroup: StructureGroup, currentStructureGroup: StructureGroup) => {
            builder.setUpdateState(state, newProps, currentProps, newTheme, currentTheme, newStructureGroup, currentStructureGroup);
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.createGeometry = true;
        },
        geometryUtils: TextureMesh.Utils
    }, materialId);
}