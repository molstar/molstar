/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Visual, VisualContext } from '../visual';
import { Structure, StructureElement } from '../../mol-model/structure';
import { Geometry, GeometryUtils } from '../../mol-geo/geometry/geometry';
import { LocationIterator } from '../../mol-geo/util/location-iterator';
import { Theme } from '../../mol-theme/theme';
import { createIdentityTransform } from '../../mol-geo/geometry/transform-data';
import { createRenderObject, RenderObjectValues, GraphicsRenderObject } from '../../mol-gl/render-object';
import { PickingId } from '../../mol-geo/geometry/picking';
import { Loci, isEveryLoci, EmptyLoci } from '../../mol-model/loci';
import { Interval } from '../../mol-data/int';
import { VisualUpdateState } from '../util';
import { ColorTheme } from '../../mol-theme/color';
import { ValueCell, deepEqual } from '../../mol-util';
import { createSizes, SizeData } from '../../mol-geo/geometry/size-data';
import { createColors } from '../../mol-geo/geometry/color-data';
import { MarkerAction } from '../../mol-util/marker-action';
import { Mat4 } from '../../mol-math/linear-algebra';
import { Overpaint } from '../../mol-theme/overpaint';
import { Transparency } from '../../mol-theme/transparency';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { Text } from '../../mol-geo/geometry/text/text';
import { SizeTheme } from '../../mol-theme/size';
import { DirectVolume } from '../../mol-geo/geometry/direct-volume/direct-volume';
import { createMarkers } from '../../mol-geo/geometry/marker-data';
import { StructureParams, StructureMeshParams, StructureTextParams, StructureDirectVolumeParams } from './params';
import { Clipping } from '../../mol-theme/clipping';

export interface  ComplexVisual<P extends StructureParams> extends Visual<Structure, P> { }

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
    createLocationIterator(structure: Structure): LocationIterator
    getLoci(pickingId: PickingId, structure: Structure, id: number): Loci
    eachLocation(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean, isMarking: boolean): boolean,
    setUpdateState(state: VisualUpdateState, newProps: PD.Values<P>, currentProps: PD.Values<P>, newTheme: Theme, currentTheme: Theme, newStructure: Structure, currentStructure: Structure): void
}

interface ComplexVisualGeometryBuilder<P extends StructureParams, G extends Geometry> extends ComplexVisualBuilder<P, G> {
    geometryUtils: GeometryUtils<G>
}

export function ComplexVisual<G extends Geometry, P extends StructureParams & Geometry.Params<G>>(builder: ComplexVisualGeometryBuilder<P, G>, materialId: number): ComplexVisual<P> {
    const { defaultProps, createGeometry, createLocationIterator, getLoci, eachLocation, setUpdateState } = builder;
    const { updateValues, updateBoundingSphere, updateRenderableState } = builder.geometryUtils;
    const updateState = VisualUpdateState.create();

    let renderObject: GraphicsRenderObject<G['kind']> | undefined;

    let newProps: PD.Values<P>;
    let newTheme: Theme;
    let newStructure: Structure;

    let currentProps: PD.Values<P> = Object.assign({}, defaultProps);
    let currentTheme: Theme = Theme.createEmpty();
    let currentStructure: Structure;

    let geometry: G;
    let locationIt: LocationIterator;

    function prepareUpdate(theme: Theme, props: Partial<PD.Values<P>>, structure: Structure) {
        if (!structure && !currentStructure) {
            throw new Error('missing structure');
        }

        newProps = Object.assign({}, currentProps, props);
        newTheme = theme;
        newStructure = structure;

        VisualUpdateState.reset(updateState);

        if (!renderObject) {
            updateState.createNew = true;
        } else if (!currentStructure || !Structure.areEquivalent(newStructure, currentStructure)) {
            updateState.createNew = true;
        }

        if (updateState.createNew) {
            updateState.createGeometry = true;
            return;
        }

        setUpdateState(updateState, newProps, currentProps, newTheme, currentTheme, newStructure, currentStructure);

        if (Structure.conformationHash(newStructure) !== Structure.conformationHash(currentStructure)) {
            updateState.updateTransform = true;
            updateState.createGeometry = true;
        }

        if (!ColorTheme.areEqual(theme.color, currentTheme.color)) {
            updateState.updateColor = true;
        }

        if (!deepEqual(newProps.unitKinds, currentProps.unitKinds)) {
            updateState.createGeometry = true;
        }

        if (updateState.createGeometry) {
            updateState.updateColor = true;
        }
    }

    function update(newGeometry?: G) {
        if (updateState.createNew) {
            locationIt = createLocationIterator(newStructure);
            if (newGeometry) {
                renderObject = createComplexRenderObject(newStructure, newGeometry, locationIt, newTheme, newProps, materialId);
            } else {
                throw new Error('expected geometry to be given');
            }
        } else {
            if (!renderObject) {
                throw new Error('expected renderObject to be available');
            }

            if (updateState.updateTransform) {
                // console.log('update transform')
                locationIt = createLocationIterator(newStructure);
                const { instanceCount, groupCount } = locationIt;
                createMarkers(instanceCount * groupCount, renderObject.values);
            }

            if (updateState.createGeometry) {
                if (newGeometry) {
                    ValueCell.update(renderObject.values.drawCount, Geometry.getDrawCount(newGeometry));
                    updateBoundingSphere(renderObject.values as RenderObjectValues<G['kind']>, newGeometry);
                } else {
                    throw new Error('expected geometry to be given');
                }
            }

            if (updateState.updateSize) {
                // not all geometries have size data, so check here
                if ('uSize' in renderObject.values) {
                    createSizes(locationIt, newTheme.size, renderObject.values as SizeData);
                }
            }

            if (updateState.updateColor) {
                createColors(locationIt, newTheme.color, renderObject.values);
            }

            updateValues(renderObject.values as RenderObjectValues<G['kind']>, newProps);
            updateRenderableState(renderObject.state, newProps);
        }

        currentProps = newProps;
        currentTheme = newTheme;
        currentStructure = newStructure;
        if (newGeometry) geometry = newGeometry;
    }

    function lociIsSuperset(loci: Loci) {
        if (isEveryLoci(loci)) return true;
        if (Structure.isLoci(loci) && Structure.areRootsEquivalent(loci.structure, currentStructure)) return true;
        if (StructureElement.Loci.is(loci) && Structure.areRootsEquivalent(loci.structure, currentStructure)) {
            if (StructureElement.Loci.isWholeStructure(loci)) return true;
        }
        return false;
    }

    function lociApply(loci: Loci, apply: (interval: Interval) => boolean, isMarking: boolean) {
        if (lociIsSuperset(loci)) {
            return apply(Interval.ofBounds(0, locationIt.groupCount * locationIt.instanceCount));
        } else {
            return eachLocation(loci, currentStructure, apply, isMarking);
        }
    }

    return {
        get groupCount() { return locationIt ? locationIt.count : 0; },
        get renderObject () { return locationIt && locationIt.count ? renderObject : undefined; },
        createOrUpdate(ctx: VisualContext, theme: Theme, props: Partial<PD.Values<P>> = {}, structure?: Structure) {
            prepareUpdate(theme, props, structure || currentStructure);
            if (updateState.createGeometry) {
                const newGeometry = createGeometry(ctx, newStructure, newTheme, newProps, geometry);
                return newGeometry instanceof Promise ? newGeometry.then(update) : update(newGeometry);
            } else {
                update();
            }
        },
        getLoci(pickingId: PickingId) {
            return renderObject ? getLoci(pickingId, currentStructure, renderObject.id) : EmptyLoci;
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
        setTransform(matrix?: Mat4, instanceMatrices?: Float32Array | null) {
            Visual.setTransform(renderObject, matrix, instanceMatrices);
        },
        setOverpaint(overpaint: Overpaint) {
            Visual.setOverpaint(renderObject, overpaint, lociApply, true);
        },
        setTransparency(transparency: Transparency) {
            Visual.setTransparency(renderObject, transparency, lociApply, true);
        },
        setClipping(clipping: Clipping) {
            Visual.setClipping(renderObject, clipping, lociApply, true);
        },
        destroy() {
            // TODO
            renderObject = undefined;
        }
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
            if (!SizeTheme.areEqual(newTheme.size, currentTheme.size)) state.updateSize = true;
        },
        geometryUtils: Mesh.Utils
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