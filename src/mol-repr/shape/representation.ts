/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Geometry, GeometryUtils } from '../../mol-geo/geometry/geometry';
import { Representation } from '../representation';
import { Shape, ShapeGroup } from '../../mol-model/shape';
import { Subject } from 'rxjs';
import { getNextMaterialId, createRenderObject, RenderObjectValues, GraphicsRenderObject } from '../../mol-gl/render-object';
import { Theme } from '../../mol-theme/theme';
import { LocationIterator } from '../../mol-geo/util/location-iterator';
import { VisualUpdateState } from '../util';
import { createMarkers } from '../../mol-geo/geometry/marker-data';
import { MarkerAction, MarkerActions } from '../../mol-util/marker-action';
import { ValueCell } from '../../mol-util';
import { createColors } from '../../mol-geo/geometry/color-data';
import { createSizes, SizeData } from '../../mol-geo/geometry/size-data';
import { Loci, isEveryLoci, EmptyLoci } from '../../mol-model/loci';
import { Interval, OrderedSet } from '../../mol-data/int';
import { PickingId } from '../../mol-geo/geometry/picking';
import { Visual } from '../visual';
import { RuntimeContext, Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';

export interface ShapeRepresentation<D, G extends Geometry, P extends Geometry.Params<G>> extends Representation<D, P> { }

export type ShapeGetter<D, G extends Geometry, P extends Geometry.Params<G>> = (ctx: RuntimeContext, data: D, props: PD.Values<P>, shape?: Shape<G>) => Shape<G> | Promise<Shape<G>>

export interface ShapeBuilder<G extends Geometry, P extends Geometry.Params<G>> {
    /** Hook to modify represetantion props */
    modifyProps?: (props: Partial<PD.Values<P>>) => Partial<PD.Values<P>>
    /** Hook to modify representation state */
    modifyState?: (state: Partial<Representation.State>) => Partial<Representation.State>
}

export function ShapeRepresentation<D, G extends Geometry, P extends Geometry.Params<G>>(getShape: ShapeGetter<D, G, P>, geometryUtils: GeometryUtils<G>, builder: ShapeBuilder<G, P> = {}): ShapeRepresentation<D, G, P> {
    let version = 0;
    const updated = new Subject<number>();
    const _state = Representation.createState();
    const materialId = getNextMaterialId();
    const renderObjects: GraphicsRenderObject<G['kind']>[] = [];
    let _renderObject: GraphicsRenderObject<G['kind']> | undefined;
    let _shape: Shape<G>;
    let _theme = Theme.createEmpty();
    let currentProps: PD.Values<P> = PD.getDefaultValues(geometryUtils.Params as P); // TODO avoid casting
    let currentParams: P;
    let locationIt: LocationIterator;

    if (builder.modifyState) Representation.updateState(_state, builder.modifyState(_state));

    const updateState = VisualUpdateState.create();

    function prepareUpdate(props: Partial<PD.Values<P>> = {}, shape?: Shape<G>) {
        VisualUpdateState.reset(updateState);

        if (!shape && !_shape) {
            // console.error('no shape given')
            return;
        } else if (shape && !_shape) {
            // console.log('first shape')
            updateState.createNew = true;
        } else if (shape && _shape && shape.id === _shape.id) {
            // console.log('same shape')
            // nothing to set
        } else if (shape && _shape && shape.id !== _shape.id) {
            // console.log('new shape')
            updateState.updateTransform = true;
            updateState.createGeometry = true;
        } else if (!shape) {
            // console.log('only props')
            // nothing to set
        } else {
            console.warn('unexpected state');
        }

        if (updateState.updateTransform) {
            updateState.updateColor = true;
            updateState.updateSize = true;
        }

        if (updateState.createGeometry) {
            updateState.updateColor = true;
            updateState.updateSize = true;
        }
    }

    function createOrUpdate(props: Partial<PD.Values<P>> = {}, data?: D) {
        if (builder.modifyProps) props = builder.modifyProps(props);

        return Task.create('ShapeRepresentation.create', async runtime => {
            const newProps = Object.assign(currentProps, props);
            const shape = data ? await getShape(runtime, data, newProps, _shape) : undefined;

            prepareUpdate(props, shape);

            if (shape) {
                _shape = shape;
                Object.assign(_theme, Shape.getTheme(_shape));
            }

            if (updateState.createNew) {
                renderObjects.length = 0; // clear list o renderObjects
                locationIt = Shape.groupIterator(_shape);
                const transform = Shape.createTransform(_shape.transforms);
                const values = geometryUtils.createValues(_shape.geometry, transform, locationIt, _theme, newProps);
                const state = geometryUtils.createRenderableState(newProps);
                if (builder.modifyState) Object.assign(state, builder.modifyState(state));
                Representation.updateState(_state, state);

                _renderObject = createRenderObject(_shape.geometry.kind, values, state, materialId);
                if (_renderObject) renderObjects.push(_renderObject); // add new renderObject to list
            } else {
                if (!_renderObject) {
                    throw new Error('expected renderObject to be available');
                }

                if (updateState.updateTransform) {
                    // console.log('update transform')
                    Shape.createTransform(_shape.transforms, _renderObject.values);
                    locationIt = Shape.groupIterator(_shape);
                    const { instanceCount, groupCount } = locationIt;
                    createMarkers(instanceCount * groupCount, _renderObject.values);
                }

                if (updateState.createGeometry) {
                    // console.log('update geometry')
                    ValueCell.update(_renderObject.values.drawCount, Geometry.getDrawCount(_shape.geometry));
                }

                if (updateState.updateTransform || updateState.createGeometry) {
                    // console.log('updateBoundingSphere')
                    geometryUtils.updateBoundingSphere(_renderObject.values as RenderObjectValues<G['kind']>, _shape.geometry);
                }

                if (updateState.updateColor) {
                    // console.log('update color')
                    createColors(locationIt, _theme.color, _renderObject.values);
                }

                if (updateState.updateSize) {
                    // not all geometries have size data, so check here
                    if ('uSize' in _renderObject.values) {
                        // console.log('update size')
                        createSizes(locationIt, _theme.size, _renderObject.values as SizeData);
                    }
                }

                geometryUtils.updateValues(_renderObject.values as RenderObjectValues<G['kind']>, newProps);
                geometryUtils.updateRenderableState(_renderObject.state, newProps);
            }

            currentProps = newProps;
            // increment version
            updated.next(version++);
        });
    }

    function lociApply(loci: Loci, apply: (interval: Interval) => boolean) {
        if (isEveryLoci(loci) || (Shape.isLoci(loci) && loci.shape === _shape)) {
            return apply(Interval.ofBounds(0, _shape.groupCount * _shape.transforms.length));
        } else {
            return eachShapeGroup(loci, _shape, apply);
        }
    }

    return {
        label: 'Shape geometry',
        get groupCount () { return locationIt ? locationIt.count : 0; },
        get props () { return currentProps; },
        get params () { return currentParams; },
        get state() { return _state; },
        get theme() { return _theme; },
        renderObjects,
        updated,
        createOrUpdate,
        getLoci(pickingId?: PickingId) {
            if (pickingId === undefined) return Shape.Loci(_shape);
            const { objectId, groupId, instanceId } = pickingId;
            if (_renderObject && _renderObject.id === objectId) {
                return ShapeGroup.Loci(_shape, [{ ids: OrderedSet.ofSingleton(groupId), instance: instanceId }]);
            }
            return EmptyLoci;
        },
        mark(loci: Loci, action: MarkerAction) {
            if (!MarkerActions.is(_state.markerActions, action)) return false;
            if (ShapeGroup.isLoci(loci) || Shape.isLoci(loci)) {
                if (loci.shape !== _shape) return false;
            } else if (!isEveryLoci(loci)) {
                return false;
            }
            return Visual.mark(_renderObject, loci, action, lociApply);
        },
        setState(state: Partial<Representation.State>) {
            if (builder.modifyState) state = builder.modifyState(state);

            if (_renderObject) {
                if (state.visible !== undefined) Visual.setVisibility(_renderObject, state.visible);
                if (state.alphaFactor !== undefined) Visual.setAlphaFactor(_renderObject, state.alphaFactor);
                if (state.pickable !== undefined) Visual.setPickable(_renderObject, state.pickable);
                if (state.overpaint !== undefined) {
                    Visual.setOverpaint(_renderObject, state.overpaint, lociApply, true);
                }
                if (state.transparency !== undefined) {
                    Visual.setTransparency(_renderObject, state.transparency, lociApply, true);
                }
                if (state.transform !== undefined) Visual.setTransform(_renderObject, state.transform);
            }

            Representation.updateState(_state, state);
        },
        setTheme(theme: Theme) {
            console.warn('The `ShapeRepresentation` theme is fixed to `ShapeGroupColorTheme` and `ShapeGroupSizeTheme`. Colors are taken from `Shape.getColor` and sizes from `Shape.getSize`');
        },
        destroy() {
            // TODO
            renderObjects.length = 0;
            _renderObject = undefined;
        }
    };
}

function eachShapeGroup(loci: Loci, shape: Shape, apply: (interval: Interval) => boolean) {
    if (!ShapeGroup.isLoci(loci)) return false;
    if (loci.shape !== shape) return false;
    let changed = false;
    const { groupCount } = shape;
    const { groups } = loci;
    for (const { ids, instance } of groups) {
        if (Interval.is(ids)) {
            const start = instance * groupCount + Interval.start(ids);
            const end = instance * groupCount + Interval.end(ids);
            if (apply(Interval.ofBounds(start, end))) changed = true;
        } else {
            for (let i = 0, _i = ids.length; i < _i; i++) {
                const idx = instance * groupCount + ids[i];
                if (apply(Interval.ofSingleton(idx))) changed = true;
            }
        }
    }
    return changed;
}