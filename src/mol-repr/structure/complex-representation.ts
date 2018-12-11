/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from 'mol-model/structure';
import { Task } from 'mol-task'
import { Loci, EmptyLoci } from 'mol-model/loci';
import { StructureRepresentation, StructureParams } from './representation';
import { ComplexVisual } from './complex-visual';
import { PickingId } from 'mol-geo/geometry/picking';
import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { RepresentationContext, RepresentationParamsGetter, Representation } from 'mol-repr/representation';
import { Theme, createEmptyTheme } from 'mol-theme/theme';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Subject } from 'rxjs';

export function ComplexRepresentation<P extends StructureParams>(label: string, ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, P>, visualCtor: () => ComplexVisual<P>): StructureRepresentation<P> {
    let version = 0
    const updated = new Subject<number>()
    const _state = Representation.createState()
    let visual: ComplexVisual<P> | undefined

    let _structure: Structure
    let _params: P
    let _props: PD.Values<P>
    let _theme = createEmptyTheme()

    function createOrUpdate(props: Partial<PD.Values<P>> = {}, structure?: Structure) {
        if (structure && structure !== _structure) {
            _params = getParams(ctx, structure)
            _structure = structure
            if (!_props) _props = PD.getDefaultValues(_params)
        }
        _props = Object.assign({}, _props, props)

        return Task.create('Creating or updating ComplexRepresentation', async runtime => {
            if (!visual) visual = visualCtor()
            const promise = visual.createOrUpdate({ webgl: ctx.webgl, runtime }, _theme, _props, structure)
            if (promise) await promise
            updated.next(version++)
        });
    }

    function getLoci(pickingId: PickingId) {
        return visual ? visual.getLoci(pickingId) : EmptyLoci
    }

    function mark(loci: Loci, action: MarkerAction) {
        return visual ? visual.mark(loci, action) : false
    }

    function setState(state: Partial<Representation.State>) {
        if (state.visible !== undefined && visual) visual.setVisibility(state.visible)
        if (state.pickable !== undefined && visual) visual.setPickable(state.pickable)
        if (state.transform !== undefined && visual) visual.setTransform(state.transform)

        Representation.updateState(_state, state)
    }

    function setTheme(theme: Theme) {
        _theme = theme
    }

    function destroy() {
        if (visual) visual.destroy()
    }

    return {
        label,
        get groupCount() {
            return visual ? visual.groupCount : 0
        },
        get renderObjects() {
            return visual && visual.renderObject ? [ visual.renderObject ] : []
        },
        get props() { return _props },
        get params() { return _params },
        get state() { return _state },
        get theme() { return _theme },
        updated,
        createOrUpdate,
        setState,
        setTheme,
        getLoci,
        mark,
        destroy
    }
}