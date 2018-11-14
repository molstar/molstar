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
import { RepresentationContext, RepresentationParamsGetter } from 'mol-repr/representation';
import { Theme, ThemeProps, createTheme } from 'mol-theme/theme';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { BehaviorSubject } from 'rxjs';

export function ComplexRepresentation<P extends StructureParams>(label: string, getParams: RepresentationParamsGetter<Structure, P>, visualCtor: () => ComplexVisual<P>): StructureRepresentation<P> {
    const updated = new BehaviorSubject(0)
    let visual: ComplexVisual<P> | undefined

    let _structure: Structure
    let _params: P
    let _props: PD.DefaultValues<P>
    let _theme: Theme

    function createOrUpdate(ctx: RepresentationContext, props: Partial<PD.DefaultValues<P>> = {}, themeProps: ThemeProps = {}, structure?: Structure) {
        if (structure && structure !== _structure) {
            _params = getParams(ctx, structure)
            _structure = structure
            if (!_props) _props = PD.getDefaultValues(_params)
        }
        _props = Object.assign({}, _props, props)
        _theme = createTheme(ctx, { structure: _structure }, props, themeProps, _theme)

        return Task.create('Creating or updating ComplexRepresentation', async runtime => {
            if (!visual) visual = visualCtor()
            await visual.createOrUpdate({ ...ctx, runtime }, _theme, _props, structure)
            updated.next(updated.getValue() + 1)
        });
    }

    function getLoci(pickingId: PickingId) {
        return visual ? visual.getLoci(pickingId) : EmptyLoci
    }

    function mark(loci: Loci, action: MarkerAction) {
        return visual ? visual.mark(loci, action) : false
    }

    function destroy() {
        if (visual) visual.destroy()
    }

    return {
        label,
        get renderObjects() {
            return visual && visual.renderObject ? [ visual.renderObject ] : []
        },
        get props() { return _props },
        get params() { return _params },
        get updated() { return updated },
        createOrUpdate,
        getLoci,
        mark,
        destroy
    }
}