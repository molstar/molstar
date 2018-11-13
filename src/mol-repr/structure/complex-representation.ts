/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from 'mol-model/structure';
import { Task } from 'mol-task'
import { Loci, EmptyLoci } from 'mol-model/loci';
import { StructureProps, StructureRepresentation } from './representation';
import { ComplexVisual } from './complex-visual';
import { PickingId } from 'mol-geo/geometry/picking';
import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { RepresentationContext } from 'mol-repr/representation';
import { Theme, ThemeProps, createTheme, ThemeRegistryContext } from 'mol-theme/theme';
import { ParamDefinition as PD } from 'mol-util/param-definition';

export function ComplexRepresentation<P extends StructureProps>(label: string, getParams: (ctx: ThemeRegistryContext, data: Structure) => PD.Params, visualCtor: () => ComplexVisual<P>): StructureRepresentation<P> {
    let visual: ComplexVisual<P> | undefined

    let _structure: Structure
    let _params: PD.Params
    let _props: P
    let _theme: Theme

    function createOrUpdate(ctx: RepresentationContext, props: Partial<P> = {}, themeProps: ThemeProps = {}, structure?: Structure) {
        if (structure && structure !== _structure) {
            _params = getParams(ctx, structure)
            _structure = structure
            if (!_props) _props = PD.getDefaultValues(_params) as P
        }
        _props = Object.assign({}, _props, props)
        _theme = createTheme(ctx, { structure: _structure }, props, themeProps, _theme)

        return Task.create('Creating or updating ComplexRepresentation', async runtime => {
            if (!visual) visual = visualCtor()
            await visual.createOrUpdate({ ...ctx, runtime }, _theme, _props, structure)
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
        createOrUpdate,
        getLoci,
        mark,
        destroy
    }
}