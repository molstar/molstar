/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from 'mol-model/structure';
import { Task } from 'mol-task'
import { PickingId } from '../../util/picking';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { MarkerAction } from '../../util/marker-data';
import { StructureProps, StructureRepresentation } from '.';
import { ComplexVisual } from './complex-visual';

export function ComplexRepresentation<P extends StructureProps>(label: string, visualCtor: () => ComplexVisual<P>): StructureRepresentation<P> {
    let visual: ComplexVisual<P> | undefined
    let _props: P

    function createOrUpdate(props: Partial<P> = {}, structure?: Structure) {
        _props = Object.assign({}, _props, props)

        return Task.create('Creating or updating ComplexRepresentation', async ctx => {
            if (!visual) visual = visualCtor()
            await visual.createOrUpdate(ctx, _props, structure)
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
        createOrUpdate,
        getLoci,
        mark,
        destroy
    }
}