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
import { getQualityProps } from '../util';
import { StructureProps, DefaultStructureProps, StructureRepresentation } from '.';
import { ComplexVisual } from './complex-visual';

export function ComplexRepresentation<P extends StructureProps>(visualCtor: () => ComplexVisual<P>): StructureRepresentation<P> {
    let visual: ComplexVisual<P> | undefined

    let _props: P
    let _structure: Structure

    function create(structure: Structure, props: Partial<P> = {}) {
        _props = Object.assign({}, DefaultStructureProps, _props, props, getQualityProps(props, structure))
        _props.colorTheme.structure = structure

        return Task.create('Creating StructureRepresentation', async ctx => {
            if (!_structure) {
                visual = visualCtor()
                await visual.create(ctx, structure, _props)
            } else {
                if (_structure.hashCode === structure.hashCode) {
                    await update(_props)
                } else {
                    if (visual && !await visual.update(ctx, _props)) {
                        await visual.create(ctx, structure, _props)
                    }
                }
            }
            _structure = structure
        });
    }

    function update(props: Partial<P>) {
        return Task.create('Updating StructureRepresentation', async ctx => {
            _props = Object.assign({}, DefaultStructureProps, _props, props, getQualityProps(props, _structure))
            _props.colorTheme.structure = _structure

            if (visual && !await visual.update(ctx, _props)) {
                await visual.create(ctx, _structure, _props)
            }
        })
    }

    function getLoci(pickingId: PickingId) {
        return visual ? visual.getLoci(pickingId) : EmptyLoci
    }

    function mark(loci: Loci, action: MarkerAction) {
        if (visual) visual.mark(loci, action)
    }

    function destroy() {
        if (visual) visual.destroy()
    }

    return {
        get renderObjects() {
            return visual && visual.renderObject ? [ visual.renderObject ] : []
        },
        get props() { return _props },
        create,
        update,
        getLoci,
        mark,
        destroy
    }
}