/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { View } from '../view';
import { TransformListController } from '../../controller/transform/list';
import { UrlEntity } from 'mol-view/state/entity';
import { ModelToStructure, StructureToBallAndStick, StructureToSpacefill, StructureToDistanceRestraint, StructureToBackbone, MmcifUrlToModel, StructureToCartoon } from 'mol-view/state/transform';
import { StateContext } from 'mol-view/state/context';
import { SpacefillProps } from 'mol-geo/representation/structure/representation/spacefill';
import { BallAndStickProps } from 'mol-geo/representation/structure/representation/ball-and-stick';
import { DistanceRestraintProps } from 'mol-geo/representation/structure/representation/distance-restraint';
import { BackboneProps } from 'mol-geo/representation/structure/representation/backbone';
import { CartoonProps } from 'mol-geo/representation/structure/representation/cartoon';

const spacefillProps: SpacefillProps = {
    doubleSided: true,
    colorTheme: { name: 'chain-id' },
    quality: 'auto',
    useFog: false
}

const ballAndStickProps: BallAndStickProps = {
    doubleSided: true,
    colorTheme: { name: 'chain-id' },
    sizeTheme: { name: 'uniform', value: 0.05 },
    linkRadius: 0.05,
    quality: 'auto',
    useFog: false
}

const distanceRestraintProps: DistanceRestraintProps = {
    doubleSided: true,
    colorTheme: { name: 'chain-id' },
    linkRadius: 0.5,
    quality: 'auto',
    useFog: false
}

const backboneProps: BackboneProps = {
    doubleSided: true,
    colorTheme: { name: 'chain-id' },
    quality: 'auto',
    useFog: false
}

const cartoonProps: CartoonProps = {
    doubleSided: true,
    colorTheme: { name: 'chain-id' },
    quality: 'auto',
    useFog: false
}

function getPdbdevUrl(pdbdevId: string) {
    return `https://pdb-dev.wwpdb.org/static/cif/${pdbdevId}.cif`
}

const exampleUrls = {
    PDBDEV_00000001: getPdbdevUrl('PDBDEV_00000001'), // ok
    PDBDEV_00000002: getPdbdevUrl('PDBDEV_00000002'), // ok
    PDBDEV_00000003: getPdbdevUrl('PDBDEV_00000003'), // ok
    PDBDEV_00000004: getPdbdevUrl('PDBDEV_00000004'), // TODO issue with cross-link extraction
    PDBDEV_00000005: getPdbdevUrl('PDBDEV_00000005'), // ok
    PDBDEV_00000006: getPdbdevUrl('PDBDEV_00000006'), // TODO only three spacefill atoms rendered
    PDBDEV_00000007: getPdbdevUrl('PDBDEV_00000007'), // TODO only three spacefill atoms rendered
    PDBDEV_00000008: getPdbdevUrl('PDBDEV_00000008'), // ok
    PDBDEV_00000010: getPdbdevUrl('PDBDEV_00000010'), // ok
    PDBDEV_00000011: getPdbdevUrl('PDBDEV_00000011'), // ok
    PDBDEV_00000012: getPdbdevUrl('PDBDEV_00000012'), // ok
    PDBDEV_00000014: getPdbdevUrl('PDBDEV_00000014'), // ok
    PDBDEV_00000016: getPdbdevUrl('PDBDEV_00000016'),
}

export class UrlLoader extends View<TransformListController, { }, { ctx: StateContext }> {
    async load(name: keyof typeof exampleUrls) {
        console.log(exampleUrls[name])
        const ctx = this.props.ctx
        const urlEntity = UrlEntity.ofUrl(ctx, exampleUrls[name])
        console.log(await urlEntity.value.getData())
        const modelEntity = await MmcifUrlToModel.apply(ctx, urlEntity)
        const structureEntity = await ModelToStructure.apply(ctx, modelEntity)

        StructureToBallAndStick.apply(ctx, structureEntity, { ...ballAndStickProps, visible: true })
        StructureToSpacefill.apply(ctx, structureEntity, { ...spacefillProps, visible: false })
        StructureToDistanceRestraint.apply(ctx, structureEntity, { ...distanceRestraintProps, visible: false })
        StructureToBackbone.apply(ctx, structureEntity, { ...backboneProps, visible: true })
        StructureToCartoon.apply(ctx, structureEntity, { ...cartoonProps, visible: false })
    }

    render() {
        const exampleOptions = Object.keys(exampleUrls).map(name => {
            return <option key={name} value={name}>{name}</option>
        })

        return <div className='molstar-transformer-wrapper'>
            <div className='molstar-panel molstar-control molstar-transformer molstar-panel-expanded'>
                <div className='molstar-panel-body'>
                    <div>
                        <div className='molstar-control-row molstar-options-group'>
                            <span>Examples</span>
                            <div>
                                <select
                                    className='molstar-form-control'
                                    value=''
                                    onChange={(e) => {
                                        if (e.target.value) {
                                            this.load(e.target.value as keyof typeof exampleUrls)}
                                        }
                                    }
                                >
                                    <option key='' value=''></option>
                                    {exampleOptions}
                                </select>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>;
    }
}