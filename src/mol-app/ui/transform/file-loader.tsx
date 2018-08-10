/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react'
import { View } from '../view';
import { FileInput } from '../controls/common';
import { TransformListController } from '../../controller/transform/list';
import { FileEntity } from 'mol-view/state/entity';
import { MmcifFileToModel, ModelToStructure, StructureToBallAndStick, StructureToSpacefill, StructureToDistanceRestraint, StructureToBackbone } from 'mol-view/state/transform';
import { StateContext } from 'mol-view/state/context';
import { SpacefillProps } from 'mol-geo/representation/structure/representation/spacefill';
import { BallAndStickProps } from 'mol-geo/representation/structure/representation/ball-and-stick';
import { DistanceRestraintProps } from 'mol-geo/representation/structure/representation/distance-restraint';
import { BackboneProps } from 'mol-geo/representation/structure/representation/backbone';

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

export class FileLoader extends View<TransformListController, {}, { ctx: StateContext }> {
    render() {
        return <div className='molstar-file-loader'>
            <FileInput
                accept='*.cif'
                onChange={async files => {
                    if (files) {
                        const ctx = this.props.ctx
                        const fileEntity = FileEntity.ofFile(ctx, files[0])
                        const modelEntity = await MmcifFileToModel.apply(ctx, fileEntity)
                        const structureEntity = await ModelToStructure.apply(ctx, modelEntity)

                        StructureToBallAndStick.apply(ctx, structureEntity, { ...ballAndStickProps, visible: true })
                        StructureToSpacefill.apply(ctx, structureEntity, { ...spacefillProps, visible: false })
                        StructureToDistanceRestraint.apply(ctx, structureEntity, { ...distanceRestraintProps, visible: false })
                        StructureToBackbone.apply(ctx, structureEntity, { ...backboneProps, visible: true })
                    }
                }}
            />
        </div>;
    }
}