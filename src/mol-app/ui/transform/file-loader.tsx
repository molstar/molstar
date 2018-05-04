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
import { MmcifFileToSpacefill } from 'mol-view/state/transform';
import { StateContext } from 'mol-view/state/context';

export class FileLoader extends View<TransformListController, {}, { ctx: StateContext }> {
    render() {
        return <div className='molstar-file-loader'>
            <FileInput
                accept='*.cif'
                onChange={files => {
                    if (files) {
                        const fileEntity = FileEntity.ofFile(this.props.ctx, files[0])
                        MmcifFileToSpacefill.apply(this.props.ctx, fileEntity)
                    }
                }}
            />
        </div>;
    }
}