/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BuilderAction } from '../base';
import { StateTransformer } from '../../../mol-state';
import { Download as DownloadData, ReadFile } from '../../transforms/data';

export const Download = BuilderAction((builder, params: StateTransformer.Params<DownloadData>, { options }) => {
    return builder.apply(DownloadData, params, options);
});

export const OpenFile = BuilderAction((builder, params: StateTransformer.Params<ReadFile>, { options }) => {
    return builder.apply(ReadFile, params, options);
});