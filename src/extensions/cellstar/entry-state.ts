/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { PluginStateObject } from '../../mol-plugin-state/objects';
import { ParamDefinition as PD } from '../../mol-util/param-definition';

import { Choice } from './helpers';


export const VolumeTypeChoice = new Choice({ 'isosurface': 'Isosurface', 'direct-volume': 'Direct volume', 'off': 'Off' }, 'isosurface');
export type VolumeType = Choice.Values<typeof VolumeTypeChoice>


export const CellstarStateParams = {
    volumeType: VolumeTypeChoice.PDSelect(),
    volumeIsovalueKind: PD.Select('relative', [['relative', 'Relative'], ['absolute', 'Absolute']]),
    volumeIsovalueValue: PD.Numeric(1),
    volumeOpacity: PD.Numeric(0.2, { min: 0, max: 1, step: 0.05 }),
    segmentOpacity: PD.Numeric(1, { min: 0, max: 1, step: 0.05 }),
    selectedSegment: PD.Numeric(-1, { step: 1 }),
    visibleSegments: PD.ObjectList({ segmentId: PD.Numeric(0) }, s => s.segmentId.toString()),
    visibleModels: PD.ObjectList({ pdbId: PD.Text('') }, s => s.pdbId.toString()),
};
export type CellstarStateData = PD.Values<typeof CellstarStateParams>;


export class CellstarState extends PluginStateObject.Create<CellstarStateData>({ name: 'Vol & Seg Entry State', typeClass: 'Data' }) { }


export const CELLSTAR_STATE_FROM_ENTRY_TRANSFORMER_NAME = 'cellstar-state-from-entry'; // defined here to avoid cyclic dependency
