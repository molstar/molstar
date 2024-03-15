/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { Choice } from '../../../mol-util/param-choice';
import { Color } from '../../../mol-util/color';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

export const VolumeTypeChoice = new Choice({ 'isosurface': 'Isosurface', 'direct-volume': 'Direct volume', 'slice': 'Slice', 'off': 'Off' }, 'isosurface');
export type VolumeType = Choice.Values<typeof VolumeTypeChoice>


export const VolsegStateParams = {
    // segmentOpacity: PD.Numeric(1, { min: 0, max: 1, step: 0.05 }),
    // segmentKey: `${kind}:${segmentationId}:${segmentId}`
    selectedSegment: PD.Text(''),
    // visibleSegments: PD.ObjectList({
    //     segmentId: PD.Numeric(0),
    //     segmentationId: PD.Text(''),
    //     kind: PD.Select('lattice', [['lattice', 'lattice'], ['mesh', 'mesh'], ['primitive', 'primitive']])
    // }, k => `${k.segmentId}:${k.segmentationId}:${k.kind}`),
    visibleSegments: PD.ObjectList({
        segmentKey: PD.Text('') }, k => k.segmentKey
    ),
    visibleModels: PD.ObjectList({ pdbId: PD.Text('') }, s => s.pdbId.toString()),
    channelsData: PD.ObjectList({
        channelId: PD.Text('0'),
        volumeIsovalueKind: PD.Select('relative', [['relative', 'Relative'], ['absolute', 'Absolute']]),
        volumeIsovalueValue: PD.Numeric(1),
        volumeType: VolumeTypeChoice.PDSelect(),
        volumeOpacity: PD.Numeric(0.2, { min: 0, max: 1, step: 0.05 }),
        label: PD.Text(''),
        color: PD.Color(Color(0x121212))
    },
    i => i.channelId
    )
};
export type VolsegStateData = PD.Values<typeof VolsegStateParams>;


export class VolsegState extends PluginStateObject.Create<VolsegStateData>({ name: 'Vol & Seg Entry State', typeClass: 'Data' }) { }


export const VOLSEG_STATE_FROM_ENTRY_TRANSFORMER_NAME = 'volseg-state-from-entry'; // defined here to avoid cyclic dependency
