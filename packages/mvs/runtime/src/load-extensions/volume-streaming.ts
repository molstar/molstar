/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { VolumeStreaming } from '@molstar/plugin/behavior/dynamic/volume-streaming/behavior';
import { CreateVolumeStreamingBehavior, CreateVolumeStreamingInfo, VolumeStreamingVisual } from '@molstar/plugin/behavior/dynamic/volume-streaming/transformers';
import { mapObjectMap } from '@molstar/core/util/object';
import { decodeColor } from '@molstar/mvs/helpers/utils';
import type { MolstarLoadingExtension } from '@molstar/mvs/load';
import { UpdateTarget } from '@molstar/mvs/load-generic';
import { ColorT } from '@molstar/mvs-builder/tree/mvs/param-types';


/** Type of `molstar_volume_streaming` custom property, used by `VolumeStreamingExtension` MVS loading extension. */
export type MolstarVolumeStreamingCustomProp = {
    /** URL of the volume streaming server, e.g. 'https://www.ebi.ac.uk/pdbe/densities'. */
    server_url?: string,
    /** Volume streaming view type ('off' | 'box' | 'selection-box' | 'camera-target' | 'cell' | 'auto'). Default value depends on structure type (X-ray/EM). */
    view?: VolumeStreaming.ViewTypes,
    /** Customization of channel parameters. */
    channel_params?: { [name in VolumeStreaming.ChannelType]?: Partial<ChannelParams_> },
    /** List of volume streaming entries (if not specified, will be retrieved automatically based on PDB ID) */
    entries?: ReturnType<typeof CreateVolumeStreamingInfo['createDefaultParams']>['entries'],
} | boolean | undefined;


/** This MVS loading extension allows turning on volume streaming for a structure by providing custom property `molstar_volume_streaming`.
 *
 *  Examples:
 *
 *  ```
 *  builder
 *      .download({ url: 'https://www.ebi.ac.uk/pdbe/entry-files/download/1cbs_updated.cif' })
 *      .parse({ format: 'mmcif' })
 *      .modelStructure({
 *          custom: {
 *              molstar_volume_streaming: true,
 *          },
 *      })
 *      .component()
 *      .representation();
 *
 *  builder
 *      .download({ url: 'https://www.ebi.ac.uk/pdbe/entry-files/download/1tqn_updated.cif' })
 *      .parse({ format: 'mmcif' })
 *      .modelStructure({
 *          custom: {
 *              molstar_volume_streaming: {
 *                  channel_params: {
 *                      '2fo-fc': { color: 'skyblue', opacity: 0.3 },
 *                      'fo-fc(+ve)': { color: 'greenyellow', wireframe: true, isoValue: { kind: 'relative', relativeValue: +2.5 } },
 *                      'fo-fc(-ve)': { color: 'orange', wireframe: true, isoValue: { kind: 'relative', relativeValue: -2.5 } },
 *                  },
 *              } satisfies MolstarVolumeStreamingCustomProp,
 *          },
 *      })
 *      .component()
 *      .representation();
 *
 *  builder
 *      .download({ url: 'https://www.ebi.ac.uk/pdbe/entry-files/download/8hra_updated.cif' })
 *      .parse({ format: 'mmcif' })
 *      .modelStructure({
 *          custom: {
 *              molstar_volume_streaming: {
 *                  server_url: 'https://www.ebi.ac.uk/pdbe/densities', // = default
 *                  entries: [{ dataId: 'EMD-34965', source: { name: 'em', params: { isoValue: { kind: 'absolute', absoluteValue: 0.015 } } } }],
 *                  view: 'auto', // default is 'auto' for EM, 'selection-box' for X-ray structures
 *                  channel_params: {
 *                      em: { color: '#ff0000', opacity: 0.4, isoValue: { kind: 'absolute', absoluteValue: 0.025 } },
 *                  },
 *              } satisfies MolstarVolumeStreamingCustomProp,
 *          },
 *      })
 *      .component()
 *      .representation();
 *  ```
 */
export const VolumeStreamingExtension: MolstarLoadingExtension<{}> = {
    id: 'wwpdb/volume-streaming',
    description: 'Allow turning on volume streaming for a structure',
    createExtensionContext: () => ({}),
    action: (updateTarget, node, context, extContext) => {
        if (node.kind !== 'structure') return;
        let params: MolstarVolumeStreamingCustomProp = node.custom?.molstar_volume_streaming;
        if (!params) return;
        if (params === true) params = {};

        const streamingInfo = UpdateTarget.apply(updateTarget, CreateVolumeStreamingInfo, {
            serverUrl: params.server_url,
            autoEntries: !params.entries,
            entries: params.entries,
            defaultView: params.view,
            defaultChannelParams: params.channel_params && mapObjectMap(params.channel_params, normalizeChannelParams),
        }, { state: { isCollapsed: true } });

        const streamingBehavior = UpdateTarget.apply(streamingInfo, CreateVolumeStreamingBehavior);

        UpdateTarget.apply(streamingBehavior, VolumeStreamingVisual, { channel: '2fo-fc' }, { state: { isGhost: true }, tags: '2fo-fc' });
        UpdateTarget.apply(streamingBehavior, VolumeStreamingVisual, { channel: 'fo-fc(+ve)' }, { state: { isGhost: true }, tags: 'fo-fc(+ve)' });
        UpdateTarget.apply(streamingBehavior, VolumeStreamingVisual, { channel: 'fo-fc(-ve)' }, { state: { isGhost: true }, tags: 'fo-fc(-ve)' });
        UpdateTarget.apply(streamingBehavior, VolumeStreamingVisual, { channel: 'em' }, { state: { isGhost: true }, tags: 'em' });
    },
};


interface ChannelParams_ extends Omit<VolumeStreaming.ChannelParams, 'color'> {
    color: ColorT | number,
}

function normalizeChannelParams(p: Partial<ChannelParams_> | undefined): Partial<VolumeStreaming.ChannelParams> | undefined {
    if (!p) return undefined;
    return {
        ...p,
        color: decodeColor(p.color),
    };
}
