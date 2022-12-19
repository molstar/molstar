/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Vec2 } from '../../mol-math/linear-algebra';
import { Volume } from '../../mol-model/volume';
import { createVolumeRepresentationParams } from '../../mol-plugin-state/helpers/volume-representation-params';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { Download } from '../../mol-plugin-state/transforms/data';
import { CreateGroup } from '../../mol-plugin-state/transforms/misc';
import { setSubtreeVisibility } from '../../mol-plugin/behavior/static/state';
import { PluginCommands } from '../../mol-plugin/commands';
import { StateObjectSelector } from '../../mol-state';
import { Color } from '../../mol-util/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';

import { BOX, VolsegEntryData, MAX_VOXELS } from './entry-root';
import { VolsegStateParams, VolumeTypeChoice } from './entry-state';
import * as ExternalAPIs from './external-api';
import { VolsegGlobalStateData } from './global-state';


const GROUP_TAG = 'volume-group';
export const VOLUME_VISUAL_TAG = 'volume-visual';

const DIRECT_VOLUME_RELATIVE_PEAK_HALFWIDTH = 0.5;


export type VolumeVisualParams = ReturnType<typeof createVolumeRepresentationParams>;

interface VolumeStats { min: number, max: number, mean: number, sigma: number };


export const SimpleVolumeParams = {
    volumeType: VolumeTypeChoice.PDSelect(),
    opacity: PD.Numeric(0.2, { min: 0, max: 1, step: 0.05 }, { hideIf: p => p.volumeType === 'off' }),
};
export type SimpleVolumeParamValues = PD.Values<typeof SimpleVolumeParams>;


export class VolsegVolumeData {
    private entryData: VolsegEntryData;
    private visualTypeParamCache: { [type: string]: any } = {};

    constructor(rootData: VolsegEntryData) {
        this.entryData = rootData;
    }

    async loadVolume() {
        const hasVolumes = this.entryData.metadata.raw.grid.volumes.volume_downsamplings.length > 0;
        if (hasVolumes) {
            const isoLevelPromise = ExternalAPIs.tryGetIsovalue(this.entryData.metadata.raw.grid.general.source_db_id ?? this.entryData.entryId);
            let group = this.entryData.findNodesByTags(GROUP_TAG)[0]?.transform.ref;
            if (!group) {
                const newGroupNode = await this.entryData.newUpdate().apply(CreateGroup, { label: 'Volume' }, { tags: [GROUP_TAG], state: { isCollapsed: true } }).commit();
                group = newGroupNode.ref;
            }
            const url = this.entryData.api.volumeUrl(this.entryData.source, this.entryData.entryId, BOX, MAX_VOXELS);
            const data = await this.entryData.newUpdate().to(group).apply(Download, { url, isBinary: true, label: `Volume Data: ${url}` }).commit();
            const parsed = await this.entryData.plugin.dataFormats.get('dscif')!.parse(this.entryData.plugin, data);
            const volumeNode: StateObjectSelector<PluginStateObject.Volume.Data> = parsed.volumes?.[0] ?? parsed.volume;
            const volumeData = volumeNode.cell!.obj!.data;

            const volumeType = VolsegStateParams.volumeType.defaultValue;
            let isovalue = await isoLevelPromise;
            if (!isovalue) {
                const stats = volumeData.grid.stats;
                const maxRelative = (stats.max - stats.mean) / stats.sigma;
                if (maxRelative > 1) {
                    isovalue = { kind: 'relative', value: 1.0 };
                } else {
                    isovalue = { kind: 'relative', value: maxRelative * 0.5 };
                }
            }

            const adjustedIsovalue = Volume.adjustedIsoValue(volumeData, isovalue.value, isovalue.kind);
            const visualParams = this.createVolumeVisualParams(volumeData, volumeType);
            this.changeIsovalueInVolumeVisualParams(visualParams, adjustedIsovalue, volumeData.grid.stats);

            await this.entryData.newUpdate()
                .to(volumeNode)
                .apply(StateTransforms.Representation.VolumeRepresentation3D, visualParams, { tags: [VOLUME_VISUAL_TAG], state: { isHidden: volumeType === 'off' } })
                .commit();
            return { isovalue: adjustedIsovalue };
        }
    }

    async setVolumeVisual(type: 'isosurface' | 'direct-volume' | 'off') {
        const visual = this.entryData.findNodesByTags(VOLUME_VISUAL_TAG)[0];
        if (!visual) return;
        const oldParams: VolumeVisualParams = visual.transform.params;
        this.visualTypeParamCache[oldParams.type.name] = oldParams.type.params;
        if (type === 'off') {
            setSubtreeVisibility(this.entryData.plugin.state.data, visual.transform.ref, true); // true means hide, ¯\_(ツ)_/¯
        } else {
            setSubtreeVisibility(this.entryData.plugin.state.data, visual.transform.ref, false); // true means hide, ¯\_(ツ)_/¯
            if (oldParams.type.name === type) return;
            const newParams: VolumeVisualParams = {
                ...oldParams,
                type: {
                    name: type,
                    params: this.visualTypeParamCache[type] ?? oldParams.type.params,
                }
            };
            const volumeStats = visual.obj?.data.sourceData.grid.stats;
            if (!volumeStats) throw new Error(`Cannot get volume stats from volume visual ${visual.transform.ref}`);
            this.changeIsovalueInVolumeVisualParams(newParams, undefined, volumeStats);
            const update = this.entryData.newUpdate().to(visual.transform.ref).update(newParams);
            await PluginCommands.State.Update(this.entryData.plugin, { state: this.entryData.plugin.state.data, tree: update, options: { doNotUpdateCurrent: true } });
        }
    }

    async updateVolumeVisual(newParams: SimpleVolumeParamValues) {
        const { volumeType, opacity } = newParams;
        const visual = this.entryData.findNodesByTags(VOLUME_VISUAL_TAG)[0];
        if (!visual) return;
        const oldVisualParams: VolumeVisualParams = visual.transform.params;
        this.visualTypeParamCache[oldVisualParams.type.name] = oldVisualParams.type.params;

        if (volumeType === 'off') {
            setSubtreeVisibility(this.entryData.plugin.state.data, visual.transform.ref, true); // true means hide, ¯\_(ツ)_/¯
        } else {
            setSubtreeVisibility(this.entryData.plugin.state.data, visual.transform.ref, false); // true means hide, ¯\_(ツ)_/¯
            const newVisualParams: VolumeVisualParams = {
                ...oldVisualParams,
                type: {
                    name: volumeType,
                    params: this.visualTypeParamCache[volumeType] ?? oldVisualParams.type.params,
                }
            };
            newVisualParams.type.params.alpha = opacity;
            const volumeStats = visual.obj?.data.sourceData.grid.stats;
            if (!volumeStats) throw new Error(`Cannot get volume stats from volume visual ${visual.transform.ref}`);
            this.changeIsovalueInVolumeVisualParams(newVisualParams, undefined, volumeStats);
            const update = this.entryData.newUpdate().to(visual.transform.ref).update(newVisualParams);
            await PluginCommands.State.Update(this.entryData.plugin, { state: this.entryData.plugin.state.data, tree: update, options: { doNotUpdateCurrent: true } });
        }
    }

    async setTryUseGpu(tryUseGpu: boolean) {
        const visuals = this.entryData.findNodesByTags(VOLUME_VISUAL_TAG);
        for (const visual of visuals) {
            const oldParams: VolumeVisualParams = visual.transform.params;
            if (oldParams.type.params.tryUseGpu === !tryUseGpu) {
                const newParams = { ...oldParams, type: { ...oldParams.type, params: { ...oldParams.type.params, tryUseGpu: tryUseGpu } } };
                const update = this.entryData.newUpdate().to(visual.transform.ref).update(newParams);
                await PluginCommands.State.Update(this.entryData.plugin, { state: this.entryData.plugin.state.data, tree: update, options: { doNotUpdateCurrent: true } });
            }
        }
    }

    private getIsovalueFromState(): Volume.IsoValue {
        const { volumeIsovalueKind, volumeIsovalueValue } = this.entryData.currentState.value;
        return volumeIsovalueKind === 'relative'
            ? Volume.IsoValue.relative(volumeIsovalueValue)
            : Volume.IsoValue.absolute(volumeIsovalueValue);
    }

    private createVolumeVisualParams(volume: Volume, type: 'isosurface' | 'direct-volume' | 'off'): VolumeVisualParams {
        if (type === 'off') type = 'isosurface';
        return createVolumeRepresentationParams(this.entryData.plugin, volume, {
            type: type,
            typeParams: { alpha: 0.2, tryUseGpu: VolsegGlobalStateData.getGlobalState(this.entryData.plugin)?.tryUseGpu },
            color: 'uniform',
            colorParams: { value: Color(0x121212) },
        });
    }

    private changeIsovalueInVolumeVisualParams(params: VolumeVisualParams, isovalue: Volume.IsoValue | undefined, stats: VolumeStats) {
        isovalue ??= this.getIsovalueFromState();
        switch (params.type.name) {
            case 'isosurface':
                params.type.params.isoValue = isovalue;
                params.type.params.tryUseGpu = VolsegGlobalStateData.getGlobalState(this.entryData.plugin)?.tryUseGpu;
                break;
            case 'direct-volume':
                const absIso = Volume.IsoValue.toAbsolute(isovalue, stats).absoluteValue;
                const fractIso = (absIso - stats.min) / (stats.max - stats.min);
                const peakHalfwidth = DIRECT_VOLUME_RELATIVE_PEAK_HALFWIDTH * stats.sigma / (stats.max - stats.min);
                params.type.params.controlPoints = [
                    Vec2.create(Math.max(fractIso - peakHalfwidth, 0), 0),
                    Vec2.create(fractIso, 1),
                    Vec2.create(Math.min(fractIso + peakHalfwidth, 1), 0),
                ];
                break;
        }
    }
}