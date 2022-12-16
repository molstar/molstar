import { Vec2 } from '../../mol-math/linear-algebra';
import { Volume } from '../../mol-model/volume';
import { createVolumeRepresentationParams } from '../../mol-plugin-state/helpers/volume-representation-params';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { Download } from '../../mol-plugin-state/transforms/data';
import { CreateGroup } from '../../mol-plugin-state/transforms/misc';
import { StateObjectSelector } from '../../mol-state';
import { Color } from '../../mol-util/color';
import { PluginCommands, setSubtreeVisibility } from '../meshes/molstar-lib-imports';

import { BOX, CellstarEntryData, MAX_VOXELS } from './entry-root';
import * as ExternalAPIs from './external-api';


type VolumeVisualParams = ReturnType<typeof createVolumeRepresentationParams>

export class CellstarVolumeData {
    private entryData: CellstarEntryData;
    private visualTypeParamCache: { [type: string]: any } = {};
    // public volume?: Volume;

    constructor(rootData: CellstarEntryData) {
        this.entryData = rootData;
    }

    async loadVolume() {
        const hasVolumes = this.entryData.metadata.raw.grid.volumes.volume_downsamplings.length > 0;
        if (hasVolumes) {
            const isoLevelPromise = ExternalAPIs.getIsovalue(this.entryData.metadata.raw.grid.general.source_db_id ?? this.entryData.entryId);
            let group = this.entryData.findNodesByTags('volume-group')[0]?.transform.ref;
            if (!group) {
                const newGroupNode = await this.entryData.newUpdate().apply(CreateGroup, { label: 'Volume' }, { tags: ['volume-group'], state: { isCollapsed: true } }).commit();
                group = newGroupNode.ref;
            }
            const url = this.entryData.api.volumeUrl(this.entryData.source, this.entryData.entryId, BOX, MAX_VOXELS);
            const data = await this.entryData.newUpdate().to(group).apply(Download, { url, isBinary: true, label: `Volume Data: ${url}` }).commit();
            const parsed = await this.entryData.plugin.dataFormats.get('dscif')!.parse(this.entryData.plugin, data);
            const volumeNode: StateObjectSelector<PluginStateObject.Volume.Data> = parsed.volumes?.[0] ?? parsed.volume;
            const volumeData = volumeNode.cell!.obj!.data;
            const isovalue = await isoLevelPromise;
            const adjustedIsovalue = Volume.adjustedIsoValue(volumeData, isovalue.value, isovalue.kind);
            // visualParams = this.initVolumeVisualTypeParamCache(volumeData, adjustedIsovalue);
            let visualParams = this.createVolumeVisualParams(volumeData, 'isosurface', adjustedIsovalue);
            visualParams = this.changeIsovalueInVolumeVisualParams(visualParams, adjustedIsovalue);
            console.log('visualParams', visualParams);
            await this.entryData.newUpdate()
                .to(volumeNode)
                .apply(StateTransforms.Representation.VolumeRepresentation3D, visualParams, { tags: ['volume-visual'] })
                .commit();
            return { isovalue: adjustedIsovalue };
        }
    }

    async setVolumeVisual(type: 'isosurface' | 'direct-volume' | 'off') {
        console.log('setVolumeVisual', type);
        const visual = this.entryData.findNodesByTags('volume-visual')[0]?.transform;
        if (!visual) return;
        const oldParams: VolumeVisualParams = visual.params;
        console.log('oldParams:', oldParams);
        this.visualTypeParamCache[oldParams.type.name] = oldParams.type.params;
        // console.log('param cache', this.visualTypeParamCache);
        if (type === 'off') {
            setSubtreeVisibility(this.entryData.plugin.state.data, visual.ref, true); // true means hide, ¯\_(ツ)_/¯
        } else {
            setSubtreeVisibility(this.entryData.plugin.state.data, visual.ref, false); // true means hide, ¯\_(ツ)_/¯
            if (oldParams.type.name === type) return;

            const newParams = {
                ...oldParams,
                type: {
                    name: type,
                    params: this.visualTypeParamCache[type] ?? oldParams.type.params,
                }
            }; // TODO or create params accordingly
            // newParams = this.changeIsovalueInVolumeVisualParams(newParams, undefined); // TODO uncomment
            // const newParams = { ...oldParams, type: { name: type, params: oldParams.type.params } };
            // console.log('volume visual params', visual.params);
            console.log('newParams:', newParams);
            let update = this.entryData.newUpdate().to(visual.ref).update(newParams);
            await PluginCommands.State.Update(this.entryData.plugin, { state: this.entryData.plugin.state.data, tree: update, options: { doNotUpdateCurrent: true } });
            let midParams = this.entryData.findNodesByTags('volume-visual')[0]?.transform.params;
            console.log('really set params:', midParams);
            midParams = this.changeIsovalueInVolumeVisualParams(midParams, undefined); // TODO uncomment

            update = this.entryData.newUpdate().to(visual.ref).update(midParams);
            await PluginCommands.State.Update(this.entryData.plugin, { state: this.entryData.plugin.state.data, tree: update, options: { doNotUpdateCurrent: true } });
            const finParams = this.entryData.findNodesByTags('volume-visual')[0]?.transform.params
            console.log('really set params:', finParams);

        }
    }

    private getIsovalueFromState(): Volume.IsoValue {
        const { volumeIsovalueKind, volumeIsovalueValue } = this.entryData.currentState.value;
        return volumeIsovalueKind === 'relative'
            ? Volume.IsoValue.relative(volumeIsovalueValue)
            : Volume.IsoValue.absolute(volumeIsovalueValue);
    }

    // private initVolumeVisualTypeParamCache(volume: Volume, isovalue: Volume.IsoValue) {
    //     this.visualTypeParamCache = {
    //         'isosurface':
    //             this.createVolumeVisualParams(volume, 'isosurface', isovalue).type.params,
    //         'direct-volume':
    //             this.createVolumeVisualParams(volume, 'direct-volume', isovalue).type.params,
    //     };
    // }

    private createVolumeVisualParams(volume: Volume, type: 'isosurface' | 'direct-volume', isovalue: Volume.IsoValue): VolumeVisualParams {
        switch (type) {
            case 'isosurface':
                return createVolumeRepresentationParams(this.entryData.plugin, volume, {
                    type: 'isosurface',
                    typeParams: { alpha: 0.2, isoValue: isovalue },
                    color: 'uniform',
                    colorParams: { value: Color(0x121212) },
                });
            case 'direct-volume':
                return createVolumeRepresentationParams(this.entryData.plugin, volume, {
                    type: 'direct-volume',
                    typeParams: { alpha: 0.2, controlPoints: [Vec2.create(0.4, 0.0), Vec2.create(0.5, 0.1), Vec2.create(0.6, 0.0)] }, // TODO smart
                    color: 'uniform',
                    colorParams: { value: Color(0x121212) },
                });
            default:
                throw new Error(`Invalid value for \`type\`: '${type}'`);
        }
    }
    private changeIsovalueInVolumeVisualParams(params: VolumeVisualParams, isovalue?: Volume.IsoValue): VolumeVisualParams {
        isovalue ??= this.getIsovalueFromState();
        switch (params.type.name) {
            case 'isosurface':
                return {
                    ...params,
                    type: {
                        name: params.type.name,
                        params: {
                            ...params.type.params,
                            isoValue: isovalue,
                        }
                    }
                }
            // params.type.params.isoValue = isovalue;
            // break;
            case 'direct-volume':
                return {
                    ...params,
                    type: {
                        name: params.type.name,
                        params: {
                            ...params.type.params,
                            controlPoints: [Vec2.create(0.4, 0.0), Vec2.create(0.5, 0.1), Vec2.create(0.6, 0.0)],
                        }
                    }
                }
            // params.type.params.controlPoints.length = 0;
            // params.type.params.controlPoints.push(Vec2.create(0.4, 0.0), Vec2.create(0.5, 0.1), Vec2.create(0.6, 0.0)); // TODO smart
            // // params.type.params.controlPoints = [Vec2.create(0.4, 0.0), Vec2.create(0.5, 0.1), Vec2.create(0.6, 0.0)]; // TODO smart
            // params.type.params.alpha = 0.5; //debug
            // break;
            default:
                return params;
        }
    }
}