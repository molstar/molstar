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


export class CellstarVolumeData {
    private entryData: CellstarEntryData;
    public volume?: Volume;

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
            await this.entryData.plugin.build()
                .to(volumeNode)
                .apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(this.entryData.plugin, volumeData, {
                    type: 'isosurface',
                    typeParams: { alpha: 0.2, isoValue: adjustedIsovalue },
                    color: 'uniform',
                    colorParams: { value: Color(0x121212) }
                }), { tags: ['volume-visual'] })
                .commit();
            return { isovalue: adjustedIsovalue };
        }
    }

    async setVolumeVisual(type: 'isosurface' | 'direct-volume' | 'off') {
        const visual = this.entryData.findNodesByTags('volume-visual')[0]?.transform;
        if (!visual) return;
        if (type === 'off') {
            setSubtreeVisibility(this.entryData.plugin.state.data, visual.ref, true); // true means hide, ¯\_(ツ)_/¯
        } else {
            setSubtreeVisibility(this.entryData.plugin.state.data, visual.ref, false); // true means hide, ¯\_(ツ)_/¯
            const oldParams: ReturnType<typeof createVolumeRepresentationParams> = visual.params;
            if (oldParams.type.name === type) return;
            console.log('volume visual params', visual.params);
            // const update = this.entryData.newUpdate().to(visual.ref).update(visual.params);
            // createVolumeRepresentationParams(this.en)
            const newParams = { ...oldParams, type: { ...oldParams.type, name: type } };
            const update = this.entryData.newUpdate().to(visual.ref).update(newParams);
            // const update = this.entryData.newUpdate().to(visual.ref).update({ type: { name: type, params: {} } });
            await PluginCommands.State.Update(this.entryData.plugin, { state: this.entryData.plugin.state.data, tree: update, options: { doNotUpdateCurrent: true } });
        }
    }

    // private createVolumeIsosurfaceParams() {
    //     const { volumeIsovalueKind, volumeIsovalueValue } = this.entryData.currentState.value;
    //     const iso: Volume.IsoValue = volumeIsovalueKind === 'relative'
    //         ? Volume.IsoValue.relative(volumeIsovalueValue)
    //         : Volume.IsoValue.absolute(volumeIsovalueValue);
    //     return createVolumeRepresentationParams(this.entryData.plugin, this.volume, {
    //         type: 'isosurface',
    //         typeParams: { alpha: 0.2, isoValue: Volume.adjustedIsoValue(this.volume, isoLevel.value, isoLevel.kind) },
    //         color: 'uniform',
    //         colorParams: { value: Color(0x121212) }
    //     }),

    // }
}