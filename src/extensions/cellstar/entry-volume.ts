import { Volume } from '../../mol-model/volume';
import { createVolumeRepresentationParams } from '../../mol-plugin-state/helpers/volume-representation-params';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { Download } from '../../mol-plugin-state/transforms/data';
import { CreateGroup } from '../../mol-plugin-state/transforms/misc';
import { StateObjectSelector } from '../../mol-state';
import { Color } from '../../mol-util/color';

import { CellStarEntryData, MAX_VOXELS } from './entry-root';
import * as ExternalAPIs from './external-api';


export class CellStarVolumeData {
    private entryData: CellStarEntryData;
    public volume?: Volume;

    constructor(rootData: CellStarEntryData) {
        this.entryData = rootData;
    }

    async showVolume() {
        const hasVolumes = this.entryData.metadata.grid.volumes.volume_downsamplings.length > 0;
        if (hasVolumes) {
            const group = await this.entryData.groupNodeMgr.showNode('Volume',
                async () => await this.entryData.newUpdate().apply(CreateGroup, { label: 'Volume' }, { state: { isCollapsed: true } }).commit(),
                false
            );
            const url = this.entryData.api.volumeUrl(this.entryData.source, this.entryData.entryId, null, MAX_VOXELS);
            const data = await this.entryData.newUpdate().to(group).apply(Download, { url, isBinary: true, label: `Volume Data: ${url}` }).commit();
            const parsed = await this.entryData.plugin.dataFormats.get('dscif')!.parse(this.entryData.plugin, data);
            const volume: StateObjectSelector<PluginStateObject.Volume.Data> = parsed.volumes?.[0] ?? parsed.volume;
            this.volume = volume.cell!.obj!.data;
            const isoLevel = await ExternalAPIs.getIsovalue(this.entryData.metadata.grid.general.source_db_id ?? this.entryData.entryId);
            await this.entryData.plugin.build()
                .to(volume)
                .apply(StateTransforms.Representation.VolumeRepresentation3D, createVolumeRepresentationParams(this.entryData.plugin, this.volume, {
                    type: 'isosurface',
                    typeParams: { alpha: 0.2, isoValue: Volume.adjustedIsoValue(this.volume, isoLevel.value, isoLevel.kind) },
                    color: 'uniform',
                    colorParams: { value: Color(0x121212) }
                }))
                .commit();
        }
    }
}