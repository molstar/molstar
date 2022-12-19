/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Download, ParseCif } from '../../mol-plugin-state/transforms/data';
import { CreateGroup } from '../../mol-plugin-state/transforms/misc';
import { TrajectoryFromMmCif } from '../../mol-plugin-state/transforms/model';
import { setSubtreeVisibility } from '../../mol-plugin/behavior/static/state';
import { StateObjectRef, StateObjectSelector } from '../../mol-state';

import { VolsegEntryData } from './entry-root';


export class VolsegModelData {
    private entryData: VolsegEntryData;

    constructor(rootData: VolsegEntryData) {
        this.entryData = rootData;
    }

    private async loadPdb(pdbId: string, parent: StateObjectSelector | StateObjectRef) {
        const url = `https://www.ebi.ac.uk/pdbe/entry-files/download/${pdbId}.bcif`;
        const dataNode = await this.entryData.plugin.build().to(parent).apply(Download, { url: url, isBinary: true }, { tags: ['fitted-model-data', `pdbid-${pdbId}`] }).commit();
        const cifNode = await this.entryData.plugin.build().to(dataNode).apply(ParseCif).commit();
        const trajectoryNode = await this.entryData.plugin.build().to(cifNode).apply(TrajectoryFromMmCif).commit();
        await this.entryData.plugin.builders.structure.hierarchy.applyPreset(trajectoryNode, 'default', { representationPreset: 'polymer-cartoon' });
        return dataNode;
    }

    async showPdbs(pdbIds: string[]) {
        const segmentsToShow = new Set(pdbIds);

        const visuals = this.entryData.findNodesByTags('fitted-model-data');
        for (const visual of visuals) {
            const theTag = visual.obj?.tags?.find(tag => tag.startsWith('pdbid-'));
            if (!theTag) continue;
            const id = theTag.split('-')[1];
            const visibility = segmentsToShow.has(id);
            setSubtreeVisibility(this.entryData.plugin.state.data, visual.transform.ref, !visibility); // true means hide, ¯\_(ツ)_/¯
            segmentsToShow.delete(id);
        }

        const segmentsToCreate = Array.from(segmentsToShow);
        if (segmentsToCreate.length === 0) return;

        let group = this.entryData.findNodesByTags('fitted-models-group')[0]?.transform.ref;
        if (!group) {
            const newGroupNode = await this.entryData.newUpdate().apply(CreateGroup, { label: 'Fitted Models' }, { tags: ['fitted-models-group'], state: { isCollapsed: true } }).commit();
            group = newGroupNode.ref;
        }

        const awaiting = [];
        for (const pdbId of segmentsToCreate) {
            awaiting.push(this.loadPdb(pdbId, group));
        }
        for (const promise of awaiting) await promise;
    }
}
