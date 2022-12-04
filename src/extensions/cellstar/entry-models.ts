import { BehaviorSubject } from 'rxjs';
import { RawData } from '../../mol-plugin-state/transforms/data';
import { CreateGroup } from '../../mol-plugin-state/transforms/misc';
import { StateObjectSelector } from '../../mol-state';
import { Asset } from '../../mol-util/assets';

import { CellStarEntryData } from './entry-root';
import { NodeManager } from './helpers';


export class CellStarModelData {
    private entryData: CellStarEntryData;
    private pdbModelNodeMgr = new NodeManager();
    currentPdb = new BehaviorSubject<string | undefined>(undefined);

    constructor(rootData: CellStarEntryData) {
        this.entryData = rootData;
    }


    private async loadPdb(pdbId: string, parent: StateObjectSelector) {
        const url = `https://www.ebi.ac.uk/pdbe/entry-files/download/${pdbId}.bcif`;
        const urlAsset = Asset.getUrlAsset(this.entryData.plugin.managers.asset, url);
        const asset = await this.entryData.plugin.runTask(this.entryData.plugin.managers.asset.resolve(urlAsset, 'binary'));
        const data = asset.data;

        const dataNode = await this.entryData.plugin.build().to(parent).apply(RawData, { data, label: `PDB Data: ${url}` }).commit();
        const trajectoryNode = await this.entryData.plugin.builders.structure.parseTrajectory(dataNode, 'mmcif');
        await this.entryData.plugin.builders.structure.hierarchy.applyPreset(trajectoryNode, 'default', { representationPreset: 'polymer-cartoon' });
        return dataNode;
    }

    async showPdb(pdbId: string | undefined) {
        console.log(pdbId, 'Nodes:', this.pdbModelNodeMgr.getNodes());
        this.pdbModelNodeMgr.hideAllNodes();
        if (pdbId) {
            // await update.commit();
            const group = await this.entryData.groupNodeMgr.showNode('FittedModels', async () => await this.entryData.newUpdate().apply(CreateGroup, { label: 'Fitted Models' }, { state: { isCollapsed: true } }).commit(), false);
            await this.pdbModelNodeMgr.showNode(pdbId, async () => await this.loadPdb(pdbId, group));
        }
        this.currentPdb.next(pdbId);
    }


}