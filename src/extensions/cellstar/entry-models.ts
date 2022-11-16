import { BehaviorSubject } from 'rxjs';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { Download, RawData } from '../../mol-plugin-state/transforms/data';
import { CreateGroup } from '../../mol-plugin-state/transforms/misc';
import { StateObjectSelector } from '../../mol-state';
import { Asset } from '../../mol-util/assets';

import { CellStarEntryData } from './entry-root';
import { NodeManager } from './helpers';


export class CellStarModelData {
    private entryData: CellStarEntryData;
    /*private*/ pdbModelNodeMgr = new NodeManager();
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
        await this.entryData.plugin.builders.structure.hierarchy.applyPreset(trajectoryNode, 'default');
        return dataNode;
    }

    async showPdb(pdbId: string | undefined) {
        console.log(pdbId, 'Nodes:', this.pdbModelNodeMgr.getNodes());
        this.pdbModelNodeMgr.hideAllNodes();
        if (pdbId) {
            const update = this.entryData.newUpdate();
            const group = await this.entryData.groupNodeMgr.showNode('Fitted Models', () => update.apply(CreateGroup, {label: 'Fitted Models'}).selector, false);
            await update.commit();
            await this.pdbModelNodeMgr.showNode(pdbId, async () => await this.loadPdb(pdbId, group));
        }
        this.currentPdb.next(pdbId);
    }


}