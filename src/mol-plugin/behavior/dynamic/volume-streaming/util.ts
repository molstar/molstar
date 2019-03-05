/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure } from 'mol-model/structure';
import { VolumeServerInfo } from './model';
import { PluginContext } from 'mol-plugin/context';
import { RuntimeContext } from 'mol-task';

export function getStreamingMethod(s?: Structure, defaultKind: VolumeServerInfo.Kind = 'x-ray'): VolumeServerInfo.Kind {
    if (!s) return defaultKind;

    const model = s.models[0];
    if (model.sourceData.kind !== 'mmCIF') return defaultKind;

    const data = model.sourceData.data.exptl.method;
    for (let i = 0; i < data.rowCount; i++) {
        const v = data.value(i).toUpperCase();
        if (v.indexOf('MICROSCOPY') >= 0) return 'em';
    }
    return 'x-ray';
}

export async function getEmdbIdAndContourLevel(plugin: PluginContext, taskCtx: RuntimeContext, pdbId: string) {
    // TODO: parametrize to a differnt URL? in plugin settings perhaps
    const summary = await plugin.fetch({ url: `https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/${pdbId}`, type: 'json' }).runInContext(taskCtx);

    const summaryEntry = summary && summary[pdbId];
    let emdbId: string;
    if (summaryEntry && summaryEntry[0] && summaryEntry[0].related_structures) {
        const emdb = summaryEntry[0].related_structures.filter((s: any) => s.resource === 'EMDB');
        if (!emdb.length) {
            throw new Error(`No related EMDB entry found for '${pdbId}'.`);
        }
        emdbId = emdb[0].accession;
    } else {
        throw new Error(`No related EMDB entry found for '${pdbId}'.`);
    }

    // TODO: parametrize to a differnt URL? in plugin settings perhaps
    const emdb = await plugin.fetch({ url: `https://www.ebi.ac.uk/pdbe/api/emdb/entry/map/${emdbId}`, type: 'json' }).runInContext(taskCtx);
    const emdbEntry = emdb && emdb[emdbId];
    let contour: number | undefined = void 0;
    if (emdbEntry && emdbEntry[0] && emdbEntry[0].map && emdbEntry[0].map.contour_level && emdbEntry[0].map.contour_level.value !== void 0) {
        contour = +emdbEntry[0].map.contour_level.value;
    }

    return { emdbId, contour };
}