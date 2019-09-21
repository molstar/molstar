/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from '../../../../mol-model/structure';
import { VolumeServerInfo } from './model';
import { PluginContext } from '../../../../mol-plugin/context';
import { RuntimeContext } from '../../../../mol-task';
import { getXMLNodeByName, XMLDocument } from '../../../../mol-util/xml-parser';

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

export function getId(s?: Structure): string {
    if (!s) return ''

    const model = s.models[0]
    if (model.sourceData.kind !== 'mmCIF') return ''

    const d = model.sourceData.data
    for (let i = 0, il = d.pdbx_database_related._rowCount; i < il; ++i) {
        if (d.pdbx_database_related.db_name.value(i).toUpperCase() === 'EMDB') {
            return d.pdbx_database_related.db_id.value(i)
        }
    }

    return s.models.length > 0 ? s.models[0].entryId : ''
}

export async function getContourLevel(provider: 'wwpdb' | 'pdbe', plugin: PluginContext, taskCtx: RuntimeContext, emdbId: string) {
    switch (provider) {
        case 'wwpdb': return getContourLevelWwpdb(plugin, taskCtx, emdbId)
        case 'pdbe': return getContourLevelPdbe(plugin, taskCtx, emdbId)
    }
}

export async function getContourLevelWwpdb(plugin: PluginContext, taskCtx: RuntimeContext, emdbId: string) {
    // TODO: parametrize to a differnt URL? in plugin settings perhaps
    const header = await plugin.fetch<XMLDocument>({ url: `https://ftp.wwpdb.org/pub/emdb/structures/${emdbId.toUpperCase()}/header/${emdbId.toLowerCase()}.xml`, type: 'xml' }).runInContext(taskCtx);

    const map = getXMLNodeByName('map', header!.root!.children!)!
    const contourLevel = parseFloat(getXMLNodeByName('contourLevel', map.children!)!.content!)

    return contourLevel;
}

export async function getContourLevelPdbe(plugin: PluginContext, taskCtx: RuntimeContext, emdbId: string) {
    emdbId = emdbId.toUpperCase()
    // TODO: parametrize to a differnt URL? in plugin settings perhaps
    const header = await plugin.fetch({ url: `https://www.ebi.ac.uk/pdbe/api/emdb/entry/map/${emdbId}`, type: 'json' }).runInContext(taskCtx);
    const emdbEntry = header && header[emdbId];
    let contourLevel: number | undefined = void 0;
    if (emdbEntry && emdbEntry[0] && emdbEntry[0].map && emdbEntry[0].map.contour_level && emdbEntry[0].map.contour_level.value !== void 0) {
        contourLevel = +emdbEntry[0].map.contour_level.value;
    }

    return contourLevel;
}

export async function getEmdbId(plugin: PluginContext, taskCtx: RuntimeContext, pdbId: string) {
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

    return emdbId
}