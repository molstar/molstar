/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, Model } from '../../../../mol-model/structure';
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

/** Returns EMD ID when available, otherwise falls back to PDB ID */
export function getEmIds(model: Model): string[] {
    const ids: string[] = []
    if (model.sourceData.kind !== 'mmCIF') return [ model.entryId ]

    const { db_id, db_name, content_type } = model.sourceData.data.pdbx_database_related
    if (!db_name.isDefined) return [ model.entryId ]

    for (let i = 0, il = db_name.rowCount; i < il; ++i) {
        if (db_name.value(i).toUpperCase() === 'EMDB' && content_type.value(i) === 'associated EM volume') {
            ids.push(db_id.value(i))
        }
    }

    return ids
}

export function getXrayIds(model: Model): string[] {
    return [ model.entryId ]
}

export function getIds(method: VolumeServerInfo.Kind, s?: Structure): string[] {
    if (!s || !s.models.length) return []
    const model = s.models[0]
    switch (method) {
        case 'em': return getEmIds(model)
        case 'x-ray': return getXrayIds(model)
    }
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

export async function getEmdbIds(plugin: PluginContext, taskCtx: RuntimeContext, pdbId: string) {
    // TODO: parametrize to a differnt URL? in plugin settings perhaps
    const summary = await plugin.fetch({ url: `https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/${pdbId}`, type: 'json' }).runInContext(taskCtx);

    const summaryEntry = summary && summary[pdbId];
    let emdbIds: string[] = [];
    if (summaryEntry && summaryEntry[0] && summaryEntry[0].related_structures) {
        const emdb = summaryEntry[0].related_structures.filter((s: any) => s.resource === 'EMDB' && s.relationship === 'associated EM volume');
        if (!emdb.length) {
            throw new Error(`No related EMDB entry found for '${pdbId}'.`);
        }
        emdbIds.push(...emdb.map((e: { accession: string }) => e.accession));
    } else {
        throw new Error(`No related EMDB entry found for '${pdbId}'.`);
    }

    return emdbIds
}