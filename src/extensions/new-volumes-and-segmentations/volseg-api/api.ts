/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { DescriptionData, type Metadata } from './data';


export const DEFAULT_VOLSEG_SERVER = 'https://molstarvolseg.ncbr.muni.cz/v2';


export class VolumeApiV2 {
    public volumeServerUrl: string;

    public constructor(volumeServerUrl: string = DEFAULT_VOLSEG_SERVER) {
        this.volumeServerUrl = volumeServerUrl.replace(/\/$/, ''); // trim trailing slash
    }
    // @app.post("/v2/{source}/{id}/descriptions/edit")
    public async editDescriptionsUrl(source: string, entryId: string, descriptionData: DescriptionData[]) {
        const url = `${this.volumeServerUrl}/${source}/${entryId}/descriptions/edit`;
        const obj = JSON.stringify({ descriptions: descriptionData });
        console.log(obj);
        debugger;
        const response = await fetch(url, {
            method: 'POST',
            // body: JSON.stringify({notification: {title: message},to : '/topics/user_'+username}),
            body: obj,
            // headers: {'Content-Type': 'application/json', 'Authorization': 'key='+API_KEY}
            headers: { 'Content-Type': 'application/json' } 
        });
    }
    public entryListUrl(maxEntries: number, keyword?: string): string {
        return `${this.volumeServerUrl}/list_entries/${maxEntries}/${keyword ?? ''}`;
    }

    public metadataUrl(source: string, entryId: string): string {
        return `${this.volumeServerUrl}/${source}/${entryId}/metadata`;
    }
    public volumeUrl(source: string, entryId: string, timeframe: number, channelId: string, box: [[number, number, number], [number, number, number]] | null, maxPoints: number): string {
        if (box) {
            const [[a1, a2, a3], [b1, b2, b3]] = box;
            return `${this.volumeServerUrl}/${source}/${entryId}/volume/box/${timeframe}/${channelId}/${a1}/${a2}/${a3}/${b1}/${b2}/${b3}?max_points=${maxPoints}`;
        } else {
            return `${this.volumeServerUrl}/${source}/${entryId}/volume/cell/${timeframe}/${channelId}?max_points=${maxPoints}`;
        }
    }

    // @app.get("/v2/{source}/{id}/segmentation/box/{segmentation}/{time}/{a1}/{a2}/{a3}/{b1}/{b2}/{b3}")
    public latticeUrl(source: string, entryId: string, segmentation: string, timeframe: number, box: [[number, number, number], [number, number, number]] | null, maxPoints: number): string {
        if (box) {
            const [[a1, a2, a3], [b1, b2, b3]] = box;
            return `${this.volumeServerUrl}/${source}/${entryId}/segmentation/box/${segmentation}/${timeframe}/${a1}/${a2}/${a3}/${b1}/${b2}/${b3}?max_points=${maxPoints}`;
        } else {
            return `${this.volumeServerUrl}/${source}/${entryId}/segmentation/cell/${segmentation}/${timeframe}?max_points=${maxPoints}`;
        }
    }
    public geometricSegmentationUrl(source: string, entryId: string, segmentation_id: string, timeframe: number) {
        return `${this.volumeServerUrl}/${source}/${entryId}/geometric_segmentation/${segmentation_id}/${timeframe}`;
    }
    // @app.get("/v2/{source}/{id}/mesh/{segmentation_id}/{time}/{segment_id}/{detail_lvl}")
    public meshUrl_Json(source: string, entryId: string, segmentation_id: string, timeframe: number, segment: number, detailLevel: number): string {
        return `${this.volumeServerUrl}/${source}/${entryId}/mesh/${segmentation_id}/${timeframe}/${segment}/${detailLevel}`;
    }

    public meshUrl_Bcif(source: string, entryId: string, segmentation_id: string, timeframe: number, segment: number, detailLevel: number): string {
        return `${this.volumeServerUrl}/${source}/${entryId}/mesh_bcif/${segmentation_id}/${timeframe}/${segment}/${detailLevel}`;
    }

    public volumeInfoUrl(source: string, entryId: string): string {
        return `${this.volumeServerUrl}/${source}/${entryId}/volume_info`;
    }

    public async getEntryList(maxEntries: number, keyword?: string): Promise<{ [source: string]: string[] }> {
        const response = await fetch(this.entryListUrl(maxEntries, keyword));
        return await response.json();
    }

    public async getMetadata(source: string, entryId: string): Promise<Metadata> {
        const url = this.metadataUrl(source, entryId);
        const response = await fetch(url);
        if (!response.ok) throw new Error(`Failed to fetch metadata from ${url}`);
        return await response.json();
    }
}
