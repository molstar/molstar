/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { type Metadata } from './data';


export const DEFAULT_VOLSEG_SERVER = 'https://molstarvolseg.ncbr.muni.cz/v2';


export class VolumeApiV2 {
    public volumeServerUrl: string;

    public constructor(volumeServerUrl: string = DEFAULT_VOLSEG_SERVER) {
        this.volumeServerUrl = volumeServerUrl.replace(/\/$/, ''); // trim trailing slash
    }

    public entryListUrl(maxEntries: number, keyword?: string): string {
        return `${this.volumeServerUrl}/list_entries/${maxEntries}/${keyword ?? ''}`;
    }

    public metadataUrl(source: string, entryId: string): string {
        return `${this.volumeServerUrl}/${source}/${entryId}/metadata`;
    }
    public volumeUrl(source: string, entryId: string, box: [[number, number, number], [number, number, number]] | null, maxPoints: number): string {
        if (box) {
            const [[a1, a2, a3], [b1, b2, b3]] = box;
            return `${this.volumeServerUrl}/${source}/${entryId}/volume/box/${a1}/${a2}/${a3}/${b1}/${b2}/${b3}?max_points=${maxPoints}`;
        } else {
            return `${this.volumeServerUrl}/${source}/${entryId}/volume/cell?max_points=${maxPoints}`;
        }
    }
    public latticeUrl(source: string, entryId: string, segmentation: number, box: [[number, number, number], [number, number, number]] | null, maxPoints: number): string {
        if (box) {
            const [[a1, a2, a3], [b1, b2, b3]] = box;
            return `${this.volumeServerUrl}/${source}/${entryId}/segmentation/box/${segmentation}/${a1}/${a2}/${a3}/${b1}/${b2}/${b3}?max_points=${maxPoints}`;
        } else {
            return `${this.volumeServerUrl}/${source}/${entryId}/segmentation/cell/${segmentation}?max_points=${maxPoints}`;
        }
    }
    public meshUrl_Json(source: string, entryId: string, segment: number, detailLevel: number): string {
        return `${this.volumeServerUrl}/${source}/${entryId}/mesh/${segment}/${detailLevel}`;
    }

    public meshUrl_Bcif(source: string, entryId: string, segment: number, detailLevel: number): string {
        return `${this.volumeServerUrl}/${source}/${entryId}/mesh_bcif/${segment}/${detailLevel}`;
    }
    public volumeInfoUrl(source: string, entryId: string): string {
        return `${this.volumeServerUrl}/${source}/${entryId}/volume_info`;
    }

    public async getEntryList(maxEntries: number, keyword?: string): Promise<{ [source: string]: string[] }> {
        const response = await fetch(this.entryListUrl(maxEntries, keyword));
        if (response.ok) {
            return await response.json();
        } else {
            console.error('Failed to fetch "Volume & Segmentation" entry list');
            return {};
        }
    }

    public async getMetadata(source: string, entryId: string): Promise<Metadata> {
        const url = this.metadataUrl(source, entryId);
        const response = await fetch(url);
        if (!response.ok) throw new Error(`Failed to fetch metadata from ${url}`);
        return await response.json();
    }
}
