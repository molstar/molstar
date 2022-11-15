import { type Metadata } from './data';

function getProcess() {
    try {
        return process; // `process` is only available in Node
    } catch {
        return undefined;
    }
}
const PROCESS = getProcess();

function createApiPrefix() {
    const hostname = PROCESS?.env.REACT_APP_API_HOSTNAME
        ? PROCESS?.env.REACT_APP_API_HOSTNAME : `${window.location.protocol}//${window.location.hostname}`;
    const port = PROCESS?.env.REACT_APP_API_PORT
        ? PROCESS?.env.REACT_APP_API_PORT : '9000';
    const prefix = PROCESS?.env.REACT_APP_API_PREFIX
        ? `/${PROCESS?.env.REACT_APP_API_PREFIX}` : ``;

    return `${hostname}:${port}${prefix}`;
}

function getGitTag() {
    return `${PROCESS?.env.REACT_APP_GIT_TAG ?? ''}`;
}

function getGitSha() {
    return `${PROCESS?.env.REACT_APP_GIT_SHA ?? ''}`;
}

const DEFAULT_API_PREFIX = createApiPrefix();
const GIT_TAG = getGitTag();
const GIT_SHA = getGitSha();
export const DEFAULT_VOLUME_SERVER_V2 = `${DEFAULT_API_PREFIX}/v2`;


export class VolumeApiV2 {
    public volumeServerUrl: string;
    public volumeServerGitTag: string;
    public volumeServerGitSha: string;
    
    public constructor(
        volumeServerUrl: string = DEFAULT_VOLUME_SERVER_V2,
        volumeServerGitTag: string = GIT_TAG,
        volumeServerGitSha: string = GIT_SHA
        ) {
        this.volumeServerUrl = volumeServerUrl.replace(/\/$/, '');  // trim trailing slash
        this.volumeServerGitTag = volumeServerGitTag;
        this.volumeServerGitSha = volumeServerGitSha;

        console.log('API V2', this.volumeServerUrl)
        console.log(`SHA: ${this.volumeServerGitSha}`)
        console.log(`GIT TAG: ${this.volumeServerGitTag}`)
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
        return await response.json();
    }

    public async getMetadata(source: string, entryId: string): Promise<Metadata> {
        const url = this.metadataUrl(source, entryId);
        const response = await fetch(url);
        if (!response.ok) throw new Error(`Failed to fetch metadata from ${url}`);
        return await response.json();
    }
}
