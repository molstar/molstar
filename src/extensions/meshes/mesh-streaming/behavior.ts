/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { distinctUntilChanged, map } from 'rxjs';

import { CIF } from '../../../mol-io/reader/cif';
import { Box3D } from '../../../mol-math/geometry';
import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { PluginBehavior } from '../../../mol-plugin/behavior';
import { PluginCommand } from '../../../mol-plugin/command';
import { PluginCommands } from '../../../mol-plugin/commands';
import { PluginContext } from '../../../mol-plugin/context';
import { UUID } from '../../../mol-util';
import { Asset } from '../../../mol-util/assets';
import { Color } from '../../../mol-util/color';
import { ColorNames } from '../../../mol-util/color/names';
import { Choice } from '../../../mol-util/param-choice';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';

import { MetadataWrapper } from '../../volumes-and-segmentations/volseg-api/utils';

import { MeshlistData } from '../mesh-extension';
import { MeshServerInfo } from './server-info';


const DEFAULT_SEGMENT_NAME = 'Untitled segment';
const DEFAULT_SEGMENT_COLOR = ColorNames.lightgray;
export const NO_SEGMENT = -1;
/** Maximum (worst) detail level available in GUI (TODO set actual maximum possible value) */
const MAX_DETAIL = 10;
const DEFAULT_DETAIL = 7; // TODO decide a reasonable default
/** Segments whose bounding box volume is above this value (relative to the overall bounding box) are considered as background segments */
export const BACKGROUND_SEGMENT_VOLUME_THRESHOLD = 0.5;


export class MeshStreaming extends PluginStateObject.CreateBehavior<MeshStreaming.Behavior>({ name: 'Mesh Streaming' }) { }

export namespace MeshStreaming {

    export namespace Params {
        export const ViewTypeChoice = new Choice({ off: 'Off', select: 'Select', all: 'All' }, 'select'); // TODO add camera target?
        export type ViewType = Choice.Values<typeof ViewTypeChoice>;

        export function create(options: MeshServerInfo.Data) {
            return {
                view: PD.MappedStatic('select', {
                    'off': PD.Group({}),
                    'select': PD.Group({
                        baseDetail: PD.Numeric(DEFAULT_DETAIL, { min: 1, max: MAX_DETAIL, step: 1 }, { description: 'Detail level for the non-selected segments (lower number = better)' }),
                        focusDetail: PD.Numeric(1, { min: 1, max: MAX_DETAIL, step: 1 }, { description: 'Detail level for the selected segment (lower number = better)' }),
                        selectedSegment: PD.Numeric(NO_SEGMENT, {}, { isHidden: true }),
                    }, { isFlat: true }),
                    'all': PD.Group({
                        detail: PD.Numeric(DEFAULT_DETAIL, { min: 1, max: MAX_DETAIL, step: 1 }, { description: 'Detail level for all segments (lower number = better)' })
                    }, { isFlat: true }),
                }, { description: '"Off" hides all segments. \n"Select" shows all segments in lower detail, clicked segment in better detail. "All" shows all segment in the same level.' }),
            };
        }

        export type Definition = ReturnType<typeof create>
        export type Values = PD.Values<Definition>

        export function copyValues(params: Values): Values {
            return {
                view: {
                    name: params.view.name,
                    params: { ...params.view.params } as any,
                }
            };
        }
        export function valuesEqual(p: Values, q: Values): boolean {
            if (p.view.name !== q.view.name) return false;
            for (const key in p.view.params) {
                if ((p.view.params as any)[key] !== (q.view.params as any)[key]) return false;
            }
            return true;
        }
        export function detailsEqual(p: Values, q: Values): boolean {
            switch (p.view.name) {
                case 'off':
                    return q.view.name === 'off';
                case 'select':
                    return q.view.name === 'select' && p.view.params.baseDetail === q.view.params.baseDetail && p.view.params.focusDetail === q.view.params.focusDetail;
                case 'all':
                    return q.view.name === 'all' && p.view.params.detail === q.view.params.detail;
                default:
                    throw new Error('Not implemented');
            }
        }
    }

    export interface VisualInfo {
        tag: string, // e.g. high-2, low-1 // ? remove if can be omitted
        segmentId: number, // ? remove if unused
        segmentName: string, // ? remove if unused
        detailType: VisualInfo.DetailType, // ? remove if unused
        detail: number, // ? remove if unused
        color: Color, // move to MeshlistData?
        visible: boolean,
        data?: MeshlistData,
    }
    export namespace VisualInfo {
        export type DetailType = 'low' | 'high';
        export const DetailTypes: DetailType[] = ['low', 'high'];
        export function tagFor(segmentId: number, detail: DetailType) {
            return `${detail}-${segmentId}`;
        }
    }


    export class Behavior extends PluginBehavior.WithSubscribers<Params.Values> {
        private id: string;
        private ref: string = '';
        public parentData: MeshServerInfo.Data;
        private metadata?: MetadataWrapper;
        public visuals?: { [tag: string]: VisualInfo };
        public backgroundSegments: { [segmentId: number]: boolean } = {};
        private focusObservable = this.plugin.behaviors.interaction.click.pipe( // QUESTION is this OK way to get focused segment?
            map(evt => evt.current.loci),
            map(loci => (loci.kind === 'group-loci') ? loci.shape.sourceData as MeshlistData : null),
            map(data => (data?.ownerId === this.id) ? data : null), // do not process shapes created by others
            distinctUntilChanged((old, current) => old?.segmentId === current?.segmentId),
        );
        private focusSubscription?: PluginCommand.Subscription = undefined;
        private backgroundSegmentsInitialized = false;

        constructor(plugin: PluginContext, data: MeshServerInfo.Data, params: Params.Values) {
            super(plugin, params);
            this.id = UUID.create22();
            this.parentData = data;
        }

        register(ref: string): void {
            this.ref = ref;
        }

        unregister(): void {
            if (this.focusSubscription) {
                this.focusSubscription.unsubscribe();
                this.focusSubscription = undefined;
            }
            // TODO empty cache here (if used)
        }

        selectSegment(segmentId: number) {
            if (this.params.view.name === 'select') {
                if (this.params.view.params.selectedSegment === segmentId) return;
                const newParams = Params.copyValues(this.params);
                if (newParams.view.name === 'select') {
                    newParams.view.params.selectedSegment = segmentId;
                }
                const state = this.plugin.state.data;
                const update = state.build().to(this.ref).update(newParams);
                PluginCommands.State.Update(this.plugin, { state, tree: update, options: { doNotUpdateCurrent: true } });
            }
        }

        async update(params: Params.Values) {
            const oldParams = this.params;
            this.params = params;

            if (!this.metadata) {
                const response = await fetch(this.getMetadataUrl());
                const rawMetadata = await response.json();
                this.metadata = new MetadataWrapper(rawMetadata);
            }

            if (!this.visuals) {
                this.initVisualInfos();
            } else if (!Params.detailsEqual(this.params, oldParams)) {
                this.updateVisualInfoDetails();
            }

            switch (params.view.name) {
                case 'off':
                    await this.disableVisuals();
                    break;
                case 'select':
                    await this.enableVisuals(params.view.params.selectedSegment);
                    break;
                case 'all':
                    await this.enableVisuals();
                    break;
                default:
                    throw new Error('Not implemented');
            }
            if (params.view.name !== 'off' && !this.backgroundSegmentsInitialized) {
                this.guessBackgroundSegments();
                this.backgroundSegmentsInitialized = true;
            }
            if (params.view.name === 'select' && !this.focusSubscription) {
                this.focusSubscription = this.subscribeObservable(this.focusObservable, data => { this.selectSegment(data?.segmentId ?? NO_SEGMENT); });
            } else if (params.view.name !== 'select' && this.focusSubscription) {
                this.focusSubscription.unsubscribe();
                this.focusSubscription = undefined;
            }
            return true;
        }

        private getMetadataUrl() {
            return `${this.parentData.serverUrl}/${this.parentData.source}/${this.parentData.entryId}/metadata`;
        }

        private getMeshUrl(segment: number, detail: number) {
            return `${this.parentData.serverUrl}/${this.parentData.source}/${this.parentData.entryId}/mesh_bcif/${segment}/${detail}`;
        }

        private initVisualInfos() {
            const visuals: { [tag: string]: VisualInfo } = {};
            for (const segid of this.metadata!.meshSegmentIds) {
                const name = this.metadata?.getSegment(segid)?.biological_annotation.name ?? DEFAULT_SEGMENT_NAME;
                const color = this.metadata?.getSegmentColor(segid) ?? DEFAULT_SEGMENT_COLOR;
                for (const detailType of VisualInfo.DetailTypes) {
                    const visual: VisualInfo = {
                        tag: VisualInfo.tagFor(segid, detailType),
                        segmentId: segid,
                        segmentName: name,
                        detailType: detailType,
                        detail: -1, // to be set at the end
                        color: color,
                        visible: false,
                        data: undefined,
                    };
                    visuals[visual.tag] = visual;
                }
            }
            this.visuals = visuals;
            this.updateVisualInfoDetails();
        }
        private updateVisualInfoDetails() {
            let highDetail: number | undefined;
            let lowDetail: number | undefined;
            switch (this.params.view.name) {
                case 'off':
                    lowDetail = undefined;
                    highDetail = undefined;
                    break;
                case 'select':
                    lowDetail = this.params.view.params.baseDetail;
                    highDetail = this.params.view.params.focusDetail;
                    break;
                case 'all':
                    lowDetail = this.params.view.params.detail;
                    highDetail = undefined;
                    break;
            }
            for (const tag in this.visuals) {
                const visual = this.visuals[tag];
                const preferredDetail = (visual.detailType === 'high') ? highDetail : lowDetail;
                if (preferredDetail !== undefined) {
                    visual.detail = this.metadata!.getSufficientMeshDetail(visual.segmentId, preferredDetail);
                }
            }
        }

        private async enableVisuals(highDetailSegment?: number) {
            for (const tag in this.visuals) {
                const visual = this.visuals[tag];
                const requiredDetailType = visual.segmentId === highDetailSegment ? 'high' : 'low';
                if (visual.detailType === requiredDetailType) {
                    visual.data = await this.getMeshData(visual);
                    visual.visible = true;
                } else {
                    visual.visible = false;
                }
            }
        }

        private async disableVisuals() {
            for (const tag in this.visuals) {
                const visual = this.visuals[tag];
                visual.visible = false;
            }
        }

        /** Fetch data in current `visual.detail`, or return already fetched data (if available in the correct detail). */
        private async getMeshData(visual: VisualInfo): Promise<MeshlistData> {
            if (visual.data && visual.data.detail === visual.detail) {
                // Do not recreate
                return visual.data;
            }
            // TODO cache
            const url = this.getMeshUrl(visual.segmentId, visual.detail);
            const urlAsset = Asset.getUrlAsset(this.plugin.managers.asset, url);
            const asset = await this.plugin.runTask(this.plugin.managers.asset.resolve(urlAsset, 'binary'));
            const parsed = await this.plugin.runTask(CIF.parseBinary(asset.data));
            if (parsed.isError) {
                throw new Error(`Failed parsing CIF file from ${url}`);
            }
            const meshlistData = await MeshlistData.fromCIF(parsed.result, visual.segmentId, visual.segmentName, visual.detail);
            meshlistData.ownerId = this.id;
            // const bbox = MeshlistData.bbox(meshlistData);
            // const bboxVolume = bbox ? MS.Box3D.volume(bbox) : 0.0;
            // console.log(`BBox ${visual.segmentId}: ${Math.round(bboxVolume! / 1e6)} M`, bbox); // DEBUG
            return meshlistData;
        }

        private async guessBackgroundSegments() {
            const bboxes: { [segid: number]: Box3D } = {};
            for (const tag in this.visuals) {
                const visual = this.visuals[tag];
                if (visual.detailType === 'low' && visual.data) {
                    const bbox = MeshlistData.bbox(visual.data);
                    if (bbox) {
                        bboxes[visual.segmentId] = bbox;
                    }
                }
            }
            const totalBbox = MeshlistData.combineBBoxes(Object.values(bboxes));
            const totalVolume = totalBbox ? Box3D.volume(totalBbox) : 0.0;
            // console.log(`BBox total: ${Math.round(totalVolume! / 1e6)} M`, totalBbox); // DEBUG

            const isBgSegment: { [segid: number]: boolean } = {};
            for (const segid in bboxes) {
                const bbox = bboxes[segid];
                const bboxVolume = Box3D.volume(bbox);
                isBgSegment[segid] = (bboxVolume > totalVolume * BACKGROUND_SEGMENT_VOLUME_THRESHOLD);
                // console.log(`BBox ${segid}: ${Math.round(bboxVolume! / 1e6)} M, ${bboxVolume / totalVolume}`, bbox); // DEBUG
            }
            this.backgroundSegments = isBgSegment;
        }

        getDescription() {
            return Params.ViewTypeChoice.prettyName(this.params.view.name);
        }

    }
}

