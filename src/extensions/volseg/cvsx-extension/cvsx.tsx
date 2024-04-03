/**
 * Copyright (c) 2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Aliaksei Chareshneu <chareshneu@mail.muni.cz>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { useEffect, useRef } from 'react';
import { BehaviorSubject, combineLatest } from 'rxjs';
import { PluginComponent } from '../../../mol-plugin-state/component';
import { CollapsableControls, CollapsableState } from '../../../mol-plugin-ui/base';
import { Button, ControlRow, ExpandGroup } from '../../../mol-plugin-ui/controls/common';
import { GetAppSvg } from '../../../mol-plugin-ui/controls/icons';
import { useBehavior } from '../../../mol-plugin-ui/hooks/use-behavior';
import { PluginContext } from '../../../mol-plugin/context';
import { SimpleVolumeParamValues, SimpleVolumeParams, VolumeVisualParams } from '../new-volumes-and-segmentations/entry-volume';
import { UpdateTransformControl } from '../../../mol-plugin-ui/state/update-transform';
import { WaitingParameterControls, WaitingSlider } from '../new-volumes-and-segmentations/ui';
import { sleep } from '../../../mol-util/sleep';
import { StateTransform } from '../../../mol-state/transform';
import { setSubtreeVisibility } from '../../../mol-plugin/behavior/static/state';
import { PluginCommands } from '../../../mol-plugin/commands';
import { StateTransforms } from '../../../mol-plugin-state/transforms';
import { Metadata } from '../new-volumes-and-segmentations/volseg-api/data';
import { MetadataWrapper } from '../new-volumes-and-segmentations/volseg-api/utils';
import { findNodesByRef, findNodesByTags } from '../common';
import { StateHierarchyMirror } from '../new-volumes-and-segmentations/entry-root';
import { DescriptionsList } from '../common-ui';

export const CVSX_VOLUME_VISUAL_TAG = 'CVSX-volume-visual';
export const CVSX_LATTICE_SEGMENTATION_VISUAL_TAG = 'CVSX-lattice-segmentation-visual';
export const CVSX_ANNOTATIONS_FILE_TAG = 'CVSX-annotations-file';
export const CVSX_METADATA_FILE_TAG = 'CVSX-metadata-file';
export const CVSX_GEOMETRIC_SEGMENTATION_FILE = 'CVSX-geometric-segmentation-file';

export class CSVXUI extends CollapsableControls<{}, {}> {
    protected defaultState(): CollapsableState {
        return {
            header: 'CSVX',
            isCollapsed: true,
            brand: { accent: 'cyan', svg: GetAppSvg }
        };
    }
    protected renderControls(): JSX.Element | null {
        return <CVSXFileControls plugin={this.plugin} />;
    }
}

// export interface CVSXProps {
//     volumes: any | undefined,
//     segmentations: any | undefined,
//     annotations: AnnotationMetadata | undefined,
//     geometricSegmentation: ShapePrimitiveData | undefined
// }

export const CVSXState = {
    // segmentOpacity: PD.Numeric(1, { min: 0, max: 1, step: 0.05 }),
    // segmentKey: `${kind}:${segmentationId}:${segmentId}`
    selectedSegment: PD.Text(''),
    // visibleSegments: PD.ObjectList({
    //     segmentId: PD.Numeric(0),
    //     segmentationId: PD.Text(''),
    //     kind: PD.Select('lattice', [['lattice', 'lattice'], ['mesh', 'mesh'], ['primitive', 'primitive']])
    // }, k => `${k.segmentId}:${k.segmentationId}:${k.kind}`),
    visibleSegments: PD.ObjectList({
        segmentKey: PD.Text('')
    }, k => k.segmentKey
    ),
};

export type CVSXStateData = PD.Values<typeof CVSXState>;

// TODO: props could be volumes, segmentations, annotations
export class CVSXStateModel extends PluginComponent {
    // state = new BehaviorSubject<{ props: CVSXProps }>({ props: { volumes: undefined, segmentations: undefined, annotations: undefined, geometricSegmentation: undefined } });
    state = {
        hierarchy: new BehaviorSubject<StateHierarchyMirror | undefined>(undefined)
    };

    currentState = new BehaviorSubject(PD.getDefaultValues(CVSXState));
    metadata = new BehaviorSubject<MetadataWrapper | undefined>(undefined);
    public currentTimeframe = new BehaviorSubject(0);
    private visualTypeParamCache: { [type: string]: any } = {};

    _updateMetadata() {
        const annotationNodes = findNodesByTags(this.plugin, CVSX_ANNOTATIONS_FILE_TAG);
        const metadataNodes = findNodesByTags(this.plugin, CVSX_METADATA_FILE_TAG);
        const meta: Metadata = {
            grid: JSON.parse(metadataNodes[0]!.obj!.data),
            annotation: JSON.parse(annotationNodes[0]!.obj!.data)
        };
        this.metadata.next(new MetadataWrapper(meta));
        console.log('Metadata updated');
        console.log(this.metadata.value);
        // this.metadata = new MetadataWrapper(meta);
    }

    async updateStateNode(params: Partial<CVSXStateData>) {
        const oldParams = this.currentState.value;
        const newParams = { ...oldParams, ...params };
        const state = this.plugin.state.data;
        console.log('State was updated');
    }

    _updateHierarchy() {
        const volumes = findNodesByTags(this.plugin, CVSX_VOLUME_VISUAL_TAG);
        const segmentations = findNodesByTags(this.plugin, CVSX_LATTICE_SEGMENTATION_VISUAL_TAG);
        // const annotationNodes = findNodesByTags(this.plugin, CVSX_ANNOTATIONS_FILE_TAG);
        // let annotations = undefined;
        // let geometricSegmentation = undefined;
        // if (annotationNodes.length > 0) annotations = JSON.parse(annotationNodes[0]!.obj!.data);
        const geometricSegmentations = findNodesByTags(this.plugin, CVSX_GEOMETRIC_SEGMENTATION_FILE);
        // TODO: mesh segmentations
        const meshSegmentations = findNodesByTags(this.plugin, '');

        // if (geometricSegmentationNodes.length > 0) geometricSegmentation = JSON.parse(geometricSegmentationNodes[0]!.obj!.data);
        this.state.hierarchy.next({
            volumes, segmentations, geometricSegmentations, meshSegmentations
        });
    }

    mount() {
        this._updateHierarchy();
        this._updateMetadata();
        const obs = combineLatest([
            this.plugin.behaviors.state.isBusy,
            this.plugin.state.data.events.cell.stateUpdated
        ]);
        this.subscribe(obs, ([busy, cell]) => {
            if (busy) return;
            this._updateHierarchy();
            this._updateMetadata();
        });
    }

    changeColor(volumeSelector: any, newColor: any) {
        // trigger state upadte
    }

    doSomething = () => {

    };

    // NOTE: currently works for all segmentations at once
    updateSegmentationOpacity = async (opacity: number) => {
        const reprs = this.state.hierarchy?.value?.segmentations;
        if (!reprs) return;
        const update = this.plugin.build().toRoot();
        console.log(reprs);
        for (const s of reprs) {
            update.to(s).update(StateTransforms.Representation.VolumeRepresentation3D, p => { p.type.params.alpha = opacity; });
        }
        return await update.commit();
    };

    updateVolumeVisual = async (newParams: SimpleVolumeParamValues, transform: StateTransform) => {
        const { volumeType, opacity } = newParams;
        const visual = findNodesByRef(this.plugin, transform.ref);
        if (!visual) return;
        const oldVisualParams: VolumeVisualParams = visual.transform.params;
        this.visualTypeParamCache[oldVisualParams.type.name] = oldVisualParams.type.params;

        if (volumeType === 'off') {
            setSubtreeVisibility(this.plugin.state.data, visual.transform.ref, true); // true means hide, ¯\_(ツ)_/¯
        } else {
            setSubtreeVisibility(this.plugin.state.data, visual.transform.ref, false); // true means hide, ¯\_(ツ)_/¯
            const newVisualParams: VolumeVisualParams = {
                ...oldVisualParams,
                type: {
                    name: volumeType,
                    params: this.visualTypeParamCache[volumeType] ?? oldVisualParams.type.params,
                }
            };
            newVisualParams.type.params.alpha = opacity;
            const volumeStats = visual.obj?.data.sourceData.grid.stats;
            if (!volumeStats) throw new Error(`Cannot get volume stats from volume visual ${visual.transform.ref}`);
            // this.changeIsovalueInVolumeVisualParams(newVisualParams, undefined, volumeStats);
            const update = this.plugin.build().toRoot().to(visual.transform.ref).update(newVisualParams);
            await PluginCommands.State.Update(this.plugin, { state: this.plugin.state.data, tree: update, options: { doNotUpdateCurrent: true } });
        }
    };

    constructor(public plugin: PluginContext) {
        super();
    }
}

function CVSXFileControls({ plugin }: { plugin: PluginContext }) {
    const _model = useRef<CVSXStateModel>();
    if (!_model.current) {
        _model.current = new CVSXStateModel(plugin);
    }
    const model = _model.current;
    useEffect(() => {
        model.mount();
        return () => model.dispose();
    }, [model]);

    const isBusy = useBehavior(plugin.behaviors.state.isBusy);

    const metadata = useBehavior(model.metadata);
    // TODO: instead get descriptions for selected segmentation id and kind later
    // possibly need to store in hierarchy top level node
    // or include in bcif info about segmentation id and kind
    // it is alredy included in the form of file name (kind) but not segmentation id
    // include new category?
    // does not work
    // get from metadata?
    // how to associate specific segmentation with metadata
    const allDescriptions = metadata?.allDescriptions;
    // const allDescriptions = model.metadata.value.allDescriptions;

    const h = useBehavior(model.state.hierarchy);
    if (!h) return null;

    return <>
        disabled: {isBusy}
        {/* TODO: create props first */}
        {/* {state.props} */}
        {/* TODO: here render UI based on props */}
        {/* check how volume controls are rendered in volseg ui */}
        <>
            {h.volumes && <ExpandGroup header='Volume data' initiallyExpanded>
                {h.volumes.map(v => {
                    console.log('v', v);
                    const transform = v.transform;
                    if (!transform) return null;
                    const volumeValues: SimpleVolumeParamValues = {
                        volumeType: transform.state.isHidden ? 'off' : transform.params?.type.name as any,
                        opacity: transform.params?.type.params.alpha,
                    };
                    return <div key={v.transform.ref}>
                        <WaitingParameterControls params={SimpleVolumeParams} values={volumeValues} onChangeValues={async next => { await sleep(20); await model.updateVolumeVisual(next, transform) }} />
                        <UpdateTransformControl state={plugin.state.data} transform={transform} customHeader='none' />
                    </div>;

                })}
            </ExpandGroup>}
            {/* <SegmentationControls model={model}></SegmentationControls> */}
            {h.segmentations && <ExpandGroup header='Segmentation data' initiallyExpanded>
                {h.segmentations.map(s => {
                    // Opacity of segmentation can get from its visual
                    return <>
                        <ControlRow key={s.transform.ref} label='Opacity' control={
                            <WaitingSlider min={0} max={1} value={s.transform.params.type.params.alpha} step={0.05} onChange={async v => await model.updateSegmentationOpacity(v)} />
                        } />
                        <DescriptionsList model={model} targetSegmentationId={ } targetKind={'lattice'}></DescriptionsList>
                    </>
                    ;
                })}
                {/* {allDescriptions && <DescriptionsList
                    model={model} targetSegmentationId={}
                ></DescriptionsList>} */}
            </ExpandGroup>}

        </>


        <Button onClick={model.doSomething}>Do something</Button>
        {/* <Button onClick={() => {const n = model.findNodesByTags(CVSX_VOLUME_TAG); console.log('NODES'); console.log(n)}}>Find nodes by Tag</Button> */}
    </>;
}
