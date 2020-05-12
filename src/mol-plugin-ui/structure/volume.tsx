/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { StructureHierarchyManager } from '../../mol-plugin-state/manager/structure/hierarchy';
import { VolumeHierarchyManager } from '../../mol-plugin-state/manager/volume/hierarchy';
import { VolumeRef, VolumeRepresentationRef } from '../../mol-plugin-state/manager/volume/hierarchy-state';
import { FocusLoci } from '../../mol-plugin/behavior/dynamic/representation';
import { VolumeStreaming } from '../../mol-plugin/behavior/dynamic/volume-streaming/behavior';
import { InitVolumeStreaming } from '../../mol-plugin/behavior/dynamic/volume-streaming/transformers';
import { State, StateSelection, StateTransform } from '../../mol-state';
import { CollapsableControls, CollapsableState, PurePluginUIComponent } from '../base';
import { ActionMenu } from '../controls/action-menu';
import { Button, ExpandGroup, IconButton } from '../controls/common';
import { ApplyActionControl } from '../state/apply-action';
import { UpdateTransformControl } from '../state/update-transform';
import { BindingsHelp } from '../viewport/help';
import { PluginCommands } from '../../mol-plugin/commands';
import { BlurOnSvg, ErrorSvg, CheckSvg, AddSvg, VisibilityOffOutlinedSvg, VisibilityOutlinedSvg, DeleteOutlinedSvg, MoreHorizSvg } from '../controls/icons';

interface VolumeStreamingControlState extends CollapsableState {
    isBusy: boolean
}

export class VolumeStreamingControls extends CollapsableControls<{}, VolumeStreamingControlState> {
    protected defaultState(): VolumeStreamingControlState {
        return {
            header: 'Volume Streaming',
            isCollapsed: false,
            isBusy: false,
            isHidden: true,
            brand: { accent: 'cyan', svg: BlurOnSvg }
        };
    }

    componentDidMount() {
        // TODO: do not hide this but instead show some help text??
        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, () => {
            this.setState({
                isHidden: !this.canEnable(),
                description: StructureHierarchyManager.getSelectedStructuresDescription(this.plugin)
            });
        });
        this.subscribe(this.plugin.state.events.cell.stateUpdated, e => {
            if (StateTransform.hasTag(e.cell.transform, VolumeStreaming.RootTag)) this.forceUpdate();
        });
        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isBusy: v });
        });
    }

    get pivot() {
        return this.plugin.managers.structure.hierarchy.selection.structures[0];
    }

    canEnable() {
        const { selection } = this.plugin.managers.structure.hierarchy;
        if (selection.structures.length !== 1) return false;
        const pivot = this.pivot.cell;
        if (!pivot.obj) return false;
        return !!InitVolumeStreaming.definition.isApplicable?.(pivot.obj, pivot.transform, this.plugin);
    }

    renderEnable() {
        const pivot = this.pivot;

        if (!pivot.cell.parent) return null;

        const root = StateSelection.findTagInSubtree(pivot.cell.parent.tree, this.pivot.cell.transform.ref, VolumeStreaming.RootTag);
        const rootCell = root && pivot.cell.parent.cells.get(root);

        const simpleApply = rootCell && rootCell.status === 'error'
            ? { header: 'Error enabling', icon: ErrorSvg, title: rootCell.errorText }
            : rootCell && rootCell.obj?.data.entries.length === 0
                ? { header: 'Error enabling', icon: ErrorSvg, title: 'No entry for streaming found' }
                : { header: 'Enable', icon: CheckSvg, title: 'Enable' };

        return <ApplyActionControl state={pivot.cell.parent} action={InitVolumeStreaming} initiallyCollapsed={true} nodeRef={pivot.cell.transform.ref} simpleApply={simpleApply} />;
    }

    renderParams() {
        const pivot = this.pivot;
        if (!pivot.cell.parent) return null;
        const bindings = pivot.volumeStreaming?.cell.transform.params?.entry.params.view.name === 'selection-box' && this.plugin.state.behaviors.cells.get(FocusLoci.id)?.params?.values?.bindings;
        return <>
            <UpdateTransformControl state={pivot.cell.parent} transform={pivot.volumeStreaming!.cell.transform} customHeader='none' noMargin />
            {bindings && <ExpandGroup header='Controls Help'>
                <BindingsHelp bindings={bindings} />
            </ExpandGroup>}
        </>;
    }

    renderControls() {
        const pivot = this.pivot;
        if (!pivot) return null;
        if (!pivot.volumeStreaming) return this.renderEnable();
        return this.renderParams();
    }
}

interface VolumeSourceControlState extends CollapsableState {
    isBusy: boolean,
    show?: 'hierarchy' | 'add-repr'
}

export class VolumeSourceControls extends CollapsableControls<{}, VolumeSourceControlState> {
    protected defaultState(): VolumeSourceControlState {
        return {
            header: 'Volume',
            isCollapsed: false,
            isBusy: false,
            isHidden: true,
            brand: { accent: 'purple', svg: BlurOnSvg }
        };
    }

    componentDidMount() {
        this.subscribe(this.plugin.managers.volume.hierarchy.behaviors.selection, sel => {
            this.setState({ isHidden: sel.hierarchy.volumes.length === 0 });
        });
        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isBusy: v });
        });
    }

    private item = (ref: VolumeRef) => {
        const selected = this.plugin.managers.volume.hierarchy.selection;

        const label = ref.cell.obj?.label || 'Volume';
        const item: ActionMenu.Item = { kind: 'item', label: label || ref.kind, selected: selected === ref, value: ref };
        return item;
    }

    get hierarchyItems() {
        const mng = this.plugin.managers.volume.hierarchy;
        const { current } = mng;
        const ret: ActionMenu.Items = [];
        for (let ref of current.volumes) {
            ret.push(this.item(ref));
        }
        return ret;
    }

    get addActions(): ActionMenu.Items {
        const mng = this.plugin.managers.volume.hierarchy;
        const current = mng.selection;

        const ret: ActionMenu.Items = [
            ...VolumeHierarchyManager.getRepresentationTypes(this.plugin, current)
                .map(t => ActionMenu.Item(t[1], () => mng.addRepresentation(current!, t[0])))
        ];

        return ret;
    }

    get isEmpty() {
        const { volumes } = this.plugin.managers.volume.hierarchy.current;
        return volumes.length === 0;
    }

    get label() {
        const selected = this.plugin.managers.volume.hierarchy.selection;
        if (!selected) return 'Nothing Selected';
        return selected?.cell.obj?.label || 'Volume';
    }

    selectCurrent: ActionMenu.OnSelect = (item) => {
        this.toggleHierarchy();
        if (!item) return;
        this.plugin.managers.volume.hierarchy.setCurrent(item.value as VolumeRef);
    }

    selectAdd: ActionMenu.OnSelect = (item) => {
        if (!item) return;
        this.setState({ show: void 0 });
        (item.value as any)();
    }

    toggleHierarchy = () => this.setState({ show: this.state.show !== 'hierarchy' ? 'hierarchy' : void 0 });
    toggleAddRepr = () => this.setState({ show: this.state.show !== 'add-repr' ? 'add-repr' : void 0 });

    renderControls() {
        const disabled = this.state.isBusy || this.isEmpty;
        const label = this.label;

        const selected = this.plugin.managers.volume.hierarchy.selection;

        return <>
            <div className='msp-flex-row' style={{ marginTop: '1px' }}>
                <Button noOverflow flex onClick={this.toggleHierarchy} disabled={disabled} title={label}>{label}</Button>
                {!this.isEmpty && <IconButton svg={AddSvg} onClick={this.toggleAddRepr} title='Apply a structure presets to the current hierarchy.' toggleState={this.state.show === 'add-repr'} disabled={disabled} />}
            </div>
            {this.state.show === 'hierarchy' && <ActionMenu items={this.hierarchyItems} onSelect={this.selectCurrent} />}
            {this.state.show === 'add-repr' && <ActionMenu items={this.addActions} onSelect={this.selectAdd} />}

            {selected && selected.representations.length > 0 && <div style={{ marginTop: '6px' }}>
                {selected.representations.map(r => <VolumeRepresentationControls key={r.cell.transform.ref} representation={r} />)}
            </div>}
        </>;
    }
}

type VolumeRepresentationEntryActions = 'update'

class VolumeRepresentationControls extends PurePluginUIComponent<{ representation: VolumeRepresentationRef }, { action?: VolumeRepresentationEntryActions }> {
    state = { action: void 0 as VolumeRepresentationEntryActions | undefined }

    componentDidMount() {
        this.subscribe(this.plugin.state.events.cell.stateUpdated, e => {
            if (State.ObjectEvent.isCell(e, this.props.representation.cell)) this.forceUpdate();
        });
    }

    remove = () => this.plugin.managers.volume.hierarchy.remove([ this.props.representation ], true);

    toggleVisible = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        e.currentTarget.blur();
        this.plugin.managers.volume.hierarchy.toggleVisibility([ this.props.representation ]);
    }

    toggleUpdate = () => this.setState({ action: this.state.action === 'update' ? void 0 : 'update' });

    highlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        if (!this.props.representation.cell.parent) return;
        PluginCommands.Interactivity.Object.Highlight(this.plugin, { state: this.props.representation.cell.parent!, ref: this.props.representation.cell.transform.ref });
    }

    clearHighlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.Interactivity.ClearHighlights(this.plugin);
    }

    focus = () => {
        const repr = this.props.representation;
        const objects = this.props.representation.cell.obj?.data.repr.renderObjects;
        if (repr.cell.state.isHidden) this.plugin.managers.volume.hierarchy.toggleVisibility([this.props.representation], 'show');
        this.plugin.managers.camera.focusRenderObjects(objects, { extraRadius: 1 });
    }

    render() {
        const repr = this.props.representation.cell;
        return <>
            <div className='msp-flex-row'>
                <Button noOverflow className='msp-control-button-label' title={`${repr.obj?.label}. Click to focus.`} onClick={this.focus} onMouseEnter={this.highlight} onMouseLeave={this.clearHighlight} style={{ textAlign: 'left' }}>
                    {repr.obj?.label}
                    <small className='msp-25-lower-contrast-text' style={{ float: 'right' }}>{repr.obj?.description}</small>
                </Button>
                <IconButton svg={repr.state.isHidden ? VisibilityOffOutlinedSvg : VisibilityOutlinedSvg} toggleState={false} onClick={this.toggleVisible} title={`${repr.state.isHidden ? 'Show' : 'Hide'} component`} small className='msp-form-control' flex />
                <IconButton svg={DeleteOutlinedSvg} onClick={this.remove} title='Remove' small />
                <IconButton svg={MoreHorizSvg} onClick={this.toggleUpdate} title='Actions' toggleState={this.state.action === 'update'} />
            </div>
            {this.state.action === 'update' && !!repr.parent && <div style={{ marginBottom: '6px' }} className='msp-accent-offset'>
                <UpdateTransformControl state={repr.parent} transform={repr.transform} customHeader='none' noMargin />
            </div>}
        </>;
    }
}