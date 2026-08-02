/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { ParticleHierarchyManager } from '../../mol-plugin-state/manager/particles/hierarchy';
import { ParticleListRef, ParticleRepresentationRef } from '../../mol-plugin-state/manager/particles/hierarchy-state';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { PluginCommands } from '../../mol-plugin/commands';
import { State } from '../../mol-state';
import { Color } from '../../mol-util/color';
import { ParamDefinition } from '../../mol-util/param-definition';
import { CollapsableControls, CollapsableState, PurePluginUIComponent } from '../base';
import { ActionMenu } from '../controls/action-menu';
import { CombinedColorControl } from '../controls/color';
import { Button, ControlGroup, IconButton } from '../controls/common';
import { AddSvg, CloseSvg, DeleteOutlinedSvg, MoreHorizSvg, ScatterPlotSvg, VisibilityOffOutlinedSvg, VisibilityOutlinedSvg } from '../controls/icons';
import { ParameterControls, ParamOnChange } from '../controls/parameters';
import { UpdateTransformControl } from '../state/update-transform';

interface ParticleSourceControlState extends CollapsableState {
    isBusy: boolean,
    show?: 'hierarchy' | 'add-repr'
}

export class ParticleSourceControls extends CollapsableControls<{}, ParticleSourceControlState> {
    protected defaultState(): ParticleSourceControlState {
        return {
            header: 'Particles',
            isCollapsed: false,
            isBusy: false,
            isHidden: true,
            brand: { accent: 'orange', svg: ScatterPlotSvg }
        };
    }

    componentDidMount() {
        this.subscribe(this.plugin.managers.particles.hierarchy.behaviors.selection, sel => {
            this.setState({ isHidden: sel.hierarchy.lists.length === 0 });
        });
        this.subscribe(this.plugin.behaviors.state.isBusy, v => {
            this.setState({ isBusy: v });
        });
    }

    private item = (ref: ParticleListRef) => {
        const selected = this.plugin.managers.particles.hierarchy.selection;

        const label = ref.cell.obj?.label || 'Particles';
        const item: ActionMenu.Item = { kind: 'item', label, selected: selected === ref, value: ref };
        return item;
    };

    get hierarchyItems() {
        const { current } = this.plugin.managers.particles.hierarchy;
        return current.lists.map(this.item);
    }

    get addActions(): ActionMenu.Items {
        const mng = this.plugin.managers.particles.hierarchy;
        const current = mng.selection;

        return ParticleHierarchyManager.getRepresentationTypes(this.plugin, current)
            .map(t => ActionMenu.Item(t[1], () => mng.addRepresentation(current!, t[0])));
    }

    get isEmpty() {
        return this.plugin.managers.particles.hierarchy.current.lists.length === 0;
    }

    get label() {
        const selected = this.plugin.managers.particles.hierarchy.selection;
        if (!selected) return 'Nothing Selected';
        return selected.cell.obj?.label || 'Particles';
    }

    selectCurrent: ActionMenu.OnSelect = (item) => {
        this.toggleHierarchy();
        if (!item) return;
        this.plugin.managers.particles.hierarchy.setCurrent(item.value as ParticleListRef);
    };

    selectAdd: ActionMenu.OnSelect = (item) => {
        if (!item) return;
        this.setState({ show: void 0 });
        (item.value as any)();
    };

    toggleHierarchy = () => this.setState({ show: this.state.show !== 'hierarchy' ? 'hierarchy' : void 0 });
    toggleAddRepr = () => this.setState({ show: this.state.show !== 'add-repr' ? 'add-repr' : void 0 });
    toggleVisibility = () => {
        const mng = this.plugin.managers.particles.hierarchy;
        const { current } = mng;
        const globalVisibility = !current.lists[0]?.representations[0]?.cell.state.isHidden;
        mng.toggleVisibility(current.lists.flatMap(l => l.representations), globalVisibility ? 'hide' : 'show');
    };

    private updateFrameQueueParams: any = void 0;
    private isUpdatingFrame = false;

    private async _updateFrameIndex() {
        if (!this.updateFrameQueueParams || this.isUpdatingFrame) return;
        const params = this.updateFrameQueueParams;
        this.updateFrameQueueParams = void 0;

        try {
            this.isUpdatingFrame = true;
            const list = this.plugin.managers.particles.hierarchy.selection;
            if (!list) return;
            await this.plugin.state.updateTransform(this.plugin.state.data, list.cell.transform.ref, params, 'Particle Frame Index');
        } finally {
            this.isUpdatingFrame = false;
            this._updateFrameIndex();
        }
    }

    updateFrameIndex = (params: any) => {
        this.updateFrameQueueParams = params;
        this._updateFrameIndex();
    };

    get frameIndex() {
        const list = this.plugin.managers.particles.hierarchy.selection;
        if (!list) return null;
        // only frames extracted from a multi-frame trajectory expose a frame index to scrub through
        if (list.cell.transform.transformer !== StateTransforms.Particles.ParticleListFromTrajectory) return null;

        const params = list.cell.params?.definition;
        if (!params) return null;

        return <ParameterControls params={params} values={list.cell.params?.values} onChangeValues={this.updateFrameIndex} isDisabled={this.state.isBusy} />;
    }

    renderControls() {
        const disabled = this.state.isBusy || this.isEmpty;
        const label = this.label;
        const selected = this.plugin.managers.particles.hierarchy.selection;
        const { current } = this.plugin.managers.particles.hierarchy;

        return <>
            <div className='msp-flex-row' style={{ marginTop: '1px' }}>
                <Button noOverflow flex onClick={this.toggleHierarchy} disabled={disabled} title={label}>{label}</Button>
                {!this.isEmpty && selected && <IconButton svg={AddSvg} onClick={this.toggleAddRepr} title='Add a representation for the current particle list.' toggleState={this.state.show === 'add-repr'} disabled={disabled} />}
                {!this.isEmpty && <IconButton svg={VisibilityOutlinedSvg} onClick={this.toggleVisibility} toggleState={false} title='Toggle visibility of all particle lists.' disabled={disabled} />}
            </div>
            {this.state.show === 'hierarchy' && <ActionMenu items={this.hierarchyItems} onSelect={this.selectCurrent} />}
            {this.state.show === 'add-repr' && <ActionMenu items={this.addActions} onSelect={this.selectAdd} />}
            {this.frameIndex}

            {current.lists.length > 0 && <div style={{ marginTop: '6px' }}>
                {current.lists.map((list) => <ParticleEntryControls list={list} key={list.cell.transform.ref} />)}
            </div>}
        </>;
    }
}

function ParticleEntryControls({ list }: { list: ParticleListRef }) {
    return <>
        <div className='msp-control-group-header' style={{ marginTop: '1px' }}>
            <div><b>{list.cell.obj?.label ?? 'n/a'}</b></div>
        </div>
        {list.representations.map(r => <ParticleRepresentationControls key={r.cell.transform.ref} list={list} representation={r} />)}
    </>;
}

type ParticleRepresentationEntryActions = 'update' | 'select-color'

class ParticleRepresentationControls extends PurePluginUIComponent<{ list: ParticleListRef, representation: ParticleRepresentationRef }, { action?: ParticleRepresentationEntryActions }> {
    state = { action: void 0 as ParticleRepresentationEntryActions | undefined };

    componentDidMount() {
        this.subscribe(this.plugin.state.events.cell.stateUpdated, e => {
            if (State.ObjectEvent.isCell(e, this.props.representation.cell)) this.forceUpdate();
        });
    }

    remove = () => this.plugin.managers.particles.hierarchy.remove([this.props.representation], true);

    toggleVisible = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        e.currentTarget.blur();
        this.plugin.managers.particles.hierarchy.toggleVisibility([this.props.representation]);
    };

    toggleColor = () => {
        this.setState({ action: this.state.action === 'select-color' ? undefined : 'select-color' });
    };

    toggleUpdate = () => this.setState({ action: this.state.action === 'update' ? void 0 : 'update' });

    highlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        if (!this.props.representation.cell.parent) return;
        PluginCommands.Interactivity.Object.Highlight(this.plugin, { state: this.props.representation.cell.parent!, ref: this.props.representation.cell.transform.ref });
    };

    clearHighlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.Interactivity.ClearHighlights(this.plugin);
    };

    focus = () => {
        const repr = this.props.representation;
        const lociList = repr.cell.obj?.data.repr.getAllLoci();
        if (repr.cell.state.isHidden) this.plugin.managers.particles.hierarchy.toggleVisibility([this.props.representation], 'show');
        if (lociList) this.plugin.managers.camera.focusLoci(lociList, { extraRadius: 1 });
    };

    private get color() {
        const repr = this.props.representation.cell;
        const isUniform = repr.transform.params?.colorTheme.name === 'uniform';
        if (!isUniform) return void 0;
        return repr.transform.params?.colorTheme.params.value;
    }

    updateColor: ParamOnChange = ({ value }) => {
        const t = this.props.representation.cell.transform;
        return this.plugin.build().to(t.ref).update({
            ...t.params,
            colorTheme: {
                name: 'uniform',
                params: { value }
            },
        }).commit();
    };

    render() {
        const repr = this.props.representation.cell;
        const color = this.color;

        return <>
            <div className='msp-flex-row'>
                {color !== void 0 && <Button style={{ backgroundColor: Color.toStyle(color), minWidth: 32, width: 32 }} onClick={this.toggleColor} />}
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
            {this.state.action === 'select-color' && color !== void 0 && <div style={{ marginBottom: '6px', marginTop: 1 }} className='msp-accent-offset'>
                <ControlGroup header='Select Color' initialExpanded={true} hideExpander={true} hideOffset={true} onHeaderClick={this.toggleColor}
                    topRightIcon={CloseSvg} noTopMargin childrenClassName='msp-viewport-controls-panel-controls'>
                    <CombinedColorControl param={ParticleColorParam} value={this.color} onChange={this.updateColor} name='color' hideNameRow />
                </ControlGroup>
            </div>}
        </>;
    }
}

const ParticleColorParam = ParamDefinition.Color(Color(0x121212));
