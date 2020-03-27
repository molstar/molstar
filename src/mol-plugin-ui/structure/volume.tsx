/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { InitVolumeStreaming } from '../../mol-plugin/behavior/dynamic/volume-streaming/transformers';
import { CollapsableControls, CollapsableState } from '../base';
import { ApplyActionControl } from '../state/apply-action';
import { UpdateTransformControl } from '../state/update-transform';
import { BindingsHelp } from '../viewport/help';
import { ExpandGroup } from '../controls/common';
import { StructureHierarchyManager } from '../../mol-plugin-state/manager/structure/hierarchy';
import { FocusLoci } from '../../mol-plugin/behavior/dynamic/representation';

interface VolumeStreamingControlState extends CollapsableState {
    isBusy: boolean
}

export class VolumeStreamingControls extends CollapsableControls<{}, VolumeStreamingControlState> {
    protected defaultState(): VolumeStreamingControlState {
        return {
            header: 'Volume Streaming',
            isCollapsed: false,
            isBusy: false,
            isHidden: true
        };
    }

    componentDidMount() {
        // TODO: do not hide this but instead show some help text??
        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, () => this.setState({ isHidden: !this.canEnable() }));
        this.subscribe(this.plugin.behaviors.state.isBusy, v => this.setState({ isBusy: v }));

        this.subscribe(this.plugin.managers.structure.hierarchy.behaviors.selection, c => this.setState({
            description: StructureHierarchyManager.getSelectedStructuresDescription(this.plugin)
        }));
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
        return <ApplyActionControl state={pivot.cell.parent} action={InitVolumeStreaming} initiallyCollapsed={true} nodeRef={pivot.cell.transform.ref} simpleApply={{ header: 'Enable', icon: 'check' }} />;
    }

    renderParams() {
        const pivot = this.pivot;
        const bindings = pivot.volumeStreaming?.cell.transform.params?.entry.params.view.name === 'selection-box' && this.plugin.state.behaviors.cells.get(FocusLoci.id)?.params?.values?.bindings;
        return <>
            <UpdateTransformControl state={pivot.cell.parent} transform={pivot.volumeStreaming!.cell.transform} customHeader='none' noMargin autoHideApply />
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