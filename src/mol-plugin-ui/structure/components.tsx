/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { CollapsableControls, CollapsableState, PurePluginUIComponent } from '../base';
import { StructureHierarchyManager } from '../../mol-plugin-state/manager/structure';
import { StructureComponentRef } from '../../mol-plugin-state/manager/structure/hierarchy';

interface StructureComponentControlState extends CollapsableState {
    isDisabled: boolean
}

export class StructureComponentControls extends CollapsableControls<{}, StructureComponentControlState> {
    protected defaultState(): StructureComponentControlState {
        return { header: 'Structure Components', isCollapsed: false, isDisabled: false };
    }

    get currentModels() {
        return this.plugin.managers.structureHierarchy.behaviors.currentModels;
    }

    componentDidMount() {
        this.subscribe(this.currentModels, () => this.forceUpdate());
    }

    renderControls() {
        const components = StructureHierarchyManager.getCommonComponentPivots(this.currentModels.value)
        return <>
            {components.map((c, i) => <StructureComponentEntry key={i} component={c} />)}
        </>;
    }
}

class StructureComponentEntry extends PurePluginUIComponent<{ component: StructureComponentRef }> {
    render() {
        return <></>;
    }
}