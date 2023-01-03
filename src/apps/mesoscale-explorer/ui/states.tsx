/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { LoadCellPackModel } from '../../../extensions/cellpack/model';
import { LoadPetworldModel } from '../../../extensions/petworld/model';
import { PluginUIComponent } from '../../../mol-plugin-ui/base';
import { ExpandGroup } from '../../../mol-plugin-ui/controls/common';
import { ApplyActionControl } from '../../../mol-plugin-ui/state/apply-action';
import { LocalStateSnapshotList, LocalStateSnapshotParams, LocalStateSnapshots, StateExportImportControls } from '../../../mol-plugin-ui/state/snapshots';

export class LoaderControls extends PluginUIComponent {
    componentDidMount() {

    }

    render() {
        return <div style={{ margin: '5px' }}>
            <ApplyActionControl state={this.plugin.state.data} action={LoadPetworldModel} nodeRef={this.plugin.state.data.tree.root.ref} />
            <ApplyActionControl state={this.plugin.state.data} action={LoadCellPackModel} nodeRef={this.plugin.state.data.tree.root.ref} />
        </div>;
    }
}

export class StateControls extends PluginUIComponent<{}> {
    render() {
        return <div style={{ margin: '5px' }}>
            <StateExportImportControls />

            <div style={{ marginBottom: '10px' }}>
                <LocalStateSnapshotList />
            </div>
            <div style={{ marginBottom: '10px' }}>
                <LocalStateSnapshots />
            </div>

            <div style={{ marginBottom: '10px' }}>
                <ExpandGroup header='Save Options' initiallyExpanded={false}>
                    <LocalStateSnapshotParams />
                </ExpandGroup>
            </div>
        </div>;
    }
}