/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { LoadCellPackModel } from '../../../extensions/cellpack/model';
import { LoadPetworldModel } from '../../../extensions/petworld/model';
import { PluginUIComponent } from '../../../mol-plugin-ui/base';
import { ApplyActionControl } from '../../../mol-plugin-ui/state/apply-action';

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