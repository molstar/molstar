/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { PluginCommands } from 'mol-plugin/command';
import { UpdateTrajectory } from 'mol-plugin/state/actions/structure';
import { PluginUIComponent } from './base';
import { LociLabelEntry } from 'mol-plugin/util/loci-label-manager';
import { IconButton } from './controls/common';
import { PluginStateObject } from 'mol-plugin/state/objects';
import { StateTransforms } from 'mol-plugin/state/transforms';
import { StateTransformer } from 'mol-state';
import { ModelFromTrajectory } from 'mol-plugin/state/transforms/model';

export class TrajectoryControls extends PluginUIComponent<{}, { show: boolean, label: string }> {
    state = { show: false, label: '' }

    private update = () => {
        const state = this.plugin.state.dataState;

        const models = state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Model)
            .filter(c => c.transform.transformer === StateTransforms.Model.ModelFromTrajectory));

        if (models.length === 0) {
            this.setState({ show: false });
            return;
        }

        let label = '', count = 0, parents = new Set<string>();
        for (const m of models) {
            if (!m.sourceRef) continue;
            const parent = state.cells.get(m.sourceRef)!.obj as PluginStateObject.Molecule.Trajectory;

            if (!parent) continue;
            if (parent.data.length > 1) {
                if (parents.has(m.sourceRef)) {
                    // do not show the controls if there are 2 models of the same trajectory present
                    this.setState({ show: false });
                }

                parents.add(m.sourceRef);
                count++;
                if (!label) {
                    const idx = (m.transform.params! as StateTransformer.Params<ModelFromTrajectory>).modelIndex;
                    label = `Model ${idx + 1} / ${parent.data.length}`;
                }
            }
        }

        if (count > 1) label = '';
        this.setState({ show: count > 0, label });
    }

    componentDidMount() {
        this.subscribe(this.plugin.state.dataState.events.changed, this.update);
    }

    reset = () => PluginCommands.State.ApplyAction.dispatch(this.plugin, {
        state: this.plugin.state.dataState,
        action: UpdateTrajectory.create({ action: 'reset' })
    });

    prev = () => PluginCommands.State.ApplyAction.dispatch(this.plugin, {
        state: this.plugin.state.dataState,
        action: UpdateTrajectory.create({ action: 'advance', by: -1 })
    });

    next = () => PluginCommands.State.ApplyAction.dispatch(this.plugin, {
        state: this.plugin.state.dataState,
        action: UpdateTrajectory.create({ action: 'advance', by: 1 })
    });

    render() {
        if (!this.state.show) return null;

        return <div className='msp-traj-controls'>
            <IconButton icon='model-first' title='First Model' onClick={this.reset} />
            <IconButton icon='model-prev' title='Previous Model' onClick={this.prev} />
            <IconButton icon='model-next' title='Next Model' onClick={this.next} />
            { !!this.state.label && <span>{this.state.label}</span> }
        </div>;
    }
}

export class StateSnapshotViewportControls extends PluginUIComponent<{}, { isBusy: boolean }> {
    state = { isBusy: false }

    componentDidMount() {
        // TODO: this needs to be diabled when the state is updating!
        this.subscribe(this.plugin.state.snapshots.events.changed, () => this.forceUpdate());
        this.subscribe(this.plugin.events.state.isUpdating, (isBusy) => this.setState({ isBusy }));
    }

    async update(id: string) {
        this.setState({ isBusy: true });
        await PluginCommands.State.Snapshots.Apply.dispatch(this.plugin, { id });
        this.setState({ isBusy: false });
    }

    change = (e: React.ChangeEvent<HTMLSelectElement>) => {
        if (e.target.value === 'none') return;
        this.update(e.target.value);
    }

    prev =  () => {
        const s = this.plugin.state.snapshots;
        const id = s.getNextId(s.state.current, -1);
        if (id) this.update(id);
    }

    next =  () => {
        const s = this.plugin.state.snapshots;
        const id = s.getNextId(s.state.current, 1);
        if (id) this.update(id);
    }

    render() {
        const snapshots = this.plugin.state.snapshots;
        const count = snapshots.state.entries.size;

        if (count < 2) {
            return null;
        }

        const current = snapshots.state.current;

        return <div className='msp-state-snapshot-viewport-controls'>
            <select className='msp-form-control' value={current || 'none'} onChange={this.change} disabled={this.state.isBusy}>
                {!current && <option key='none' value='none'></option>}
                {snapshots.state.entries.valueSeq().map((e, i) => <option key={e!.snapshot.id} value={e!.snapshot.id}>{`[${i! + 1}/${count}]`} {e!.name || new Date(e!.timestamp).toLocaleString()}</option>)}
            </select>
            <IconButton icon='model-prev' title='Previous State' onClick={this.prev} disabled={this.state.isBusy} />
            <IconButton icon='model-next' title='Next State' onClick={this.next} disabled={this.state.isBusy} />
        </div>;
    }
}

export class LociLabelControl extends PluginUIComponent<{}, { entries: ReadonlyArray<LociLabelEntry> }> {
    state = { entries: [] }

    componentDidMount() {
        this.subscribe(this.plugin.behaviors.labels.highlight, e => this.setState({ entries: e.entries }));
    }

    render() {
        if (this.state.entries.length === 0) return null;

        return <div className='msp-highlight-info'>
            {this.state.entries.map((e, i) => <div key={'' + i}>{e}</div>)}
        </div>;
    }
}