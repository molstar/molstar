/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { PluginCommands } from '../mol-plugin/commands';
import { UpdateTrajectory } from '../mol-plugin-state/actions/structure';
import { PluginUIComponent } from './base';
import { LociLabelEntry } from '../mol-plugin-state/manager/loci-label';
import { IconButton } from './controls/common';
import { PluginStateObject } from '../mol-plugin-state/objects';
import { StateTransforms } from '../mol-plugin-state/transforms';
import { StateTransformer } from '../mol-state';
import { ModelFromTrajectory } from '../mol-plugin-state/transforms/model';
import { AnimationControls } from './state/animation';
import { StructureSelectionControls } from './structure/selection';
import { StructureMeasurementsControls } from './structure/measurements';
import { Icon } from './controls/icons';
import { StructureComponentControls } from './structure/components';
import { StructureSourceControls } from './structure/source';
import { VolumeStreamingControls } from './structure/volume';
import { ObjectControls } from './structure/objects';

export class TrajectoryViewportControls extends PluginUIComponent<{}, { show: boolean, label: string }> {
    state = { show: false, label: '' }

    private update = () => {
        const state = this.plugin.state.data;

        const models = state.selectQ(q => q.ofTransformer(StateTransforms.Model.ModelFromTrajectory));

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
                    return;
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
        this.subscribe(this.plugin.state.data.events.changed, this.update);
        this.subscribe(this.plugin.behaviors.state.isAnimating, this.update);
    }

    reset = () => PluginCommands.State.ApplyAction(this.plugin, {
        state: this.plugin.state.data,
        action: UpdateTrajectory.create({ action: 'reset' })
    });

    prev = () => PluginCommands.State.ApplyAction(this.plugin, {
        state: this.plugin.state.data,
        action: UpdateTrajectory.create({ action: 'advance', by: -1 })
    });

    next = () => PluginCommands.State.ApplyAction(this.plugin, {
        state: this.plugin.state.data,
        action: UpdateTrajectory.create({ action: 'advance', by: 1 })
    });

    // stopAnimation = () => {
    //     this.plugin.state.animation.stop();
    // }

    // playAnimation = () => {
    //     const anim = this.plugin.state.animation;
    //     if (anim.state.params.current === AnimateModelIndex.name) {
    //         anim.start();
    //     } else {
    //         anim.play(AnimateModelIndex, ParamDefinition.getDefaultValues(AnimateModelIndex.params(this.plugin) as any as ParamDefinition.Params))
    //     }
    // }

    render() {
        const isAnimating = this.plugin.behaviors.state.isAnimating.value;

        if (!this.state.show || (isAnimating && !this.state.label)) return null;

        return <div className='msp-traj-controls'>
            {/* <IconButton icon={isAnimating ? 'stop' : 'play'} title={isAnimating ? 'Stop' : 'Play'} onClick={isAnimating ? this.stopAnimation : this.playAnimation} /> */}
            {!isAnimating && <IconButton icon='model-first' title='First Model' onClick={this.reset} disabled={isAnimating} />}
            {!isAnimating && <IconButton icon='model-prev' title='Previous Model' onClick={this.prev} disabled={isAnimating} />}
            {!isAnimating && <IconButton icon='model-next' title='Next Model' onClick={this.next} disabled={isAnimating} />}
            {!!this.state.label && <span>{this.state.label}</span> }
        </div>;
    }
}

export class StateSnapshotViewportControls extends PluginUIComponent<{}, { isBusy: boolean, show: boolean }> {
    state = { isBusy: false, show: true }

    componentDidMount() {
        // TODO: this needs to be diabled when the state is updating!
        this.subscribe(this.plugin.state.snapshots.events.changed, () => this.forceUpdate());
        this.subscribe(this.plugin.behaviors.state.isUpdating, isBusy => this.setState({ isBusy }));
        this.subscribe(this.plugin.behaviors.state.isAnimating, isBusy => this.setState({ isBusy }))

        window.addEventListener('keyup', this.keyUp, false);
    }

    componentWillUnmount() {
        super.componentWillUnmount();
        window.removeEventListener('keyup', this.keyUp, false);
    }

    keyUp = (e: KeyboardEvent) => {
        if (!e.ctrlKey || this.state.isBusy || e.target !== document.body) return;
        const snapshots = this.plugin.state.snapshots;
        if (e.keyCode === 37) { // left
            if (snapshots.state.isPlaying) snapshots.stop();
            this.prev();
        } else if (e.keyCode === 38) { // up
            if (snapshots.state.isPlaying) snapshots.stop();
            if (snapshots.state.entries.size === 0) return;
            const e = snapshots.state.entries.get(0);
            this.update(e.snapshot.id);
        } else if (e.keyCode === 39) { // right
            if (snapshots.state.isPlaying) snapshots.stop();
            this.next();
        } else if (e.keyCode === 40) { // down
            if (snapshots.state.isPlaying) snapshots.stop();
            if (snapshots.state.entries.size === 0) return;
            const e = snapshots.state.entries.get(snapshots.state.entries.size - 1);
            this.update(e.snapshot.id);
        }
    };

    async update(id: string) {
        this.setState({ isBusy: true });
        await PluginCommands.State.Snapshots.Apply(this.plugin, { id });
        this.setState({ isBusy: false });
    }

    change = (e: React.ChangeEvent<HTMLSelectElement>) => {
        if (e.target.value === 'none') return;
        this.update(e.target.value);
    }

    prev = () => {
        const s = this.plugin.state.snapshots;
        const id = s.getNextId(s.state.current, -1);
        if (id) this.update(id);
    }

    next = () => {
        const s = this.plugin.state.snapshots;
        const id = s.getNextId(s.state.current, 1);
        if (id) this.update(id);
    }

    togglePlay = () => {
        this.plugin.state.snapshots.togglePlay();
    }

    render() {
        const snapshots = this.plugin.state.snapshots;
        const count = snapshots.state.entries.size;

        if (count < 2 || !this.state.show) {
            return null;
        }

        const current = snapshots.state.current;
        const isPlaying = snapshots.state.isPlaying;

        return <div className='msp-state-snapshot-viewport-controls'>
            <select className='msp-form-control' value={current || 'none'} onChange={this.change} disabled={this.state.isBusy || isPlaying}>
                {!current && <option key='none' value='none'></option>}
                {snapshots.state.entries.valueSeq().map((e, i) => <option key={e!.snapshot.id} value={e!.snapshot.id}>{`[${i! + 1}/${count}]`} {e!.name || new Date(e!.timestamp).toLocaleString()}</option>)}
            </select>
            <IconButton icon={isPlaying ? 'stop' : 'play'} title={isPlaying ? 'Pause' : 'Cycle States'} onClick={this.togglePlay}
                disabled={isPlaying ? false : this.state.isBusy} />
            {!isPlaying && <>
                <IconButton icon='left-open' title='Previous State' onClick={this.prev} disabled={this.state.isBusy || isPlaying} />
                <IconButton icon='right-open' title='Next State' onClick={this.next} disabled={this.state.isBusy || isPlaying} />
            </>}
        </div>;
    }
}

export class AnimationViewportControls extends PluginUIComponent<{}, { isEmpty: boolean, isExpanded: boolean, isUpdating: boolean, isAnimating: boolean, isPlaying: boolean }> {
    state = { isEmpty: true, isExpanded: false, isUpdating: false, isAnimating: false, isPlaying: false };

    componentDidMount() {
        this.subscribe(this.plugin.state.snapshots.events.changed, () => {
            if (this.plugin.state.snapshots.state.isPlaying) this.setState({ isPlaying: true, isExpanded: false });
            else this.setState({ isPlaying: false });
        });
        this.subscribe(this.plugin.behaviors.state.isUpdating, isUpdating => {
            if (isUpdating) this.setState({ isUpdating: true, isExpanded: false, isEmpty: this.plugin.state.data.tree.transforms.size < 2 });
            else this.setState({ isUpdating: false, isEmpty: this.plugin.state.data.tree.transforms.size < 2 });
        });
        this.subscribe(this.plugin.behaviors.state.isAnimating, isAnimating => {
            if (isAnimating) this.setState({ isAnimating: true, isExpanded: false });
            else this.setState({ isAnimating: false });
        });
    }
    toggleExpanded = () => this.setState({ isExpanded: !this.state.isExpanded });
    stop = () => {
        this.plugin.state.animation.stop();
        this.plugin.state.snapshots.stop();
    }

    render() {
        // if (!this.state.show) return null;
        const isPlaying = this.plugin.state.snapshots.state.isPlaying;
        if (isPlaying || this.state.isEmpty) return null;

        const isAnimating = this.state.isAnimating;

        return <div className='msp-animation-viewport-controls'>
            <div>
                <div className='msp-semi-transparent-background' />
                <IconButton icon={isAnimating || isPlaying ? 'stop' : 'tape'} title={isAnimating ? 'Stop' : 'Select Animation'}
                    onClick={isAnimating || isPlaying ? this.stop : this.toggleExpanded}
                    disabled={isAnimating|| isPlaying ? false : this.state.isUpdating || this.state.isPlaying || this.state.isEmpty} />
            </div>
            {(this.state.isExpanded && !this.state.isUpdating) && <div className='msp-animation-viewport-controls-select'>
                <AnimationControls onStart={this.toggleExpanded} />
            </div>}
        </div>;
    }
}

export class LociLabels extends PluginUIComponent<{}, { entries: ReadonlyArray<LociLabelEntry> }> {
    state = { entries: [] }

    componentDidMount() {
        this.subscribe(this.plugin.behaviors.labels.highlight, e => this.setState({ entries: e.entries }));
    }

    render() {
        if (this.state.entries.length === 0) {
            return null;
        }

        return <div className='msp-highlight-info'>
            {this.state.entries.map((e, i) => <div key={'' + i} dangerouslySetInnerHTML={{ __html: e }} />)}
        </div>;
    }
}

export class DefaultStructureTools extends PluginUIComponent {
    render() {
        return <>
            <div className='msp-section-header'><Icon name='tools' /> Structure Tools</div>

            <StructureSourceControls />
            <StructureSelectionControls />
            <StructureComponentControls />
            <StructureMeasurementsControls />
            <VolumeStreamingControls />
            <ObjectControls />
        </>;
    }
}
