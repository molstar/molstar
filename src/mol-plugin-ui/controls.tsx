/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { UpdateTrajectory } from '../mol-plugin-state/actions/structure';
import { LociLabel } from '../mol-plugin-state/manager/loci-label';
import { PluginStateObject } from '../mol-plugin-state/objects';
import { StateTransforms } from '../mol-plugin-state/transforms';
import { ModelFromTrajectory } from '../mol-plugin-state/transforms/model';
import { PluginCommands } from '../mol-plugin/commands';
import { StateTransformer } from '../mol-state';
import { PluginUIComponent } from './base';
import { IconButton } from './controls/common';
import { Icon, NavigateBeforeSvg, NavigateNextSvg, SkipPreviousSvg, StopSvg, PlayArrowSvg, SubscriptionsOutlinedSvg, BuildSvg } from './controls/icons';
import { AnimationControls } from './state/animation';
import { StructureComponentControls } from './structure/components';
import { StructureMeasurementsControls } from './structure/measurements';
import { StructureSelectionActionsControls } from './structure/selection';
import { StructureSourceControls } from './structure/source';
import { VolumeStreamingControls, VolumeSourceControls } from './structure/volume';
import { PluginConfig } from '../mol-plugin/config';
import { StructureSuperpositionControls } from './structure/superposition';

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

    render() {
        const isAnimating = this.plugin.behaviors.state.isAnimating.value;

        if (!this.state.show || (isAnimating && !this.state.label)) return null;

        return <div className='msp-traj-controls'>
            {!isAnimating && <IconButton svg={SkipPreviousSvg} title='First Model' onClick={this.reset} disabled={isAnimating} />}
            {!isAnimating && <IconButton svg={NavigateBeforeSvg} title='Previous Model' onClick={this.prev} disabled={isAnimating} />}
            {!isAnimating && <IconButton svg={NavigateNextSvg} title='Next Model' onClick={this.next} disabled={isAnimating} />}
            {!!this.state.label && <span>{this.state.label}</span> }
        </div>;
    }
}

export class StateSnapshotViewportControls extends PluginUIComponent<{}, { isBusy: boolean, show: boolean }> {
    state = { isBusy: false, show: true }

    componentDidMount() {
        // TODO: this needs to be diabled when the state is updating!
        this.subscribe(this.plugin.managers.snapshot.events.changed, () => this.forceUpdate());
        this.subscribe(this.plugin.behaviors.state.isBusy, isBusy => this.setState({ isBusy }));
        this.subscribe(this.plugin.behaviors.state.isAnimating, isBusy => this.setState({ isBusy }));

        window.addEventListener('keyup', this.keyUp, false);
    }

    componentWillUnmount() {
        super.componentWillUnmount();
        window.removeEventListener('keyup', this.keyUp, false);
    }

    keyUp = (e: KeyboardEvent) => {
        if (!e.ctrlKey || this.state.isBusy || e.target !== document.body) return;
        const snapshots = this.plugin.managers.snapshot;
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
        const s = this.plugin.managers.snapshot;
        const id = s.getNextId(s.state.current, -1);
        if (id) this.update(id);
    }

    next = () => {
        const s = this.plugin.managers.snapshot;
        const id = s.getNextId(s.state.current, 1);
        if (id) this.update(id);
    }

    togglePlay = () => {
        this.plugin.managers.snapshot.togglePlay();
    }

    render() {
        const snapshots = this.plugin.managers.snapshot;
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
            <IconButton svg={isPlaying ? StopSvg : PlayArrowSvg} title={isPlaying ? 'Pause' : 'Cycle States'} onClick={this.togglePlay}
                disabled={isPlaying ? false : this.state.isBusy} />
            {!isPlaying && <>
                <IconButton svg={NavigateBeforeSvg} title='Previous State' onClick={this.prev} disabled={this.state.isBusy || isPlaying} />
                <IconButton svg={NavigateNextSvg} title='Next State' onClick={this.next} disabled={this.state.isBusy || isPlaying} />
            </>}
        </div>;
    }
}

export class AnimationViewportControls extends PluginUIComponent<{}, { isEmpty: boolean, isExpanded: boolean, isBusy: boolean, isAnimating: boolean, isPlaying: boolean }> {
    state = { isEmpty: true, isExpanded: false, isBusy: false, isAnimating: false, isPlaying: false };

    componentDidMount() {
        this.subscribe(this.plugin.managers.snapshot.events.changed, () => {
            if (this.plugin.managers.snapshot.state.isPlaying) this.setState({ isPlaying: true, isExpanded: false });
            else this.setState({ isPlaying: false });
        });
        this.subscribe(this.plugin.behaviors.state.isBusy, isBusy => {
            if (isBusy) this.setState({ isBusy: true, isExpanded: false, isEmpty: this.plugin.state.data.tree.transforms.size < 2 });
            else this.setState({ isBusy: false, isEmpty: this.plugin.state.data.tree.transforms.size < 2 });
        });
        this.subscribe(this.plugin.behaviors.state.isAnimating, isAnimating => {
            if (isAnimating) this.setState({ isAnimating: true, isExpanded: false });
            else this.setState({ isAnimating: false });
        });
    }
    toggleExpanded = () => this.setState({ isExpanded: !this.state.isExpanded });
    stop = () => {
        this.plugin.managers.animation.stop();
        this.plugin.managers.snapshot.stop();
    }

    render() {
        const isPlaying = this.plugin.managers.snapshot.state.isPlaying;
        if (isPlaying || this.state.isEmpty || this.plugin.managers.animation.isEmpty || !this.plugin.config.get(PluginConfig.Viewport.ShowAnimation)) return null;

        const isAnimating = this.state.isAnimating;

        return <div className='msp-animation-viewport-controls'>
            <div>
                <div className='msp-semi-transparent-background' />
                <IconButton svg={isAnimating || isPlaying ? StopSvg : SubscriptionsOutlinedSvg} transparent title={isAnimating ? 'Stop' : 'Select Animation'}
                    onClick={isAnimating || isPlaying ? this.stop : this.toggleExpanded} toggleState={this.state.isExpanded}
                    disabled={isAnimating || isPlaying ? false : this.state.isBusy || this.state.isPlaying || this.state.isEmpty} />
            </div>
            {(this.state.isExpanded && !this.state.isBusy) && <div className='msp-animation-viewport-controls-select'>
                <AnimationControls onStart={this.toggleExpanded} />
            </div>}
        </div>;
    }
}

export class SelectionViewportControls extends PluginUIComponent {
    componentDidMount() {
        this.subscribe(this.plugin.behaviors.interaction.selectionMode, () => this.forceUpdate());
    }

    onMouseMove = (e: React.MouseEvent) => {
        // ignore mouse moves when no button is held
        if (e.buttons === 0) e.stopPropagation();
    }

    render() {
        if (!this.plugin.selectionMode) return null;
        return <div className='msp-selection-viewport-controls' onMouseMove={this.onMouseMove}>
            <StructureSelectionActionsControls />
        </div>;
    }
}

export class LociLabels extends PluginUIComponent<{}, { labels: ReadonlyArray<LociLabel> }> {
    state = { labels: [] }

    componentDidMount() {
        this.subscribe(this.plugin.behaviors.labels.highlight, e => this.setState({ labels: e.labels }));
    }

    render() {
        if (this.state.labels.length === 0) {
            return null;
        }

        return <div className='msp-highlight-info'>
            {this.state.labels.map((e, i) => <div key={'' + i} dangerouslySetInnerHTML={{ __html: e }} />)}
        </div>;
    }
}

export class CustomStructureControls extends PluginUIComponent<{ initiallyCollapsed?: boolean }> {
    componentDidMount() {
        this.subscribe(this.plugin.state.behaviors.events.changed, () => this.forceUpdate());
    }

    render() {
        const controls: JSX.Element[] = [];
        this.plugin.customStructureControls.forEach((Controls, key) => {
            controls.push(<Controls initiallyCollapsed={this.props.initiallyCollapsed} key={key} />);
        });
        return controls.length > 0 ? <>{controls}</> : null;
    }
}

export class DefaultStructureTools extends PluginUIComponent {
    render() {
        return <>
            <div className='msp-section-header'><Icon svg={BuildSvg} />Structure Tools</div>

            <StructureSourceControls />
            <StructureMeasurementsControls />
            <StructureSuperpositionControls />
            <StructureComponentControls />
            <VolumeStreamingControls />
            <VolumeSourceControls />

            <CustomStructureControls />
        </>;
    }
}
