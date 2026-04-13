/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Russ Taylor <russ@reliasolve.com>
 */

/**
 * Kinemage right-panel controls (right-panel only).
 *
 * Shows kinemage views, animate buttons, and group/subgroup/master toggles in the right inspector.
 * Controls directly operate on the loaded kinemage runtime data and call exported helpers
 * to rebuild visuals.
 */

import * as React from 'react';
import { CollapsableState, CollapsableControls } from '../../mol-plugin-ui/base';
import { Camera } from '../../mol-canvas3d/camera';
import { applyViewSnapshot, createShapesForKinemage, destroyShapesForKinemage, getLoadedKinemageData } from './behavior';

interface KinemageControlState extends CollapsableState {
    isBusy: boolean
}

function nameFromString(s: string | undefined) {
  // If this is undefined, return undefined.
  if (!s) return undefined;

  // Return up to the first 30 characters of the string.
  return s.length > 30 ? s.substring(0, 30) + '...' : s;
}

export class KinemageControls extends CollapsableControls<{}, KinemageControlState> {
    protected defaultState(): KinemageControlState {
        return {
            header: 'Kinemage',
            isCollapsed: false,
            isBusy: false,
            // make the control invisible until a Kinemage is loaded
            isHidden: true,
            brand: { accent: 'cyan', svg: undefined as any }
        };
    }

    componentDidMount() {
        // Listen for shape/state changes: when state tree cells are created or removed the visuals changed.
        // This is sufficient to update the UI after applyKinemageInfoToState creates/destroys transforms.
        //this.subscribe(this.plugin.state.data.events.cell.removed, () => this.updateVisibility());
        // also track cell state updates that may change labels / visibility
        this.subscribe(this.plugin.state.data.events.cell.stateUpdated, () => this.forceUpdate());
        this.updateVisibility();
    }

    private updateVisibility() {
        // Keep the card visible; the content will show "No Kinemage data" when nothing is loaded.
        // If you prefer the card to hide when no kinemage is present, change this implementation.
        this.setState({ isHidden: false });
    }

    private getKinemageList() {
        const data = getLoadedKinemageData();
        return (data && data.kinemages) ? data.kinemages : [];
    }

    private async applyView(kin: any, viewKey: string) {
        const snap = (kin as any).viewSnapshots?.[viewKey];
        if (snap) {
            await applyViewSnapshot(this.plugin, snap as Partial<Camera.Snapshot>);
        }
    }

    private async toggleVisibility(kin: any, target: { type: 'group'|'subgroup'|'master', key: string }) {
        try {
            if (target.type === 'group') {
                const g = kin.groupDict[target.key];
                if (g) g.off = !g.off;
            } else if (target.type === 'subgroup') {
                const s = kin.subgroupDict[target.key];
                if (s) s.off = !s.off;
            } else {
                const m = kin.masterDict[target.key];
                if (m) m.visible = !m.visible;
            }

            // rebuild shapes for this kinemage
            const update = this.plugin.state.data.build();
            await destroyShapesForKinemage(this.plugin, kin);
            await createShapesForKinemage(this.plugin, update, kin);
            await update.commit();
            this.updateVisibility();
        } catch (e) {
            console.error('Failed to toggle kinemage visibility', e);
        }
    }

    private async triggerAnimateForKin(kin: any, mode: 'animate'|'2animate') {
        try {
            if (mode === 'animate') {
                kin.activeAnimateGroup = (kin.activeAnimateGroup + 1) % Math.max(1, kin.groupsAnimate.length);
            } else {
                kin.activeAnimateGroup2 = (kin.activeAnimateGroup2 + 1) % Math.max(1, kin.groupsAnimate2.length);
            }

            const update = this.plugin.state.data.build();
            await destroyShapesForKinemage(this.plugin, kin);
            await createShapesForKinemage(this.plugin, update, kin);
            await update.commit();
            this.updateVisibility();
        } catch (e) {
            console.error('Failed to trigger animate', e);
        }
    }

    renderControls() {
        const kins = this.getKinemageList();
        if (kins.length === 0) return <div className='msp-row-text'>No Kinemage data</div>;

        const blocks: React.ReactNode[] = [];
        for (const kin of kins) {
            const title = kin.pdbfile || nameFromString(kin.caption) || 'Kinemage';
            const kinBlock: React.ReactNode[] = [];
            kinBlock.push(<div key={'title-' + title} className='msp-row-text'><b>{title}</b></div>);

            // views
            for (const [viewKey, viewObj] of Object.entries(kin.viewDict || {})) {
                const label = viewObj.name || `View ${viewKey}`;
                kinBlock.push(
                    <div key={'view-' + title + '-' + viewKey} className='msp-row'>
                        <button className='msp-button' onClick={() => this.applyView(kin, viewKey)}>Apply View</button>
                        <span style={{ marginLeft: 8 }}>{label}</span>
                    </div>
                );
            }

            // animate
            if (kin.groupsAnimate && kin.groupsAnimate.length > 0) {
                kinBlock.push(
                    <div key={'anim-' + title} className='msp-row'>
                        <button className='msp-button' onClick={() => this.triggerAnimateForKin(kin, 'animate')}>Animate</button>
                        <span style={{ marginLeft: 8 }}>Animate</span>
                    </div>
                );
            }
            if (kin.groupsAnimate2 && kin.groupsAnimate2.length > 0) {
                kinBlock.push(
                    <div key={'anim2-' + title} className='msp-row'>
                        <button className='msp-button' onClick={() => this.triggerAnimateForKin(kin, '2animate')}>Animate2</button>
                        <span style={{ marginLeft: 8 }}>Animate2</span>
                    </div>
                );
            }

            // groups
            for (const [groupKey, groupInfo] of Object.entries(kin.groupDict || {})) {
                if ((groupInfo as any).nobutton) continue;
                const visible = !(groupInfo as any).off;
                kinBlock.push(
                    <div key={'group-' + title + '-' + groupKey} className='msp-row'>
                        <label style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                            <input type='checkbox' checked={visible} onChange={() => this.toggleVisibility(kin, { type: 'group', key: groupKey })} />
                            <span>{groupKey}</span>
                        </label>
                    </div>
                );
            }

            // subgroups that don't belong to a group (standalone)
            for (const [subgroupKey, subgroupInfo] of Object.entries(kin.subgroupDict || {})) {
                // if parent group present, those groups' subgroups are already shown when iterating groups
                if (subgroupKey.indexOf(':') !== -1) {
                    // subgroups with parent group; skip here (shown under parent group)
                    continue;
                }
                if ((subgroupInfo as any).nobutton) continue;
                const visible = !(subgroupInfo as any).off;
                kinBlock.push(
                    <div key={'subgroup-' + title + '-' + subgroupKey} className='msp-row'>
                        <label style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                            <input type='checkbox' checked={visible} onChange={() => this.toggleVisibility(kin, { type: 'subgroup', key: subgroupKey })} />
                            <span>{subgroupKey}</span>
                        </label>
                    </div>
                );
            }

            // masters
            for (const [masterKey, masterInfo] of Object.entries(kin.masterDict || {})) {
                const visible = !!(masterInfo && (masterInfo as any).visible);
                kinBlock.push(
                    <div key={'master-' + title + '-' + masterKey} className='msp-row'>
                        <label style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                            <input type='checkbox' checked={visible} onChange={() => this.toggleVisibility(kin, { type: 'master', key: masterKey })} />
                            <span>{masterKey}</span>
                        </label>
                    </div>
                );
            }

            blocks.push(<div key={'kin-block-' + title} style={{ marginBottom: 8 }}>{kinBlock}</div>);
        }

        return <div>{blocks}</div>;
    }
}
