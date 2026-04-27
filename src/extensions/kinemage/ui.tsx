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
 * to rebuild visuals. No State Tree JSON nodes are created for these UI items.
 */

import * as React from 'react';
import { CollapsableState, CollapsableControls } from '../../mol-plugin-ui/base';
import { Camera } from '../../mol-canvas3d/camera';
import { applyViewSnapshot, rebuildShapesForKinemage } from './behavior';
import { Kinemage } from './reader/schema';

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
      // default hidden until a kinemage is present
      isHidden: true,
      brand: { accent: 'cyan', svg: undefined as any }
    };
  }

  componentDidMount() {
    // Listen for shape/state changes: when state tree cells are created or removed the visuals changed.
    this.subscribe(this.plugin.state.data.events.cell.created, (e: any) => this.onCellCreated(e));
    this.subscribe(this.plugin.state.data.events.cell.removed, () => this.onCellRemoved());
    // also track cell state updates that may change labels / visibility
    this.subscribe(this.plugin.state.data.events.cell.stateUpdated, () => this.forceUpdate());

    // ensure initial visibility reflects current state
    this.updateVisibility();
  }

  private onCellCreated(e: any) {
    this.updateVisibility();
  }

  private onCellRemoved() {
    this.updateVisibility();
  }

  private updateVisibility() {
    const kinemages = this.getKinemageList();
    this.setState({ isHidden: kinemages.length === 0 });
  }

  private getKinemageList(): Array<{ kinData: Kinemage, ref: string }> {
    const result: Array<{ kinData: Kinemage, ref: string }> = [];

    try {
      const cells = (this.plugin.state.data as any).cells as Map<string, any>;
      for (const [ref, entry] of cells) {
        const obj = (entry as any).obj;
        // Look for Format.Json nodes that contain kinData
        if (obj && obj.data && (obj.data as any).kinData) {
          result.push({ kinData: (obj.data as any).kinData, ref });
        }
      }
    } catch (e) {
      console.warn('Failed to enumerate kinemage nodes', e);
    }

    return result;
  }

  private async applyView(kinData: Kinemage, viewKey: string) {
    const snap = (kinData as any).viewSnapshots?.[viewKey];
    if (snap) {
      await applyViewSnapshot(this.plugin, snap as Partial<Camera.Snapshot>);
    }
  }

  private async toggleVisibility(kinData: Kinemage, kinRef: string, target: { type: 'group' | 'subgroup' | 'master', key: string }) {
    try {
      if (target.type === 'group') {
        const g = kinData.groupDict[target.key];
        if (g) g.off = !g.off;
      } else if (target.type === 'subgroup') {
        const s = kinData.subgroupDict[target.key];
        if (s) s.off = !s.off;
      } else {
        const m = kinData.masterDict[target.key];
        if (m) m.visible = !m.visible;
      }

      // Rebuild shapes for this kinemage using the state ref
      await rebuildShapesForKinemage(this.plugin, { ref: kinRef } as any);
      this.updateVisibility();
    } catch (e) {
      console.error('Failed to toggle kinemage visibility', e);
    }
  }

  private async triggerAnimateForKin(kinData: Kinemage, kinRef: string, mode: 'animate' | '2animate') {
    try {
      if (mode === 'animate') {
        kinData.activeAnimateGroup = (kinData.activeAnimateGroup + 1) % Math.max(1, kinData.groupsAnimate.length);

        // Make only the active animate group visible, hide the others (if any)
        for (let i = 0; i < kinData.groupsAnimate.length; i++) {
          const groupName = kinData.groupsAnimate[i];
          const groupInfo = kinData.groupDict[groupName];
          if (groupInfo) {
            groupInfo.off = (i !== kinData.activeAnimateGroup);
          }
        }
      } else {
        kinData.activeAnimateGroup2 = (kinData.activeAnimateGroup2 + 1) % Math.max(1, kinData.groupsAnimate2.length);

        // Make only the active animate group visible, hide the others (if any)
        for (let i = 0; i < kinData.groupsAnimate2.length; i++) {
          const groupName = kinData.groupsAnimate2[i];
          const groupInfo = kinData.groupDict[groupName];
          if (groupInfo) {
            groupInfo.off = (i !== kinData.activeAnimateGroup2);
          }
        }
      }

      // Rebuild shapes for this kinemage using the state ref
      await rebuildShapesForKinemage(this.plugin, { ref: kinRef } as any);
      this.updateVisibility();
    } catch (e) {
      console.error('Failed to trigger animate', e);
    }
  }

  renderControls() {
    const kins = this.getKinemageList();
    if (kins.length === 0) return <div className='msp-row-text'>No Kinemage data</div>;

    const blocks: React.ReactNode[] = [];
    for (const { kinData, ref } of kins) {
      const title = kinData.pdbfile || nameFromString(kinData.caption) || 'Kinemage';
      const kinBlock: React.ReactNode[] = [];
      kinBlock.push(<div key={'title-' + title} className='msp-row-text'><b>{title}</b></div>);

      // views
      for (const [viewKey, viewObj] of Object.entries(kinData.viewDict || {})) {
        const label = viewObj.name || `View ${viewKey}`;
        kinBlock.push(
          <div key={'view-' + title + '-' + viewKey} className='msp-row'>
            <button className='msp-button' onClick={() => this.applyView(kinData, viewKey)}>Apply View</button>
            <span style={{ marginLeft: 8 }}>{label}</span>
          </div>
        );
      }

      // animate
      if (kinData.groupsAnimate && kinData.groupsAnimate.length > 0) {
        kinBlock.push(
          <div key={'anim-' + title} className='msp-row'>
            <button className='msp-button' onClick={() => this.triggerAnimateForKin(kinData, ref, 'animate')}>Animate</button>
            <span style={{ marginLeft: 8 }}>Animate</span>
          </div>
        );
      }
      if (kinData.groupsAnimate2 && kinData.groupsAnimate2.length > 0) {
        kinBlock.push(
          <div key={'anim2-' + title} className='msp-row'>
            <button className='msp-button' onClick={() => this.triggerAnimateForKin(kinData, ref, '2animate')}>Animate2</button>
            <span style={{ marginLeft: 8 }}>Animate2</span>
          </div>
        );
      }

      // groups
      for (const [groupKey, groupInfo] of Object.entries(kinData.groupDict || {})) {
        if (!(groupInfo as any).nobutton) {
          const visible = !(groupInfo as any).off;
          // If this group is in animate or animate2, then add '*' before its groupKey name to indicate that it's an animation group
          const isAnimate = kinData.groupsAnimate?.includes(groupKey) || kinData.groupsAnimate2?.includes(groupKey);
          const label = isAnimate ? `* ${groupKey}` : groupKey;
          kinBlock.push(
            <div key={'group-' + title + '-' + groupKey} className='msp-row'>
              <label style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                <input type='checkbox' checked={visible} onChange={() => this.toggleVisibility(kinData, ref, { type: 'group', key: groupKey })} />
                <span>{label}</span>
              </label>
            </div>
          );
        }
        // If this group is not dominant, find any subgroups of this group and show them here (indented) unless they have nobutton set
        if (!(groupInfo as any).dominant) {
          for (const [subgroupKey, subgroupInfo] of Object.entries(kinData.subgroupDict || {})) {
            if (subgroupKey.startsWith(groupKey + ':')) {
              if ((subgroupInfo as any).nobutton) continue;
              const visible = !(subgroupInfo as any).off;
              kinBlock.push(
                <div key={'subgroup-' + title + '-' + subgroupKey} className='msp-row' style={{ paddingLeft: 16 }}>
                  <label style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
                    <input type='checkbox' checked={visible} onChange={() => this.toggleVisibility(kinData, ref, { type: 'subgroup', key: subgroupKey })} />
                    <span>{subgroupKey.split(':')[1]}</span>
                  </label>
                </div>
              );
            }
          }
        }
      }

      // subgroups that don't belong to a group (standalone)
      for (const [subgroupKey, subgroupInfo] of Object.entries(kinData.subgroupDict || {})) {
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
              <input type='checkbox' checked={visible} onChange={() => this.toggleVisibility(kinData, ref, { type: 'subgroup', key: subgroupKey })} />
              <span>{subgroupKey}</span>
            </label>
          </div>
        );
      }

      // masters
      for (const [masterKey, masterInfo] of Object.entries(kinData.masterDict || {})) {
        const visible = !!(masterInfo && (masterInfo as any).visible);
        kinBlock.push(
          <div key={'master-' + title + '-' + masterKey} className='msp-row'>
            <label style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
              <input type='checkbox' checked={visible} onChange={() => this.toggleVisibility(kinData, ref, { type: 'master', key: masterKey })} />
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