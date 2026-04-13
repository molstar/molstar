/**
 * Kinemage right-panel controls (right-panel only).
 *
 * Shows kinemage views, animate buttons, and group/subgroup/master toggles in the right inspector.
 * Controls directly operate on the loaded kinemage runtime data and call exported helpers
 * to rebuild visuals. No State Tree JSON nodes are created for these UI items.
 */

import * as React from 'react';
import { CollapsableState, CollapsableControls } from '../../mol-plugin-ui/base';
//import { PluginContext } from '../../mol-plugin/context';
//import { Task } from '../../mol-task';
import { Camera } from '../../mol-canvas3d/camera';
import { applyViewSnapshot, createShapesForKinemage, destroyShapesForKinemage, getLoadedKinemageData } from './behavior';
//import { PluginCommands } from '../../mol-plugin/commands';

interface KinemageControlState extends CollapsableState {
  isBusy: boolean
}

export class KinemageControls extends CollapsableControls<{}, KinemageControlState> {
  protected defaultState(): KinemageControlState {
    return {
      header: 'Kinemage',
      isCollapsed: false,
      isBusy: false,
      // make the control visible even when no structure is loaded
      isHidden: false,
      brand: { accent: 'cyan', svg: undefined as any }
    };
  }

  componentDidMount() {
    // Show/hide depending on whether there are any kinemage entries loaded
    //this.subscribe(this.plugin.managers.dragAndDrop.behaviors.onDrop, () => this.updateVisibility());
    // also track plugin state changes that might affect focus/visuals
    //this.subscribe(this.plugin.state.events.updated, () => this.forceUpdate());
    //this.updateVisibility();
  }

  /*
  private updateVisibility() {
    // Keep the card visible; the content will show "No Kinemage data" when nothing is loaded.
    // If you prefer the card to hide when no kinemage is present, change this implementation.
    this.setState({ isHidden: false });
  }
  */

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

  private async toggleVisibility(kin: any, target: { type: 'group' | 'subgroup' | 'master', key: string }) {
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
    } catch (e) {
      console.error('Failed to toggle kinemage visibility', e);
    }
  }

  private async triggerAnimateForKin(kin: any, mode: 'animate' | '2animate') {
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
    } catch (e) {
      console.error('Failed to trigger animate', e);
    }
  }

  renderControls() {
    const kins = this.getKinemageList();
    if (kins.length === 0) return <div className='msp-row-text'>No Kinemage data</div>;

    const blocks: React.ReactNode[] = [];
    for (const kin of kins) {
      const title = kin.pdbfile || kin.caption || 'Kinemage';
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
            <span style={{ marginLeft: 8 }}>Animate (change vis)</span>
          </div>
        );
      }
      if (kin.groupsAnimate2 && kin.groupsAnimate2.length > 0) {
        kinBlock.push(
          <div key={'anim2-' + title} className='msp-row'>
            <button className='msp-button' onClick={() => this.triggerAnimateForKin(kin, '2animate')}>Animate2</button>
            <span style={{ marginLeft: 8 }}>Animate2 (change vis)</span>
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
