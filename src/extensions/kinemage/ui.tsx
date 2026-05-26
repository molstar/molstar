/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Russ Taylor <russ@reliasolve.com>
 */

/**
 * Kinemage right-panel controls (right-panel only).
 *
 * Shows kinemage views, animate buttons, and group/subgroup/master toggles in the right inspector.
 * Controls update visibility controller parameters which trigger rebuilds via the state tree.
 */

import * as React from 'react';
import { CollapsableState, CollapsableControls } from '../../mol-plugin-ui/base';
import { Camera } from '../../mol-canvas3d/camera';
import { applyViewSnapshot } from './behavior';
import { Kinemage } from './reader/schema';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { KinemageShapePointsProvider, KinemageShapeLinesProvider, KinemageShapeMeshProvider, KinemageShapeSpheresProvider } from './behavior';

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

  private getKinemageList(): Array<{ kinData: Kinemage, ref: string, visControllerRef: string }> {
    const result: Array<{ kinData: Kinemage, ref: string, visControllerRef: string }> = [];

    try {
      const cells = (this.plugin.state.data as any).cells as Map<string, any>;
      for (const [ref, entry] of cells) {
        const obj = (entry as any).obj;
        // Look for Format.Json nodes that contain kinData and visibilityState (visibility controller)
        if (obj && obj.data && (obj.data as any).kinData && (obj.data as any).visibilityState) {
          result.push({
            kinData: (obj.data as any).kinData,
            ref,
            visControllerRef: ref
          });
        }
      }
    } catch (e) {
      console.warn('Failed to enumerate kinemage nodes', e);
    }

    return result;
  }

  private getAllDescendants(nodeRef: string): string[] {
    const result: string[] = [];
    const tree = this.plugin.state.data.tree;
    const queue = [nodeRef];

    while (queue.length > 0) {
      const current = queue.shift()!;
      const children = tree.children.get(current);
      if (children) {
        for (const childRef of children.values()) {
          result.push(childRef);
          queue.push(childRef);
        }
      }
    }

    return result;
  }

  private async applyView(kinData: Kinemage, viewKey: string) {
    const snap = (kinData as any).viewSnapshots?.[viewKey];
    if (snap) {
      await applyViewSnapshot(this.plugin, snap as Partial<Camera.Snapshot>);
    }
  }

  private async rebuildShapes(visControllerRef: string, kinData: Kinemage) {
    const update = this.plugin.state.data.build();

    // Delete all descendants (shape providers and representations)
    const descendants = this.getAllDescendants(visControllerRef);
    for (const nodeRef of descendants) {
      update.delete(nodeRef);
    }

    await update.commit();

    // Recreate shapes
    const rebuildUpdate = this.plugin.state.data.build();

    // Generate all shape types that have data, each as child of the visibility controller
    if (kinData.dotLists.length > 0) {
      rebuildUpdate
        .to(visControllerRef)
        .apply(KinemageShapePointsProvider, {}, { state: { isGhost: true } })
        .apply(StateTransforms.Representation.ShapeRepresentation3D);
    }
    if (kinData.vectorLists.length > 0) {
      rebuildUpdate
        .to(visControllerRef)
        .apply(KinemageShapeLinesProvider, {}, { state: { isGhost: true } })
        .apply(StateTransforms.Representation.ShapeRepresentation3D);
    }
    if (kinData.ribbonLists.length > 0) {
      rebuildUpdate
        .to(visControllerRef)
        .apply(KinemageShapeMeshProvider, {}, { state: { isGhost: true } })
        .apply(StateTransforms.Representation.ShapeRepresentation3D, { doubleSided: true });
    }
    if (kinData.ballLists.length > 0) {
      rebuildUpdate
        .to(visControllerRef)
        .apply(KinemageShapeSpheresProvider, {}, { state: { isGhost: true } })
        .apply(StateTransforms.Representation.ShapeRepresentation3D);
    }

    await rebuildUpdate.commit();
  }

  private async toggleVisibility(visControllerRef: string, kinData: Kinemage, target: { type: 'group' | 'subgroup' | 'master', key: string }) {
    try {
      const cell = this.plugin.state.data.cells.get(visControllerRef);
      if (!cell || !cell.transform || !cell.transform.params) return;

      const currentParams = cell.transform.params;
      const newGroupVisibility = { ...currentParams.groupVisibility };
      const newSubgroupVisibility = { ...currentParams.subgroupVisibility };
      const newMasterVisibility = { ...currentParams.masterVisibility };

      if (target.type === 'group') {
        newGroupVisibility[target.key] = !newGroupVisibility[target.key];
      } else if (target.type === 'subgroup') {
        newSubgroupVisibility[target.key] = !newSubgroupVisibility[target.key];
      } else {
        newMasterVisibility[target.key] = !newMasterVisibility[target.key];
      }

      const update = this.plugin.state.data.build();

      // Update the visibility controller
      update.to(visControllerRef).update({
        groupVisibility: newGroupVisibility,
        subgroupVisibility: newSubgroupVisibility,
        masterVisibility: newMasterVisibility
      });

      await update.commit();

      // Rebuild all shapes to reflect new visibility
      await this.rebuildShapes(visControllerRef, kinData);
    } catch (e) {
      console.error('Failed to toggle kinemage visibility', e);
    }
  }

  private async triggerAnimateForKin(visControllerRef: string, kinData: Kinemage, mode: 'animate' | '2animate') {
    try {
      const cell = this.plugin.state.data.cells.get(visControllerRef);
      if (!cell || !cell.transform || !cell.transform.params) return;

      const currentParams = cell.transform.params;
      const animateGroups = mode === 'animate' ? kinData.groupsAnimate : kinData.groupsAnimate2;
      const currentActive = mode === 'animate' ? currentParams.activeAnimateGroup : currentParams.activeAnimateGroup2;
      const nextActive = (currentActive + 1) % Math.max(1, animateGroups.length);

      // Update group visibility to show only the active animate group
      const newGroupVisibility = { ...currentParams.groupVisibility };

      // Hide all animate groups, then show only the active one
      for (let i = 0; i < animateGroups.length; i++) {
        newGroupVisibility[animateGroups[i]] = (i === nextActive);
      }

      const update = this.plugin.state.data.build();

      // Update the visibility controller
      if (mode === 'animate') {
        update.to(visControllerRef).update({
          groupVisibility: newGroupVisibility,
          activeAnimateGroup: nextActive
        });
      } else {
        update.to(visControllerRef).update({
          groupVisibility: newGroupVisibility,
          activeAnimateGroup2: nextActive
        });
      }

      await update.commit();

      // Rebuild all shapes to reflect new visibility
      await this.rebuildShapes(visControllerRef, kinData);
    } catch (e) {
      console.error('Failed to trigger animate', e);
    }
  }

  private isVisible(visControllerRef: string, target: { type: 'group' | 'subgroup' | 'master', key: string }): boolean {
    try {
      const cell = this.plugin.state.data.cells.get(visControllerRef);
      if (!cell || !cell.transform || !cell.transform.params) return true;

      const params = cell.transform.params;
      if (target.type === 'group') {
        return params.groupVisibility[target.key] !== false;
      } else if (target.type === 'subgroup') {
        return params.subgroupVisibility[target.key] !== false;
      } else {
        return params.masterVisibility[target.key] !== false;
      }
    } catch (e) {
      return true;
    }
  }

  renderControls() {
    const kins = this.getKinemageList();
    if (kins.length === 0) return <div style={{ padding: '6px' }}>No Kinemage data</div>;

    const blocks: React.ReactNode[] = [];
    for (const { kinData, visControllerRef } of kins) {
      const title = kinData.pdbfile || nameFromString(kinData.caption) || 'Kinemage';
      const kinBlock: React.ReactNode[] = [];

      // Title
      kinBlock.push(
        <div key={'title-' + title} style={{ padding: '6px', fontWeight: 'bold', borderBottom: '1px solid rgba(255,255,255,0.1)' }}>
          {title}
        </div>
      );

      // views
      const viewEntries = Object.entries(kinData.viewDict || {});
      if (viewEntries.length > 0) {
        for (const [viewKey, viewObj] of viewEntries) {
          const label = `View ${viewObj.name || `View ${viewKey}`}`;
          kinBlock.push(
            <div key={'view-' + title + '-' + viewKey} style={{ padding: '2px 6px' }}>
              <button
                className='msp-btn msp-btn-block'
                onClick={() => this.applyView(kinData, viewKey)}
                title={`Apply view: ${label}`}
              >
                {label}
              </button>
            </div>
          );
        }
      }

      // animate
      if (kinData.groupsAnimate && kinData.groupsAnimate.length > 0) {
        kinBlock.push(
          <div key={'anim-' + title} style={{ padding: '2px 6px' }}>
            <button
              className='msp-btn msp-btn-block'
              onClick={() => this.triggerAnimateForKin(visControllerRef, kinData, 'animate')}
              title='Cycle through animation frames'
            >
              Animate
            </button>
          </div>
        );
      }
      if (kinData.groupsAnimate2 && kinData.groupsAnimate2.length > 0) {
        kinBlock.push(
          <div key={'anim2-' + title} style={{ padding: '2px 6px' }}>
            <button
              className='msp-btn msp-btn-block'
              onClick={() => this.triggerAnimateForKin(visControllerRef, kinData, '2animate')}
              title='Cycle through second animation frames'
            >
              Animate2
            </button>
          </div>
        );
      }

      // groups
      for (const [groupKey, groupInfo] of Object.entries(kinData.groupDict || {})) {
        if (!(groupInfo as any).nobutton) {
          const visible = this.isVisible(visControllerRef, { type: 'group', key: groupKey });
          // If this group is in animate or animate2, then add '*' before its groupKey name to indicate that it's an animation group
          const isAnimate = (kinData.groupsAnimate?.includes(groupKey) ?? false) || (kinData.groupsAnimate2?.includes(groupKey) ?? false);
          const label = isAnimate ? `* ${groupKey}` : groupKey;
          kinBlock.push(
            <div key={'group-' + title + '-' + groupKey} style={{ padding: '2px 6px' }}>
              <label style={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}>
                <input
                  type='checkbox'
                  checked={visible}
                  onChange={() => this.toggleVisibility(visControllerRef, kinData, { type: 'group', key: groupKey })}
                  style={{ marginRight: '6px' }}
                />
                <span title={label}>{label}</span>
              </label>
            </div>
          );
        }
        // If this group is not dominant, find any subgroups of this group and show them here (indented) unless they have nobutton set
        if (!(groupInfo as any).dominant) {
          for (const [subgroupKey, subgroupInfo] of Object.entries(kinData.subgroupDict || {})) {
            if (subgroupKey.startsWith(groupKey + ':')) {
              if ((subgroupInfo as any).nobutton) continue;
              const visible = this.isVisible(visControllerRef, { type: 'subgroup', key: subgroupKey });
              const subgroupLabel = subgroupKey.split(':')[1];
              kinBlock.push(
                <div key={'subgroup-' + title + '-' + subgroupKey} style={{ padding: '2px 6px', paddingLeft: '24px' }}>
                  <label style={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}>
                    <input
                      type='checkbox'
                      checked={visible}
                      onChange={() => this.toggleVisibility(visControllerRef, kinData, { type: 'subgroup', key: subgroupKey })}
                      style={{ marginRight: '6px' }}
                    />
                    <span title={subgroupLabel}>{subgroupLabel}</span>
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
        const visible = this.isVisible(visControllerRef, { type: 'subgroup', key: subgroupKey });
        kinBlock.push(
          <div key={'subgroup-' + title + '-' + subgroupKey} style={{ padding: '2px 6px' }}>
            <label style={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}>
              <input
                type='checkbox'
                checked={visible}
                onChange={() => this.toggleVisibility(visControllerRef, kinData, { type: 'subgroup', key: subgroupKey })}
                style={{ marginRight: '6px' }}
              />
              <span title={subgroupKey}>{subgroupKey}</span>
            </label>
          </div>
        );
      }

      // masters
      for (const [masterKey] of Object.entries(kinData.masterDict || {})) {
        const visible = this.isVisible(visControllerRef, { type: 'master', key: masterKey });
        kinBlock.push(
          <div key={'master-' + title + '-' + masterKey} style={{ padding: '2px 6px' }}>
            <label style={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}>
              <input
                type='checkbox'
                checked={visible}
                onChange={() => this.toggleVisibility(visControllerRef, kinData, { type: 'master', key: masterKey })}
                style={{ marginRight: '6px' }}
              />
              <span title={masterKey}>{masterKey}</span>
            </label>
          </div>
        );
      }

      blocks.push(<div key={'kin-block-' + title} className='msp-control-group-wrapper'>{kinBlock}</div>);
    }

    return <>{blocks}</>;
  }
}