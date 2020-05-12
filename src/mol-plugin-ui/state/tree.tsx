/**
 * Copyright (c) 2018 - 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react';
import { debounceTime, filter } from 'rxjs/operators';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { PluginCommands } from '../../mol-plugin/commands';
import { State, StateAction, StateObject, StateObjectCell, StateTransform } from '../../mol-state';
import { StateTreeSpine } from '../../mol-state/tree/spine';
import { PluginUIComponent, _Props, _State } from '../base';
import { ActionMenu } from '../controls/action-menu';
import { Button, ControlGroup, IconButton } from '../controls/common';
import { Icon, HomeOutlinedSvg, ArrowRightSvg, ArrowDropDownSvg, DeleteOutlinedSvg, VisibilityOffOutlinedSvg, VisibilityOutlinedSvg, CloseSvg } from '../controls/icons';
import { ApplyActionControl } from './apply-action';
import { UpdateTransformControl } from './update-transform';

export class StateTree extends PluginUIComponent<{ state: State }, { showActions: boolean }> {
    state = { showActions: true };

    componentDidMount() {
        this.subscribe(this.plugin.state.events.cell.created, e => {
            if (e.cell.transform.parent === StateTransform.RootRef) this.forceUpdate();
        });

        this.subscribe(this.plugin.state.events.cell.removed, e => {
            if (e.parent === StateTransform.RootRef) this.forceUpdate();
        });
    }

    static getDerivedStateFromProps(props: { state: State }, state: { showActions: boolean }) {
        const n = props.state.tree.root.ref;
        const children = props.state.tree.children.get(n);
        const showActions = children.size === 0;
        if (state.showActions === showActions) return null;
        return { showActions };
    }

    render() {
        const ref = this.props.state.tree.root.ref;
        if (this.state.showActions) {
            return <div style={{ margin: '10px', cursor: 'default' }}>
                <p>Nothing to see here yet.</p>
                <p>Structures and Volumes can be loaded from the <Icon svg={HomeOutlinedSvg} /> tab.</p>
            </div>;
        }
        return <StateTreeNode cell={this.props.state.cells.get(ref)!} depth={0} />;
    }
}

class StateTreeNode extends PluginUIComponent<{ cell: StateObjectCell, depth: number }, { isCollapsed: boolean, isNull: boolean, showLabel: boolean }> {
    is(e: State.ObjectEvent) {
        return e.ref === this.ref && e.state === this.props.cell.parent;
    }

    get ref() {
        return this.props.cell.transform.ref;
    }

    componentDidMount() {
        this.subscribe(this.plugin.state.events.cell.stateUpdated, e => {
            if (this.props.cell === e.cell && this.is(e) && e.state.cells.has(this.ref)) {
                if (this.state.isCollapsed !== !!e.cell.state.isCollapsed
                    || this.state.isNull !== StateTreeNode.isNull(e.cell)
                    || this.state.showLabel !== StateTreeNode.showLabel(e.cell)) {
                    this.forceUpdate();
                }
            }
        });

        this.subscribe(this.plugin.state.events.cell.created, e => {
            if (this.props.cell.parent === e.state && this.ref === e.cell.transform.parent) {
                this.forceUpdate();
            }
        });

        this.subscribe(this.plugin.state.events.cell.removed, e => {
            if (this.props.cell.parent === e.state && this.ref === e.parent) {
                this.forceUpdate();
            }
        });
    }

    state = {
        isCollapsed: !!this.props.cell.state.isCollapsed,
        isNull: StateTreeNode.isNull(this.props.cell),
        showLabel: StateTreeNode.showLabel(this.props.cell)
    }

    static getDerivedStateFromProps(props: _Props<StateTreeNode>, state: _State<StateTreeNode>): _State<StateTreeNode> | null {
        const isNull = StateTreeNode.isNull(props.cell);
        const showLabel = StateTreeNode.showLabel(props.cell);
        if (!!props.cell.state.isCollapsed === state.isCollapsed && state.isNull === isNull && state.showLabel === showLabel) return null;
        return { isCollapsed: !!props.cell.state.isCollapsed, isNull, showLabel };
    }

    private static hasDecorator(cell: StateObjectCell) {
        const children = cell.parent!.tree.children.get(cell.transform.ref);
        if (children.size !== 1) return false;
        return !!cell.parent?.tree.transforms.get(children.first()).transformer.definition.isDecorator;
    }

    private static isNull(cell?: StateObjectCell) {
        return !cell || !cell.parent || cell.obj === StateObject.Null || !cell.parent.tree.transforms.has(cell.transform.ref);
    }

    private static showLabel(cell: StateObjectCell) {
        return (cell.transform.ref !== StateTransform.RootRef) && (cell.status !== 'ok' || (!cell.state.isGhost && !StateTreeNode.hasDecorator(cell)));
    }

    render() {
        if (this.state.isNull) {
            return null;
        }

        const cell = this.props.cell!;
        const children = cell.parent!.tree.children.get(this.ref);

        if (!this.state.showLabel) {
            if (children.size === 0) return null;
            return <>
                {children.map(c => <StateTreeNode cell={cell.parent!.cells.get(c!)!} key={c} depth={this.props.depth} />)}
            </>;
        }

        const newDepth = this.props.depth + 1;
        return <>
            <StateTreeNodeLabel cell={cell} depth={this.props.depth} />
            {children.size === 0
                ? void 0
                : <div style={{ display: this.state.isCollapsed ? 'none' : 'block' }}>
                    {children.map(c => <StateTreeNode cell={cell.parent!.cells.get(c!)!} key={c} depth={newDepth} />)}
                </div>
            }
        </>;
    }
}

interface StateTreeNodeLabelState {
    isCurrent: boolean,
    isCollapsed: boolean,
    action?: 'options' | 'apply',
    currentAction?: StateAction
}

class StateTreeNodeLabel extends PluginUIComponent<{ cell: StateObjectCell, depth: number }, StateTreeNodeLabelState> {

    is(e: State.ObjectEvent) {
        return e.ref === this.ref && e.state === this.props.cell.parent;
    }

    get ref() {
        return this.props.cell.transform.ref;
    }

    componentDidMount() {
        this.subscribe(this.plugin.state.events.cell.stateUpdated.pipe(filter(e => this.is(e)), debounceTime(33)), e => {
            this.forceUpdate();
        });

        this.subscribe(this.props.cell.parent!.behaviors.currentObject, e => {
            if (!this.is(e)) {
                if (this.state.isCurrent && e.state.transforms.has(this.ref)) {
                    this._setCurrent(this.props.cell.parent!.current === this.ref, this.state.isCollapsed);
                }
                return;
            }

            if (e.state.transforms.has(this.ref)) {
                this._setCurrent(this.props.cell.parent!.current === this.ref, !!this.props.cell.state.isCollapsed);
            }
        });
    }

    private _setCurrent(isCurrent: boolean, isCollapsed: boolean) {
        if (isCurrent) {
            this.setState({ isCurrent, action: 'options', currentAction: void 0, isCollapsed });
        } else {
            this.setState({ isCurrent, action: void 0, currentAction: void 0, isCollapsed });
        }
    }

    state: StateTreeNodeLabelState = {
        isCurrent: this.props.cell.parent!.current === this.ref,
        isCollapsed: !!this.props.cell.state.isCollapsed,
        action: void 0,
        currentAction: void 0 as StateAction | undefined
    }

    static getDerivedStateFromProps(props: _Props<StateTreeNodeLabel>, state: _State<StateTreeNodeLabel>): _State<StateTreeNodeLabel> | null {
        const isCurrent = props.cell.parent!.current === props.cell.transform.ref;
        const isCollapsed = !!props.cell.state.isCollapsed;

        if (state.isCollapsed === isCollapsed && state.isCurrent === isCurrent) return null;
        return { isCurrent, isCollapsed, action: void 0, currentAction: void 0 };
    }

    setCurrent = (e?: React.MouseEvent<HTMLElement>) => {
        e?.preventDefault();
        e?.currentTarget.blur();
        PluginCommands.State.SetCurrentObject(this.plugin, { state: this.props.cell.parent!, ref: this.ref });
    }

    setCurrentRoot = (e?: React.MouseEvent<HTMLElement>) => {
        e?.preventDefault();
        e?.currentTarget.blur();
        PluginCommands.State.SetCurrentObject(this.plugin, { state: this.props.cell.parent!, ref: StateTransform.RootRef });
    }

    remove = (e?: React.MouseEvent<HTMLElement>) => {
        e?.preventDefault();
        PluginCommands.State.RemoveObject(this.plugin, { state: this.props.cell.parent!, ref: this.ref, removeParentGhosts: true });
    }

    toggleVisible = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.ToggleVisibility(this.plugin, { state: this.props.cell.parent!, ref: this.ref });
        e.currentTarget.blur();
    }

    toggleExpanded = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.ToggleExpanded(this.plugin, { state: this.props.cell.parent!, ref: this.ref });
        e.currentTarget.blur();
    }

    highlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.Interactivity.Object.Highlight(this.plugin, { state: this.props.cell.parent!, ref: this.ref });
        e.currentTarget.blur();
    }

    clearHighlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.Interactivity.ClearHighlights(this.plugin);
        e.currentTarget.blur();
    }

    hideApply = () => {
        this.setCurrentRoot();
    }

    get actions() {
        const cell = this.props.cell;
        const actions = [...cell.parent!.actions.fromCell(cell, this.plugin)];
        if (actions.length === 0) return;

        actions.sort((a, b) => a.definition.display.name < b.definition.display.name ? -1 : a.definition.display.name === b.definition.display.name ? 0 : 1);

        return [
            ActionMenu.Header('Apply Action'),
            ...actions.map(a => ActionMenu.Item(a.definition.display.name, () => this.setState({ action: 'apply', currentAction: a })))
        ];
    }

    selectAction: ActionMenu.OnSelect = item => {
        if (!item) return;
        (item?.value as any)();
    }

    updates(margin: string) {
        const cell = this.props.cell;
        const decoratorChain = StateTreeSpine.getDecoratorChain(cell.parent!, cell.transform.ref);

        const decorators = [];
        for (let i = decoratorChain.length - 1; i >= 0; i--) {
            const d = decoratorChain[i];
            decorators!.push(<UpdateTransformControl key={`${d.transform.transformer.id}-${i}`} state={cell.parent!} transform={d.transform} noMargin wrapInExpander expanderHeaderLeftMargin={margin} />);
        }

        return <div className='msp-tree-updates-wrapper'>{decorators}</div>;
    }

    render() {
        const cell = this.props.cell;
        const n = cell.transform;
        if (!cell) return null;

        const isCurrent = this.is(cell.parent!.behaviors.currentObject.value);

        const disabled = cell.status !== 'error' && cell.status !== 'ok';

        let label: React.ReactNode;
        if (cell.status === 'error' || !cell.obj) {
            const name = cell.status === 'error' ? cell.errorText : n.transformer.definition.display.name;
            label = <Button className='msp-btn-tree-label msp-no-hover-outline' noOverflow title={name} onClick={this.state.isCurrent ? this.setCurrentRoot : this.setCurrent} disabled={disabled}>
                {cell.status === 'error' && <b>[{cell.status}]</b>} <span>{name}</span>
            </Button>;
        } else {
            const obj = cell.obj as PluginStateObject.Any;
            const title = `${obj.label} ${obj.description ? obj.description : ''}`;
            label = <Button className={`msp-btn-tree-label msp-type-class-${obj.type.typeClass}`} noOverflow disabled={disabled} title={title} onClick={this.state.isCurrent ? this.setCurrentRoot : this.setCurrent}>
                <span>{obj.label}</span> {obj.description ? <small>{obj.description}</small> : void 0}
            </Button>;
        }

        const children = cell.parent!.tree.children.get(this.ref);
        const cellState = cell.state;

        const expand = <IconButton svg={cellState.isCollapsed ? ArrowRightSvg : ArrowDropDownSvg} flex='20px' disabled={disabled} onClick={this.toggleExpanded} transparent className='msp-no-hover-outline' style={{ visibility: children.size > 0 ? 'visible' : 'hidden' }} />;
        const remove = !cell.state.isLocked ? <IconButton svg={DeleteOutlinedSvg} onClick={this.remove} disabled={disabled} small toggleState={false} /> : void 0;
        const visibility = <IconButton svg={cellState.isHidden ? VisibilityOffOutlinedSvg : VisibilityOutlinedSvg} toggleState={false} disabled={disabled} small onClick={this.toggleVisible} />;

        const marginStyle: React.CSSProperties = {
            marginLeft: `${this.props.depth * 8}px`
        };

        const row = <div className={`msp-flex-row msp-tree-row${isCurrent ? ' msp-tree-row-current' : ''}`} onMouseEnter={this.highlight} onMouseLeave={this.clearHighlight} style={marginStyle}>
            {expand}
            {label}
            {remove}
            {visibility}
        </div>;

        if (!isCurrent) return row;

        if (this.state.action === 'apply' && this.state.currentAction) {
            return <div style={{ marginBottom: '1px' }}>
                {row}
                <ControlGroup header={`Apply ${this.state.currentAction.definition.display.name}`} initialExpanded={true} hideExpander={true} hideOffset={false} onHeaderClick={this.hideApply} topRightIcon={CloseSvg} headerLeftMargin={`${this.props.depth * 8 + 21}px`}>
                    <ApplyActionControl onApply={this.hideApply} state={this.props.cell.parent!} action={this.state.currentAction} nodeRef={this.props.cell.transform.ref} hideHeader noMargin />
                </ControlGroup>
            </div>;
        }

        if (this.state.action === 'options') {
            const actions = this.actions;
            const updates = this.updates(`${this.props.depth * 8 + 21}px`);
            return <div style={{ marginBottom: '1px' }}>
                {row}
                {updates}
                {actions && <div style={{ marginLeft: `${this.props.depth * 8 + 21}px`, marginTop: '-1px' }}>
                    <ActionMenu items={actions} onSelect={this.selectAction} />
                </div>}
            </div>;
        }

        return row;
    }
}