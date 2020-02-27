/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import { Icon } from './common';
import { Subscription, BehaviorSubject, Observable } from 'rxjs';
import { ParamDefinition } from '../../mol-util/param-definition';

export class ActionMenu {
    private _command: BehaviorSubject<ActionMenu.Command>;

    get commands(): Observable<ActionMenu.Command> { return this._command; }

    hide() {
        this._command.next(HideCmd)
    }

    toggle(params: { items: ActionMenu.Spec, header?: string, current?: ActionMenu.Item, onSelect: (value: any) => void }) {
        this._command.next({ type: 'toggle', ...params });
    }

    constructor(defaultCommand?: ActionMenu.Command) {
        this._command = new BehaviorSubject<ActionMenu.Command>(defaultCommand || { type: 'hide' });
    }
}

const HideCmd: ActionMenu.Command = { type: 'hide' };

export namespace ActionMenu {
    export type Command = 
        | { type: 'toggle', items: Spec, header?: string, current?: Item, onSelect: (value: any) => void }
        | { type: 'hide' }

    function isToggleOff(a: Command, b: Command) {
        if (a.type === 'hide' || b.type === 'hide') return false;
        return a.onSelect === b.onSelect && a.items === b.items;
    }


    export type ToggleProps = {
        style?: React.HTMLAttributes<HTMLButtonElement>,
        className?: string,
        menu: ActionMenu,
        disabled?: boolean,
        items: ActionMenu.Spec,
        header?: string,
        label?: string,
        current?: ActionMenu.Item,
        onSelect: (value: any) => void
    }

    export class Toggle extends React.PureComponent<ToggleProps, { isSelected: boolean }> {
        private sub: Subscription | undefined = void 0;

        state = { isSelected: false };

        componentDidMount() {
            this.sub = this.props.menu.commands.subscribe(command => {
                if (command.type === 'hide') {
                    this.hide();
                } else if (command.type === 'toggle') {
                    const cmd = this.props;
                    if (command.items === cmd.items && command.onSelect === cmd.onSelect) {
                        this.setState({ isSelected: !this.state.isSelected });
                    } else {
                        this.hide();
                    }
                }
            });
        }

        componentWillUnmount() {
            if (!this.sub) return;
            this.sub.unsubscribe();
            this.sub = void 0;
        }

        hide = () => this.setState({ isSelected: false });

        onClick = (e: React.MouseEvent<HTMLButtonElement>) => {
            e.currentTarget.blur();
            this.props.menu.toggle(this.props);
        }

        render() {
            const props = this.props;
            const label = props.label || props.header;
            return <button onClick={this.onClick} 
                disabled={props.disabled} style={props.style} className={props.className}>
                    {this.state.isSelected ? <b>{label}</b> : label}
            </button>;
        }
    }

    export class Options extends React.PureComponent<{ menu: ActionMenu }, { command: Command, isVisible: boolean }> {
        private sub: Subscription | undefined = void 0;

        state = { isVisible: false, command: HideCmd };

        componentDidMount() {
            this.sub = this.props.menu.commands.subscribe(command => {
                if (command.type === 'hide' || isToggleOff(command, this.state.command)) {
                    this.setState({ isVisible: false, command: HideCmd });
                } else {
                    this.setState({ isVisible: true, command })
                }
            });
        }

        componentWillUnmount() {
            if (!this.sub) return;
            this.sub.unsubscribe();
            this.sub = void 0;
        }

        onSelect: OnSelect = item => {
            const cmd = this.state.command;
            this.hide();
            if (cmd.type === 'toggle') cmd.onSelect(item.value);
        }

        hide = () => {
            this.props.menu.hide();
        }

        render() {
            if (!this.state.isVisible || this.state.command.type !== 'toggle') return null;
            return <div className='msp-action-menu-options' style={{ marginTop: '1px' }}>
                {this.state.command.header && <div className='msp-control-group-header' style={{ position: 'relative' }}>
                    <button className='msp-btn msp-btn-block' onClick={this.hide}>
                        <Icon name='off' style={{ position: 'absolute', right: '2px', top: 0 }} />
                        <b>{this.state.command.header}</b>
                    </button>
                </div>}
                <Section menu={this.props.menu} items={this.state.command.items} onSelect={this.onSelect} current={this.state.command.current} />
            </div>
        }
    }

    type SectionProps = { menu: ActionMenu, header?: string, items: Spec, onSelect: OnSelect, current: Item | undefined }
    type SectionState = { items: Spec, current: Item | undefined, isExpanded: boolean }

    class Section extends React.PureComponent<SectionProps, SectionState> {
        state = {
            items: this.props.items,
            current: this.props.current,
            isExpanded: !!this.props.current && !!findCurrent(this.props.items, this.props.current.value)
        }

        toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
            this.setState({ isExpanded: !this.state.isExpanded });
            e.currentTarget.blur();
        }

        static getDerivedStateFromProps(props: SectionProps, state: SectionState) {
            if (props.items === state.items && props.current === state.current) return null;
            return { items: props.items, current: props.current, isExpanded: props.current && !!findCurrent(props.items, props.current.value) }
        }

        render() {
            const { header, items, onSelect, current, menu } = this.props;

            if (typeof items === 'string') return null;
            if (isItem(items)) return <Action menu={menu} item={items} onSelect={onSelect} current={current} />
            return <div>
                {header && <div className='msp-control-group-header' style={{ marginTop: '1px' }}>
                    <button className='msp-btn msp-btn-block' onClick={this.toggleExpanded}>
                        <span className={`msp-icon msp-icon-${this.state.isExpanded ? 'collapse' : 'expand'}`} />
                        {header}
                    </button>
                </div>}
                <div className='msp-control-offset'>
                    {(!header || this.state.isExpanded) && items.map((x, i) => {
                        if (typeof x === 'string') return null;
                        if (isItem(x)) return <Action menu={menu} key={i} item={x} onSelect={onSelect} current={current} />
                        return <Section menu={menu} key={i} header={typeof x[0] === 'string' ? x[0] : void 0} items={x} onSelect={onSelect} current={current} />
                    })}
                </div>
            </div>;
        }
    }

    const Action: React.FC<{ menu: ActionMenu, item: Item, onSelect: OnSelect, current: Item | undefined }> = ({ menu, item, onSelect, current }) => {
        const isCurrent = current === item;
        return <div className='msp-control-row'>
            <button onClick={isCurrent ? () => menu.hide() : () => onSelect(item)}>
                {item.icon && <Icon name={item.icon} />}
                {isCurrent ? <b>{item.name}</b> : item.name}
            </button>
        </div>;
    }

    type OnSelect = (item: Item) => void

    function isItem(x: any): x is Item {
        const v = x as Item;
        return v && !!v.name && typeof v.value !== 'undefined';
    }

    export type OptionsParams = { items: Spec, header?: string, onSelect: (value: any) => void }
    export type Spec = string | Item | [Spec]
    export type Item = { name: string, icon?: string, value: unknown }

    export function Item(name: string, value: unknown): Item
    export function Item(name: string, icon: string, value: unknown): Item
    export function Item(name: string, iconOrValue: any, value?: unknown): Item {
        if (value) return { name, icon: iconOrValue, value };
        return { name, value: iconOrValue };
    }

    function createSpecFromSelectParamSimple(param: ParamDefinition.Select<any>) {
        const spec: Item[] = [];
        for (const [v, l] of param.options) {
            spec.push(ActionMenu.Item(l, v));
        }
        return spec as Spec;
    }

    function createSpecFromSelectParamCategories(param: ParamDefinition.Select<any>) {
        const cats = new Map<string, (Item | string)[]>();
        const spec: (Item | (Item | string)[] | string)[] = [];
        for (const [v, l, c] of param.options) {
            if (!!c) {
                let cat = cats.get(c);
                if (!cat) {
                    cat = [c];
                    cats.set(c, cat);
                    spec.push(cat);
                }
                cat.push(ActionMenu.Item(l, v));
            } else {
                spec.push(ActionMenu.Item(l, v));
            }
        }
        return spec as Spec;
    }

    export function createSpecFromSelectParam(param: ParamDefinition.Select<any>) {
        for (const o of param.options) {
            if (!!o[2]) return createSpecFromSelectParamCategories(param);
        }
        return createSpecFromSelectParamSimple(param);
    }

    export function findCurrent(spec: Spec, value: any): Item | undefined {
        if (typeof spec === 'string') return;
        if (isItem(spec)) return spec.value === value ? spec : void 0;
        for (const s of spec) {
            const found = findCurrent(s, value);
            if (found) return found;
        }
    }
}