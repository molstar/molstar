/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import { Icon } from './common';
import { Subscription, BehaviorSubject, Observable } from 'rxjs';

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
        header: string,
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

        render() {
            const props = this.props;
            return <button onClick={() => props.menu.toggle(props)} 
                disabled={props.disabled} style={props.style} className={props.className}>
                    {this.state.isSelected ? <b>{props.header}</b> : props.header}
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
            return <div className='msp-action-menu-options'>
                {this.state.command.header && <div className='msp-control-group-header' style={{ position: 'relative' }}>
                    <button className='msp-btn msp-btn-block' onClick={this.hide}>
                        <Icon name='off' style={{ position: 'absolute', right: '2px', top: 0 }} />
                        <b>{this.state.command.header}</b>
                    </button>
                </div>}
                <Section items={this.state.command.items} onSelect={this.onSelect} />
            </div>
        }
    }

    class Section extends React.PureComponent<{ header?: string, items: Spec, onSelect: OnSelect }, { isExpanded: boolean }> {
        state = { isExpanded: false }

        toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
            this.setState({ isExpanded: !this.state.isExpanded });
            e.currentTarget.blur();
        }

        render() {
            const { header, items, onSelect } = this.props;
            if (typeof items === 'string') return null;
            if (isItem(items)) return <Action item={items} onSelect={onSelect} />
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
                        if (isItem(x)) return <Action key={i} item={x} onSelect={onSelect} />
                        return <Section key={i} header={typeof x[0] === 'string' ? x[0] : void 0} items={x} onSelect={onSelect} />
                    })}
                </div>
            </div>;
        }
    }

    const Action: React.FC<{ item: Item, onSelect: OnSelect }> = ({ item, onSelect }) => {
        return <div className='msp-control-row'>
            <button onClick={() => onSelect(item)}>
                {item.icon && <Icon name={item.icon} />}
                {item.name}
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
}