/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import { Icon } from './common';
import { Observable, Subscription } from 'rxjs';

export namespace ActionMenu {
    export class Options extends React.PureComponent<{ toggle: Observable<OptionsParams | undefined>, hide?: Observable<any> }, { options: OptionsParams | undefined, isVisible: boolean }> { 
        private subs: Subscription[] = [];

        state = { isVisible: false, options: void 0 as OptionsParams | undefined }

        componentDidMount() {
            this.subs.push(this.props.toggle.subscribe(options => {
                if (options && this.state.options?.items === options.items && this.state.options?.onSelect === options.onSelect) {
                    this.setState({ isVisible: !this.state.isVisible});
                } else {
                    this.setState({ isVisible: !!options, options: options })
                }
            }));

            if (this.props.hide) {
                this.subs.push(this.props.hide.subscribe(() => this.hide()));
            }
        }

        componentWillUnmount() {
            if (!this.subs) return;
            for (const s of this.subs) s.unsubscribe();
            this.subs = [];
        }

        onSelect: OnSelect = item => {
            this.setState({ isVisible: false, options: void 0 });
            this.state.options?.onSelect(item.value);
        }

        hide = () => this.setState({ isVisible: false, options: void 0 });

        render() {
            if (!this.state.isVisible || !this.state.options) return null;
            return <div className='msp-action-menu-options'>
                {this.state.options.header && <div className='msp-control-group-header'>
                    <button className='msp-btn msp-btn-block' onClick={this.hide}>
                        {this.state.options.header}
                    </button>
                </div>}
                <Section items={this.state.options!.items} onSelect={this.onSelect} />
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
            return <div className='msp-control-offset'>
                {header && <div className='msp-control-group-header'>
                    <button className='msp-btn msp-btn-block' onClick={this.toggleExpanded}>
                        <span className={`msp-icon msp-icon-${this.state.isExpanded ? 'collapse' : 'expand'}`} /> 
                        {header}
                    </button>
                </div>}
                {(!header || this.state.isExpanded) && items.map((x, i) => {
                    if (typeof x === 'string') return null;
                    if (isItem(x)) return <Action key={i} item={x} onSelect={onSelect} />
                    return <Section key={i} header={typeof x[0] === 'string' ? x[0] : void 0} items={x} onSelect={onSelect} />
                })}
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