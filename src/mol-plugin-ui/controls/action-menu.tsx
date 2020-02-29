/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import { Icon } from './common';
import { ParamDefinition } from '../../mol-util/param-definition';

export class ActionMenu extends React.PureComponent<ActionMenu.Props> {
    hide = () => this.props.onSelect(void 0)

    render() {
        const cmd = this.props;

        return <div className='msp-action-menu-options' style={{ marginTop: '1px' }}>
            {cmd.header && <div className='msp-control-group-header' style={{ position: 'relative' }}>
                <button className='msp-btn msp-btn-block' onClick={this.hide}>
                    <Icon name='off' style={{ position: 'absolute', right: '2px', top: 0 }} />
                    <b>{cmd.header}</b>
                </button>
            </div>}
            <Section items={cmd.items} onSelect={cmd.onSelect} current={cmd.current} />
        </div>
    }
}

export namespace ActionMenu {
    export type Props = { items: Items, onSelect: OnSelect, header?: string, current?: Item | undefined }

    export type OnSelect = (item: Item | undefined) => void

    export type Items = string | Item | [Items]
    export type Item = { label: string, icon?: string, value: unknown }

    export function Item(label: string, value: unknown): Item
    export function Item(label: string, icon: string, value: unknown): Item
    export function Item(label: string, iconOrValue: any, value?: unknown): Item {
        if (value) return { label, icon: iconOrValue, value };
        return { label, value: iconOrValue };
    }

    function createSpecFromSelectParamSimple(param: ParamDefinition.Select<any>) {
        const items: Item[] = [];
        for (const [v, l] of param.options) {
            items.push(ActionMenu.Item(l, v));
        }
        return items as Items;
    }

    function createSpecFromSelectParamCategories(param: ParamDefinition.Select<any>) {
        const cats = new Map<string, (Item | string)[]>();
        const items: (Item | (Item | string)[] | string)[] = [];
        for (const [v, l, c] of param.options) {
            if (!!c) {
                let cat = cats.get(c);
                if (!cat) {
                    cat = [c];
                    cats.set(c, cat);
                    items.push(cat);
                }
                cat.push(ActionMenu.Item(l, v));
            } else {
                items.push(ActionMenu.Item(l, v));
            }
        }
        return items as Items;
    }

    export function createSpecFromSelectParam(param: ParamDefinition.Select<any>) {
        for (const o of param.options) {
            if (!!o[2]) return createSpecFromSelectParamCategories(param);
        }
        return createSpecFromSelectParamSimple(param);
    }

    export function findCurrent(items: Items, value: any): Item | undefined {
        if (typeof items === 'string') return;
        if (isItem(items)) return items.value === value ? items : void 0;
        for (const s of items) {
            const found = findCurrent(s, value);
            if (found) return found;
        }
    }

    export function getFirstItem(items: Items): Item | undefined {
        if (typeof items === 'string') return;
        if (isItem(items)) return items;
        for (const s of items) {
            const found = getFirstItem(s);
            if (found) return found;
        }
    }
}

type SectionProps = { header?: string, items: ActionMenu.Items, onSelect: ActionMenu.OnSelect, current: ActionMenu.Item | undefined }
type SectionState = { items: ActionMenu.Items, current: ActionMenu.Item | undefined, isExpanded: boolean }

class Section extends React.PureComponent<SectionProps, SectionState> {
    state = {
        items: this.props.items,
        current: this.props.current,
        isExpanded: !!this.props.current && !!ActionMenu.findCurrent(this.props.items, this.props.current.value)
    }

    toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isExpanded: !this.state.isExpanded });
        e.currentTarget.blur();
    }

    static getDerivedStateFromProps(props: SectionProps, state: SectionState) {
        if (props.items === state.items && props.current === state.current) return null;
        return { items: props.items, current: props.current, isExpanded: props.current && !!ActionMenu.findCurrent(props.items, props.current.value) }
    }

    render() {
        const { header, items, onSelect, current } = this.props;

        if (typeof items === 'string') return null;
        if (isItem(items)) return <Action item={items} onSelect={onSelect} current={current} />

        const hasCurrent = header && current && !!ActionMenu.findCurrent(items, current.value)

        return <div>
            {header && <div className='msp-control-group-header' style={{ marginTop: '1px' }}>
                <button className='msp-btn msp-btn-block' onClick={this.toggleExpanded}>
                    <span className={`msp-icon msp-icon-${this.state.isExpanded ? 'collapse' : 'expand'}`} />
                    {hasCurrent ? <b>{header}</b> : header}
                </button>
            </div>}
            <div className='msp-control-offset'>
                {(!header || this.state.isExpanded) && items.map((x, i) => {
                    if (typeof x === 'string') return null;
                    if (isItem(x)) return <Action key={i} item={x} onSelect={onSelect} current={current} />
                    return <Section key={i} header={typeof x[0] === 'string' ? x[0] : void 0} items={x} onSelect={onSelect} current={current} />
                })}
            </div>
        </div>;
    }
}

const Action: React.FC<{ item: ActionMenu.Item, onSelect: ActionMenu.OnSelect, current: ActionMenu.Item | undefined }> = ({ item, onSelect, current }) => {
    const isCurrent = current === item;
    return <div className='msp-control-row'>
        <button onClick={() => onSelect(item)}>
            {item.icon && <Icon name={item.icon} />}
            {isCurrent ? <b>{item.label}</b> : item.label}
        </button>
    </div>;
}

function isItem(x: any): x is ActionMenu.Item {
    const v = x as ActionMenu.Item;
    return v && !!v.label && typeof v.value !== 'undefined';
}