/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import { Icon, IconName } from './icons';
import { ParamDefinition } from '../../mol-util/param-definition';
import { ControlGroup } from './common';

export class ActionMenu extends React.PureComponent<ActionMenu.Props> {
    hide = () => this.props.onSelect(void 0)

    render() {
        const cmd = this.props;

        return <div className='msp-action-menu-options' style={{ marginTop: cmd.header ? void 0 : '1px' }}>
            {/* {cmd.header && <div className='msp-control-group-header' style={{ position: 'relative' }}>
                <button className='msp-btn msp-btn-block' onClick={this.hide}>
                    <Icon name='off' style={{ position: 'absolute', right: '2px', top: 0 }} />
                    <b>{cmd.header}</b>
                </button>
            </div>} */}
            {cmd.header && <ControlGroup header={cmd.header} initialExpanded={true} hideExpander={true} hideOffset={false} onHeaderClick={this.hide} topRightIcon='off'></ControlGroup>}
            <Section items={cmd.items} onSelect={cmd.onSelect} current={cmd.current} multiselect={this.props.multiselect} />
        </div>
    }
}

export namespace ActionMenu {
    export type Props = { items: Items, onSelect: OnSelect | OnSelectMany, header?: string, current?: Item, multiselect?: boolean }

    export type OnSelect = (item: Item | undefined) => void
    export type OnSelectMany = (itemOrItems: Item[] | undefined) => void

    export type Items =  string | Item | Items[]
    export type Item = { label: string, icon?: IconName, disabled?: boolean, selected?: boolean, value: unknown }

    export function Item(label: string, value: unknown): Item
    export function Item(label: string, icon: string, value: unknown): Item
    export function Item(label: string, iconOrValue: any, value?: unknown): Item {
        if (value) return { label, icon: iconOrValue, value };
        return { label, value: iconOrValue };
    }

    export function createItems<T>(xs: ArrayLike<T>, options?: { filter?: (t: T) => boolean, label?: (t: T) => string, value?: (t: T) => any, category?: (t: T) => string | undefined }) {
        const { label, value, category } = options || { };
        let cats: Map<string, (ActionMenu.Item | string)[]> | undefined = void 0;
        const items: (ActionMenu.Item | (ActionMenu.Item | string)[] | string)[] = [];
        for (let i = 0; i < xs.length; i++) {
            const x = xs[i];

            if (options?.filter && !options.filter(x)) continue;

            const catName = category?.(x);
            const l = label ? label(x) : '' + x;
            const v = value ? value(x) : x;

            if (!!catName) {
                if (!cats) cats = new Map<string, (ActionMenu.Item | string)[]>();

                let cat = cats.get(catName);
                if (!cat) {
                    cat = [catName];
                    cats.set(catName, cat);
                    items.push(cat);
                }
                cat.push(ActionMenu.Item(l, v));
            } else {
                items.push(ActionMenu.Item(l, v));
            }
        }
        return items as ActionMenu.Items;
    }
    
    type Opt = ParamDefinition.Select<any>['options'][0];
    const _selectOptions = { value: (o: Opt) => o[0], label: (o: Opt) => o[1], category: (o: Opt) => o[2] };

    export function createItemsFromSelectParam(param: ParamDefinition.Select<any>) {
        return createItems(param.options, _selectOptions);
    }

    export function hasSelectedItem(items: Items): boolean {
        if (typeof items === 'string') return false;
        if (isItem(items)) return !!items.selected;
        for (const s of items) {
            const found = hasSelectedItem(s);
            if (found) return true;
        }
        return false;
    }

    export function findItem(items: Items, value: any): Item | undefined {
        if (typeof items === 'string') return;
        if (isItem(items)) return items.value === value ? items : void 0;
        for (const s of items) {
            const found = findItem(s, value);
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

type SectionProps = { header?: string, items: ActionMenu.Items, onSelect: ActionMenu.OnSelect | ActionMenu.OnSelectMany, current: ActionMenu.Item | undefined, multiselect: boolean | undefined }
type SectionState = { items: ActionMenu.Items, current: ActionMenu.Item | undefined, isExpanded: boolean }

class Section extends React.PureComponent<SectionProps, SectionState> {
    state = {
        items: this.props.items,
        current: this.props.current,
        isExpanded: this.hasCurrent
    }

    get hasCurrent() {
        return this.props.multiselect
            ? ActionMenu.hasSelectedItem(this.props.items)
            : !!this.props.current && !!ActionMenu.findItem(this.props.items, this.props.current.value);
    }

    toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isExpanded: !this.state.isExpanded });
        e.currentTarget.blur();
    }

    static getDerivedStateFromProps(props: SectionProps, state: SectionState) {
        if (props.items === state.items && props.current === state.current) return null;
        return {
            items: props.items,
            current: props.current,
            isExpanded: props.multiselect
                ? ActionMenu.hasSelectedItem(props.items)
                : props.current && !!ActionMenu.findItem(props.items, props.current.value) }
    }

    selectAll = () => {
        const items = collectItems(this.props.items, []).filter(i => !i.selected);
        this.props.onSelect(items as any);
    }

    selectNone = () => {
        const items = collectItems(this.props.items, []).filter(i => !!i.selected);
        this.props.onSelect(items as any);
    }

    get multiselectHeader() {
        const { header } = this.props;
        const hasCurrent = this.hasCurrent;

        return <div className='msp-control-group-header msp-flex-row' style={{ marginTop: '1px' }}>
            <button className='msp-btn msp-form-control msp-flex-item' onClick={this.toggleExpanded}>
                <Icon name={this.state.isExpanded ? 'collapse' : 'expand'} />
                {hasCurrent ? <b>{header}</b> : header}
            </button>
            <button className='msp-btn msp-form-control msp-flex-item' onClick={this.selectAll} style={{ flex: '0 0 50px', textAlign: 'right' }}>
                <Icon name='check' />
                All
            </button>
            <button className='msp-btn msp-form-control msp-flex-item' onClick={this.selectNone} style={{ flex: '0 0 50px', textAlign: 'right' }}>
                <Icon name='cancel' />
                None
            </button>
        </div>;
    }

    get basicHeader() {
        const { header } = this.props;
        const hasCurrent = this.hasCurrent;

        return <div className='msp-control-group-header' style={{ marginTop: '1px' }}>
            <button className='msp-btn msp-btn-block msp-form-control' onClick={this.toggleExpanded}>
                <Icon name={this.state.isExpanded ? 'collapse' : 'expand'} />
                {hasCurrent ? <b>{header}</b> : header}
            </button>
        </div>;
    }

    render() {
        const { header, items, onSelect, current } = this.props;

        if (typeof items === 'string') return null;
        if (isItem(items)) return <Action item={items} onSelect={onSelect} current={current} multiselect={this.props.multiselect} />

        return <div>
            {header && (this.props.multiselect && this.state.isExpanded ? this.multiselectHeader : this.basicHeader)}
            <div className='msp-control-offset'>
                {(!header || this.state.isExpanded) && items.map((x, i) => {
                    if (typeof x === 'string') return null;
                    if (isItem(x)) return <Action key={i} item={x} onSelect={onSelect} current={current} multiselect={this.props.multiselect} />
                    return <Section key={i} header={typeof x[0] === 'string' ? x[0] : void 0} items={x} onSelect={onSelect} current={current} multiselect={this.props.multiselect} />
                })}
            </div>
        </div>;
    }
}

const Action: React.FC<{
    item: ActionMenu.Item, 
    onSelect: ActionMenu.OnSelect | ActionMenu.OnSelectMany, 
    multiselect: boolean | undefined, 
    current: ActionMenu.Item | undefined }> = ({ item, onSelect, current, multiselect }) => {
    const isCurrent = current === item;
    return <div className='msp-control-row'>
        <button onClick={() => onSelect(multiselect ? [item] : item as any)} disabled={item.disabled}>
            {item.icon && <Icon name={item.icon} style={{ fontSize: '80%', marginRight: '6px' }} />}
            {isCurrent || item.selected ? <b>{item.label}</b> : item.label}
        </button>
    </div>;
}

function isItem(x: any): x is ActionMenu.Item {
    const v = x as ActionMenu.Item;
    return v && !!v.label && typeof v.value !== 'undefined';
}

function collectItems(items: ActionMenu.Items, target: ActionMenu.Item[]) {
    if (typeof items === 'string') return target;
    if (isItem(items)) {
        target.push(items);
        return target;
    }
    for (const i of items) {
        collectItems(i, target);
    }
    return target;
}