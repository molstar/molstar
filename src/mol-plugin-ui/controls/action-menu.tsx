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
        return <div className='msp-action-menu-options' style={{ /* marginTop: cmd.header ? void 0 : '1px', */ maxHeight: '300px', overflow: 'hidden', overflowY: 'auto' }}>
            {cmd.header && <ControlGroup header={cmd.header} initialExpanded={true} hideExpander={true} hideOffset={false} onHeaderClick={this.hide} topRightIcon='off'></ControlGroup>}
            <Section items={cmd.items} onSelect={cmd.onSelect} current={cmd.current} multiselect={this.props.multiselect} noOffset={this.props.noOffset} />
        </div>
    }
}

export namespace ActionMenu {
    export type Props = { items: Items, onSelect: OnSelect | OnSelectMany, header?: string, current?: Item, multiselect?: boolean, noOffset?: boolean }

    export type OnSelect = (item: Item | undefined) => void
    export type OnSelectMany = (itemOrItems: Item[] | undefined) => void

    export type Items =  Header | Item | Items[]
    export type Header = { kind: 'header', label: string, isIndependent?: boolean, initiallyExpanded?: boolean }
    export type Item = { kind: 'item', label: string, icon?: IconName, disabled?: boolean, selected?: boolean, value: unknown, addOn?: JSX.Element }

    export function Header(label: string, options?: { isIndependent?: boolean, initiallyExpanded?: boolean }): Header {
        return options ? { kind: 'header', label, ...options } : { kind: 'header', label };
    }

    export function Item(label: string, value: unknown): Item
    export function Item(label: string, icon: IconName, value: unknown): Item
    export function Item(label: string, iconOrValue: any, value?: unknown): Item {
        if (value) return { kind: 'item', label, icon: iconOrValue, value };
        return { kind: 'item', label, value: iconOrValue };
    }

    export interface CreateItemsParams<T> {
        filter?: (t: T) => boolean,
        label?: (t: T) => string,
        value?: (t: T) => any,
        category?: (t: T) => string | undefined,
        icon?: (t: T) => IconName | undefined,
        selected?: (t: T) => boolean | undefined,
        addOn?: (t: T) => JSX.Element | undefined
    }

    export function createItems<T>(xs: ArrayLike<T>, params?: CreateItemsParams<T>) {
        const { label, value, category, selected, icon, addOn } = params || { };
        let cats: Map<string, (ActionMenu.Item | ActionMenu.Header)[]> | undefined = void 0;
        const items: (ActionMenu.Item | (ActionMenu.Item | ActionMenu.Header)[] | string)[] = [];
        for (let i = 0; i < xs.length; i++) {
            const x = xs[i];

            if (params?.filter && !params.filter(x)) continue;

            const catName = category?.(x);
            const l = label ? label(x) : '' + x;
            const v = value ? value(x) : x;

            let cat: (ActionMenu.Item | ActionMenu.Header)[] | undefined;
            if (!!catName) {
                if (!cats) cats = new Map<string, (ActionMenu.Item | ActionMenu.Header)[]>();

                cat = cats.get(catName);
                if (!cat) {
                    cat = [ActionMenu.Header(catName)];
                    cats.set(catName, cat);
                    items.push(cat);
                }
            } else {
                cat = items as any;
            }

            const ao = addOn?.(x);

            cat!.push({ kind: 'item', label: l, value: v, icon: icon ? icon(x) : void 0, selected: selected ? selected(x) : void 0, addOn: ao });
        }
        return items as ActionMenu.Items;
    }

    type Opt = ParamDefinition.Select<any>['options'][0];
    const _selectOptions = { value: (o: Opt) => o[0], label: (o: Opt) => o[1], category: (o: Opt) => o[2] };

    export function createItemsFromSelectOptions<O extends ParamDefinition.Select<any>['options']>(options: O, params?: CreateItemsParams<O[0]>) {
        return createItems(options, params ? { ..._selectOptions, ...params } : _selectOptions);
    }

    export function hasSelectedItem(items: Items): boolean {
        if (isHeader(items)) return false;
        if (isItem(items)) return !!items.selected;
        for (const s of items) {
            const found = hasSelectedItem(s);
            if (found) return true;
        }
        return false;
    }

    export function findItem(items: Items, value: any): Item | undefined {
        if (isHeader(items)) return;
        if (isItem(items)) return items.value === value ? items : void 0;
        for (const s of items) {
            const found = findItem(s, value);
            if (found) return found;
        }
    }

    export function getFirstItem(items: Items): Item | undefined {
        if (isHeader(items)) return;
        if (isItem(items)) return items;
        for (const s of items) {
            const found = getFirstItem(s);
            if (found) return found;
        }
    }
}

type SectionProps = { items: ActionMenu.Items, onSelect: ActionMenu.OnSelect | ActionMenu.OnSelectMany, current: ActionMenu.Item | undefined, multiselect: boolean | undefined, noOffset?: boolean }
type SectionState = { items: ActionMenu.Items, current: ActionMenu.Item | undefined, isExpanded: boolean, hasCurrent: boolean, header?: ActionMenu.Header }

class Section extends React.PureComponent<SectionProps, SectionState> {
    static createState(props: SectionProps): SectionState {
        const header = isItems(props.items) && isHeader(props.items[0]) ? props.items[0] : void 0;

        const hasCurrent = header?.isIndependent
            ? false
            : props.multiselect
                ? ActionMenu.hasSelectedItem(props.items)
                : (!!props.current && !!ActionMenu.findItem(props.items, props.current.value)) || ActionMenu.hasSelectedItem(props.items);

        return {
            items: props.items,
            current: props.current,
            header,
            hasCurrent,
            isExpanded: hasCurrent || !!header?.initiallyExpanded
        };
    }

    state = Section.createState(this.props)

    toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isExpanded: !this.state.isExpanded });
        e.currentTarget.blur();
    }

    static getDerivedStateFromProps(props: SectionProps, state: SectionState) {
        if (props.items === state.items && props.current === state.current) return null;
        return Section.createState(props);
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
        const { header, hasCurrent } = this.state;

        return <div className='msp-control-group-header msp-flex-row' style={{ marginTop: '1px' }}>
            <button className='msp-btn msp-form-control msp-flex-item msp-no-overflow' onClick={this.toggleExpanded}>
                <Icon name={this.state.isExpanded ? 'collapse' : 'expand'} />
                {hasCurrent ? <b>{header?.label}</b> : header?.label}
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
        const { header, hasCurrent } = this.state;

        return <div className='msp-control-group-header' style={{ marginTop: '1px' }}>
            <button className='msp-btn msp-btn-block msp-form-control' onClick={this.toggleExpanded}>
                <Icon name={this.state.isExpanded ? 'collapse' : 'expand'} />
                {hasCurrent ? <b>{header?.label}</b> : header?.label}
            </button>
        </div>;
    }

    render() {
        const { items, onSelect, current } = this.props;

        if (isHeader(items)) return null;
        if (isItem(items)) return <Action item={items} onSelect={onSelect} current={current} multiselect={this.props.multiselect} />

        const { header } = this.state;
        return <>
            {header && (this.props.multiselect && this.state.isExpanded ? this.multiselectHeader : this.basicHeader)}
            <div className={this.props.noOffset ? void 0 : 'msp-control-offset'}>
                {(!header || this.state.isExpanded) && items.map((x, i) => {
                    if (isHeader(x)) return null;
                    if (isItem(x)) return <Action key={i} item={x} onSelect={onSelect} current={current} multiselect={this.props.multiselect} />
                    return <Section key={i} items={x} onSelect={onSelect} current={current} multiselect={this.props.multiselect} />
                })}
            </div>
        </>;
    }
}

const Action: React.FC<{
    item: ActionMenu.Item,
    onSelect: ActionMenu.OnSelect | ActionMenu.OnSelectMany,
    multiselect: boolean | undefined,
    current: ActionMenu.Item | undefined }> = ({ item, onSelect, current, multiselect }) => {
    const isCurrent = current === item;

    const style: React.CSSProperties | undefined = item.addOn ? { position: 'relative' } : void 0;

    return <button className='msp-btn msp-btn-block msp-form-control msp-action-menu-button msp-no-overflow' onClick={() => onSelect(multiselect ? [item] : item as any)} disabled={item.disabled} style={style}>
        {item.icon && <Icon name={item.icon} />}
        {isCurrent || item.selected ? <b>{item.label}</b> : item.label}
        {item.addOn}
    </button>;
}

function isItems(x: any): x is ActionMenu.Items[] {
    return !!x && Array.isArray(x);
}

function isItem(x: any): x is ActionMenu.Item {
    const v = x as ActionMenu.Item;
    return v && v.kind === 'item';
}

function isHeader(x: any): x is ActionMenu.Header {
    const v = x as ActionMenu.Header;
    return v && v.kind === 'header';
}

function collectItems(items: ActionMenu.Items, target: ActionMenu.Item[]) {
    if (isHeader(items)) return target;
    if (isItem(items)) {
        target.push(items);
        return target;
    }
    for (const i of items) {
        collectItems(i, target);
    }
    return target;
}