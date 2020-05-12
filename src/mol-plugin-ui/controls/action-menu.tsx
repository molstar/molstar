/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { ParamDefinition } from '../../mol-util/param-definition';
import { Button, ControlGroup } from './common';
import { CloseSvg, ArrowDropDownSvg, ArrowRightSvg, CheckSvg } from './icons';

export class ActionMenu extends React.PureComponent<ActionMenu.Props> {
    hide = () => this.props.onSelect(void 0)

    render() {
        const cmd = this.props;
        const section = <Section items={cmd.items} onSelect={cmd.onSelect} current={cmd.current} multiselect={this.props.multiselect} noOffset={this.props.noOffset} noAccent={this.props.noAccent} />;
        return <div className={`msp-action-menu-options${cmd.header ? '' : ' msp-action-menu-options-no-header'}`}>
            {cmd.header && <ControlGroup header={cmd.header} title={cmd.title} initialExpanded={true} hideExpander={true} hideOffset onHeaderClick={this.hide} topRightIcon={CloseSvg}>
                {section}
            </ControlGroup>}
            {!cmd.header && section}
        </div>;
    }
}

export namespace ActionMenu {
    export type Props = {
        items: Items,
        onSelect: OnSelect | OnSelectMany,
        header?: string,
        title?: string,
        current?: Item,
        multiselect?: boolean,
        noOffset?: boolean,
        noAccent?: boolean
    }

    export type OnSelect = (item: Item | undefined, e?: React.MouseEvent<HTMLButtonElement>) => void
    export type OnSelectMany = (itemOrItems: Item[] | undefined, e?: React.MouseEvent<HTMLButtonElement>) => void

    export type Items =  Header | Item | Items[]
    export type Header = { kind: 'header', label: string, isIndependent?: boolean, initiallyExpanded?: boolean, description?: string }
    export type Item = { kind: 'item', label: string, icon?: React.FC, disabled?: boolean, selected?: boolean, value: unknown, addOn?: JSX.Element, description?: string }

    export function Header(label: string, options?: { isIndependent?: boolean, initiallyExpanded?: boolean, description?: string }): Header {
        return options ? { kind: 'header', label, ...options } : { kind: 'header', label };
    }

    export function Item(label: string, value: unknown, options?: { icon?: React.FC, description?: string }): Item {
        return { kind: 'item', label, value, ...options };
    }

    export interface CreateItemsParams<T> {
        filter?: (t: T) => boolean,
        label?: (t: T) => string,
        value?: (t: T) => any,
        category?: (t: T) => string | undefined,
        icon?: (t: T) => React.FC | undefined,
        selected?: (t: T) => boolean | undefined,
        addOn?: (t: T) => JSX.Element | undefined
        description?: (t: T) => string | undefined,
    }

    export function createItems<T>(xs: ArrayLike<T>, params?: CreateItemsParams<T>): Items[] {
        const { label, value, category, selected, icon, addOn, description } = params || { };
        let cats: Map<string, (ActionMenu.Item | ActionMenu.Header)[]> | undefined = void 0;
        const items: (ActionMenu.Item | (ActionMenu.Item | ActionMenu.Header)[])[] = [];
        for (let i = 0; i < xs.length; i++) {
            const x = xs[i];

            if (params?.filter && !params.filter(x)) continue;

            const catName = category?.(x);
            const l = label ? label(x) : '' + x;
            const v = value ? value(x) : x;
            const d = description ? description(x) :
                typeof x === 'string' ? x : undefined;

            let cat: (ActionMenu.Item | ActionMenu.Header)[] | undefined;
            if (!!catName) {
                if (!cats) cats = new Map<string, (ActionMenu.Item | ActionMenu.Header)[]>();

                cat = cats.get(catName);
                if (!cat) {
                    cat = [ActionMenu.Header(catName, { description: catName })];
                    cats.set(catName, cat);
                    items.push(cat);
                }
            } else {
                cat = items as any;
            }

            const ao = addOn?.(x);

            cat!.push({ kind: 'item', label: l, value: v, icon: icon ? icon(x) : void 0, selected: selected ? selected(x) : void 0, addOn: ao, description: d });
        }
        return items;
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

    // export type SelectProps<T> = {
    //     items: Items,
    //     onSelect: (item: Item) => void,
    //     disabled?: boolean,
    //     label?: string,
    //     current?: Item,
    //     style?: React.CSSProperties
    // }

    // export class Select<T> extends React.PureComponent<SelectProps<T>, { isExpanded: boolean }> {
    //     state = { isExpanded: false };

    //     toggleExpanded = () => this.setState({ isExpanded: !this.state.isExpanded })
    //     onSelect: OnSelect = (item) => {
    //         this.setState({ isExpanded: false });
    //         if (!item) return;
    //         this.onSelect(item);
    //     }

    //     render() {
    //         const current = this.props.current;
    //         const label = this.props.label || current?.label || '';

    //         return <div className='msp-action-menu-select' style={this.props.style}>
    //             <ToggleButton disabled={this.props.disabled} style={{ textAlign: 'left' }} className='msp-no-overflow'
    //                 label={label} title={label as string} toggle={this.toggleExpanded} isSelected={this.state.isExpanded} />
    //             {this.state.isExpanded && <ActionMenu items={this.props.items} current={this.props.current} onSelect={this.onSelect} />}
    //         </div>
    //     }
    // }
}

type SectionProps = {
    items: ActionMenu.Items,
    onSelect: ActionMenu.OnSelect | ActionMenu.OnSelectMany,
    current: ActionMenu.Item | undefined,
    multiselect: boolean | undefined,
    noOffset?: boolean,
    noAccent?: boolean
}
type SectionState = { isExpanded: boolean, hasCurrent: boolean, header?: ActionMenu.Header }

class Section extends React.PureComponent<SectionProps, SectionState> {
    static createState(props: SectionProps, isExpanded?: boolean): SectionState {
        const header = isItems(props.items) && isHeader(props.items[0]) ? props.items[0] : void 0;

        const hasCurrent = header?.isIndependent
            ? false
            : props.multiselect
                ? ActionMenu.hasSelectedItem(props.items)
                : (!!props.current && !!ActionMenu.findItem(props.items, props.current.value)) || ActionMenu.hasSelectedItem(props.items);

        return {
            header,
            hasCurrent,
            isExpanded: hasCurrent || (isExpanded ?? !!header?.initiallyExpanded)
        };
    }

    state = Section.createState(this.props)

    toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isExpanded: !this.state.isExpanded });
        e.currentTarget.blur();
    }

    componentDidUpdate(prevProps: SectionProps) {
        if (this.props.items !== prevProps.items || this.props.current !== prevProps.current) {
            // keep previously expanded section if the header label is the same
            const isExpanded = (
                isItems(this.props.items) && isItems(prevProps.items) &&
                isHeader(this.props.items[0]) && isHeader(prevProps.items[0]) &&
                this.props.items[0].label === prevProps.items[0].label
            ) ? this.state.isExpanded : undefined;
            this.setState(Section.createState(this.props, isExpanded));
        }
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

        return <div className='msp-flex-row msp-control-group-header'>
            <Button icon={this.state.isExpanded ? ArrowDropDownSvg : ArrowRightSvg} flex noOverflow onClick={this.toggleExpanded} title={`Click to ${this.state.isExpanded ? 'collapse' : 'expand'}.${header?.description ? ` ${header?.description}` : ''}`}>
                {hasCurrent ? <b>{header?.label}</b> : header?.label}
            </Button>
            <Button icon={CheckSvg} flex onClick={this.selectAll} style={{ flex: '0 0 50px', textAlign: 'right' }}>
                All
            </Button>
            <Button icon={CloseSvg} flex onClick={this.selectNone} style={{ flex: '0 0 50px', textAlign: 'right' }}>
                None
            </Button>
        </div>;
    }

    get basicHeader() {
        const { header, hasCurrent } = this.state;

        return <div className='msp-control-group-header' style={{ marginTop: '1px' }}>
            <Button noOverflow icon={this.state.isExpanded ? ArrowDropDownSvg : ArrowRightSvg} onClick={this.toggleExpanded} title={`Click to ${this.state.isExpanded ? 'collapse' : 'expand'}. ${header?.description ? header?.description : ''}`}>
                {hasCurrent ? <b>{header?.label}</b> : header?.label}
            </Button>
        </div>;
    }

    render() {
        const { items, onSelect, current } = this.props;

        if (isHeader(items)) return null;
        if (isItem(items)) return <Action item={items} onSelect={onSelect} current={current} multiselect={this.props.multiselect} />;

        const { header } = this.state;
        return <>
            {header && (this.props.multiselect && this.state.isExpanded ? this.multiselectHeader : this.basicHeader)}
            <div className={this.props.noOffset ? void 0 : this.props.noAccent ? 'msp-control-offset' : 'msp-accent-offset'}>
                {(!header || this.state.isExpanded) && items.map((x, i) => {
                    if (isHeader(x)) return null;
                    if (isItem(x)) return <Action key={i} item={x} onSelect={onSelect} current={current} multiselect={this.props.multiselect} />;
                    return <Section key={i} items={x} onSelect={onSelect} current={current} multiselect={this.props.multiselect} noAccent />;
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

    return <Button icon={item.icon} noOverflow className='msp-action-menu-button' onClick={e => onSelect(multiselect ? [item] : item as any, e)} disabled={item.disabled} style={style} title={item.description}>
        {isCurrent || item.selected ? <b>{item.label}</b> : item.label}
        {item.addOn}
    </Button>;
};

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