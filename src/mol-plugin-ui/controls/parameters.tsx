/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { Mat4, Vec2, Vec3 } from '../../mol-math/linear-algebra';
import { PluginContext } from '../../mol-plugin/context';
import { Color } from '../../mol-util/color';
import { ColorListName, ColorListOptions, ColorListOptionsScale, ColorListOptionsSet, getColorListFromName } from '../../mol-util/color/lists';
import { Legend as LegendData } from '../../mol-util/legend';
import { memoize1, memoizeLatest } from '../../mol-util/memoize';
import { getPrecision } from '../../mol-util/number';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ParamMapping } from '../../mol-util/param-mapping';
import { camelCaseToWords } from '../../mol-util/string';
import { PluginUIComponent } from '../base';
import { ActionMenu } from './action-menu';
import { ColorOptions, ColorValueOption, CombinedColorControl } from './color';
import { Button, ControlGroup, ControlRow, ExpandGroup, IconButton, TextInput, ToggleButton } from './common';
import { Icon, HelpOutlineSvg, CheckSvg, ClearSvg, BookmarksOutlinedSvg, MoreHorizSvg, ArrowDropDownSvg, ArrowRightSvg, ArrowDownwardSvg, ArrowUpwardSvg, DeleteOutlinedSvg } from './icons';
import { legendFor } from './legend';
import LineGraphComponent from './line-graph/line-graph-component';
import { Slider, Slider2 } from './slider';
import { Asset } from '../../mol-util/assets';

export type ParameterControlsCategoryFilter = string | null | (string | null)[]

export interface ParameterControlsProps<P extends PD.Params = PD.Params> {
    params: P,
    values: any,
    onChange?: ParamsOnChange<PD.ValuesFor<P>>,
    onChangeValues?: (values: PD.ValuesFor<P>, prev: PD.ValuesFor<P>) => void,
    isDisabled?: boolean,
    onEnter?: () => void
}

export class ParameterControls<P extends PD.Params> extends React.PureComponent<ParameterControlsProps<P>> {
    onChange: ParamOnChange = (params) => {
        this.props.onChange?.(params, this.props.values);
        if (this.props.onChangeValues) {
            const values = { ...this.props.values, [params.name]: params.value };
            this.props.onChangeValues(values, this.props.values);
        }
    }

    renderGroup(group: ParamInfo[]) {
        if (group.length === 0) return null;

        const values = this.props.values;
        let ctrls: JSX.Element[] | null = null;
        let category: string | undefined = void 0;

        for (const [key, p, Control] of group) {
            if (p.hideIf?.(values)) continue;

            if (!ctrls) ctrls = [];
            category = p.category;
            ctrls.push(<Control param={p} key={key} onChange={this.onChange} onEnter={this.props.onEnter} isDisabled={this.props.isDisabled} name={key} value={values[key]} />);
        }

        if (!ctrls) return null;

        if (category) {
            return [<ExpandGroup key={category} header={category}>{ctrls}</ExpandGroup>];
        }
        return ctrls;
    }

    renderPart(groups: ParamInfo[][]) {
        let parts: JSX.Element[] | null = null;
        for (const g of groups) {
            const ctrls = this.renderGroup(g);
            if (!ctrls) continue;
            if (!parts) parts = [];
            for (const c of ctrls) parts.push(c);
        }
        return parts;
    }

    paramGroups = memoizeLatest((params: PD.Params) => classifyParams(params));

    render() {
        const groups = this.paramGroups(this.props.params);

        const essentials = this.renderPart(groups.essentials);
        const advanced = this.renderPart(groups.advanced);

        if (essentials && advanced) {
            return <>
                {essentials}
                <ExpandGroup header='Advanced Options'>
                    {advanced}
                </ExpandGroup>
            </>;
        } else if (essentials) {
            return essentials;
        } else {
            return advanced;
        }
    }
}

export class ParameterMappingControl<S, T> extends PluginUIComponent<{ mapping: ParamMapping<S, T, PluginContext> }> {
    setSettings = (p: { param: PD.Base<any>, name: string, value: any }, old: any) => {
        const values = { ...old, [p.name]: p.value };
        const t = this.props.mapping.update(values, this.plugin);
        this.props.mapping.apply(t, this.plugin);
    }

    componentDidMount() {
        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => this.forceUpdate());
    }

    render() {
        const t = this.props.mapping.getTarget(this.plugin);
        const values = this.props.mapping.getValues(t, this.plugin);
        const params = this.props.mapping.params(this.plugin) as any as PD.Params;
        return <ParameterControls params={params} values={values} onChange={this.setSettings} />;
    }
}

type ParamInfo = [string, PD.Any, ParamControl];
function classifyParams(params: PD.Params) {
    function addParam(k: string, p: PD.Any, group: typeof essentials) {
        const ctrl = controlFor(p);
        if (!ctrl) return;

        if (!p.category) group.params[0].push([k, p, ctrl]);
        else {
            if (!group.map) group.map = new Map();
            let c = group.map.get(p.category);
            if (!c) {
                c = [];
                group.map.set(p.category, c);
                group.params.push(c);
            }
            c.push([k, p, ctrl]);
        }
    }

    function sortGroups(x: ParamInfo[], y: ParamInfo[]) {
        const a = x[0], b = y[0];
        if (!a || !a[1].category) return -1;
        if (!b || !b[1].category) return 1;
        return a[1].category < b[1].category ? -1 : 1;
    }

    const keys = Object.keys(params);

    const essentials: { params: ParamInfo[][], map: Map<string, ParamInfo[]> | undefined } = { params: [[]], map: void 0 };
    const advanced: typeof essentials = { params: [[]], map: void 0 };

    for (const k of keys) {
        const p = params[k];
        if (p.isHidden) continue;

        if (p.isEssential) addParam(k, p, essentials);
        else addParam(k, p, advanced);
    }

    essentials.params.sort(sortGroups);
    advanced.params.sort(sortGroups);

    return { essentials: essentials.params, advanced: advanced.params };
}

function controlFor(param: PD.Any): ParamControl | undefined {
    switch (param.type) {
        case 'value': return void 0;
        case 'boolean': return BoolControl;
        case 'number': return typeof param.min !== 'undefined' && typeof param.max !== 'undefined'
            ? NumberRangeControl : NumberInputControl;
        case 'converted': return ConvertedControl;
        case 'conditioned': return ConditionedControl;
        case 'multi-select': return MultiSelectControl;
        case 'color': return CombinedColorControl;
        case 'color-list': return ColorListControl;
        case 'vec3': return Vec3Control;
        case 'mat4': return Mat4Control;
        case 'url': return UrlControl;
        case 'file': return FileControl;
        case 'file-list': return FileListControl;
        case 'select': return SelectControl;
        case 'text': return TextControl;
        case 'interval': return typeof param.min !== 'undefined' && typeof param.max !== 'undefined'
            ? BoundedIntervalControl : IntervalControl;
        case 'group': return GroupControl;
        case 'mapped': return MappedControl;
        case 'line-graph': return LineGraphControl;
        case 'script': return ScriptControl;
        case 'object-list': return ObjectListControl;
        default:
            const _: never = param;
            console.warn(`${_} has no associated UI component`);
            return void 0;
    }
}

export class ParamHelp<L extends LegendData> extends React.PureComponent<{ legend?: L, description?: string }> {
    render() {
        const { legend, description } = this.props;
        const Legend = legend && legendFor(legend);

        return <div className='msp-help-text'>
            <div>
                <div className='msp-help-description'><Icon svg={HelpOutlineSvg} inline />{description}</div>
                {Legend && <div className='msp-help-legend'><Legend legend={legend} /></div>}
            </div>
        </div>;
    }
}

export type ParamsOnChange<P> = (params: { param: PD.Base<any>, name: string, value: any }, values: Readonly<P>) => void
export type ParamOnChange = (params: { param: PD.Base<any>, name: string, value: any }) => void
export interface ParamProps<P extends PD.Base<any> = PD.Base<any>> {
    name: string,
    value: P['defaultValue'],
    param: P,
    isDisabled?: boolean,
    onChange: ParamOnChange,
    onEnter?: () => void
}
export type ParamControl = React.ComponentClass<ParamProps<any>>

function renderSimple(options: { props: ParamProps<any>, state: { showHelp: boolean }, control: JSX.Element, addOn: JSX.Element | null, toggleHelp: () => void }) {
    const { props, state, control, toggleHelp, addOn } = options;

    let _className = [];
    if (props.param.shortLabel) _className.push('msp-control-label-short');
    if (props.param.twoColumns) _className.push('msp-control-col-2');
    const className = _className.join(' ');

    const label = props.param.label || camelCaseToWords(props.name);
    const help = props.param.help
        ? props.param.help(props.value)
        : { description: props.param.description, legend: props.param.legend };
    const desc = props.param.description;
    const hasHelp = help.description || help.legend;
    return <>
        <ControlRow
            className={className}
            title={desc}
            label={<>
                {label}
                {hasHelp &&
                    <button className='msp-help msp-btn-link msp-btn-icon msp-control-group-expander' onClick={toggleHelp}
                        title={desc || `${state.showHelp ? 'Hide' : 'Show'} help`}
                        style={{ background: 'transparent', textAlign: 'left', padding: '0' }}>
                        <Icon svg={HelpOutlineSvg} />
                    </button>
                }
            </>}
            control={control}
        />
        {hasHelp && state.showHelp && <div className='msp-control-offset'>
            <ParamHelp legend={help.legend} description={help.description} />
        </div>}
        {addOn}
    </>;
}

export abstract class SimpleParam<P extends PD.Any> extends React.PureComponent<ParamProps<P>, { showHelp: boolean }> {
    state = { showHelp: false };

    protected update(value: P['defaultValue']) {
        this.props.onChange({ param: this.props.param, name: this.props.name, value });
    }

    abstract renderControl(): JSX.Element;
    renderAddOn(): JSX.Element | null { return null; }

    toggleHelp = () => this.setState({ showHelp: !this.state.showHelp });

    render() {
        return renderSimple({
            props: this.props,
            state: this.state,
            control: this.renderControl(),
            toggleHelp: this.toggleHelp,
            addOn: this.renderAddOn()
        });
    }
}

export class BoolControl extends SimpleParam<PD.BooleanParam> {
    onClick = (e: React.MouseEvent<HTMLButtonElement>) => { this.update(!this.props.value); e.currentTarget.blur(); }
    renderControl() {
        return <button onClick={this.onClick} disabled={this.props.isDisabled}>
            <Icon svg={this.props.value ? CheckSvg : ClearSvg} />
            {this.props.value ? 'On' : 'Off'}
        </button>;
    }
}

export class LineGraphControl extends React.PureComponent<ParamProps<PD.LineGraph>, { isExpanded: boolean, isOverPoint: boolean, message: string }> {
    state = {
        isExpanded: false,
        isOverPoint: false,
        message: `${this.props.param.defaultValue.length} points`,
    }

    onHover = (point?: Vec2) => {
        this.setState({ isOverPoint: !this.state.isOverPoint });
        if (point) {
            this.setState({ message: `(${point[0].toFixed(2)}, ${point[1].toFixed(2)})` });
            return;
        }
        this.setState({ message: `${this.props.value.length} points` });
    }

    onDrag = (point: Vec2) => {
        this.setState({ message: `(${point[0].toFixed(2)}, ${point[1].toFixed(2)})` });
    }

    onChange = (value: PD.LineGraph['defaultValue']) => {
        this.props.onChange({ name: this.props.name, param: this.props.param, value: value });
    }

    toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isExpanded: !this.state.isExpanded });
        e.currentTarget.blur();
    }

    render() {
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        return <>
            <ControlRow label={label} control={<button onClick={this.toggleExpanded} disabled={this.props.isDisabled}>{`${this.state.message}`}</button>} />
            <div className='msp-control-offset' style={{ display: this.state.isExpanded ? 'block' : 'none' }}>
                <LineGraphComponent
                    data={this.props.param.defaultValue}
                    onChange={this.onChange}
                    onHover={this.onHover}
                    onDrag={this.onDrag} />
            </div>
        </>;
    }
}

export class NumberInputControl extends React.PureComponent<ParamProps<PD.Numeric>> {
    state = { value: '0' };

    update = (value: number) => {
        const p = getPrecision(this.props.param.step || 0.01);
        value = parseFloat(value.toFixed(p));
        this.props.onChange({ param: this.props.param, name: this.props.name, value });
    }

    render() {
        const placeholder = this.props.param.label || camelCaseToWords(this.props.name);
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        const p = getPrecision(this.props.param.step || 0.01);
        return <ControlRow
            title={this.props.param.description}
            label={label}
            control={<TextInput numeric
                value={parseFloat(this.props.value.toFixed(p))} onEnter={this.props.onEnter} placeholder={placeholder}
                isDisabled={this.props.isDisabled} onChange={this.update} />} />;
    }
}

export class NumberRangeControl extends SimpleParam<PD.Numeric> {
    onChange = (v: number) => { this.update(v); }
    renderControl() {
        const value = typeof this.props.value === 'undefined' ? this.props.param.defaultValue : this.props.value;
        return <Slider value={value} min={this.props.param.min!} max={this.props.param.max!}
            step={this.props.param.step} onChange={this.onChange} disabled={this.props.isDisabled} onEnter={this.props.onEnter} />;
    }
}

export class TextControl extends SimpleParam<PD.Text> {
    onChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        const value = e.target.value;
        if (value !== this.props.value) {
            this.update(value);
        }
    }

    onKeyPress = (e: React.KeyboardEvent<HTMLInputElement>) => {
        if ((e.keyCode === 13 || e.charCode === 13)) {
            if (this.props.onEnter) this.props.onEnter();
        }
        e.stopPropagation();
    }

    renderControl() {
        const placeholder = this.props.param.label || camelCaseToWords(this.props.name);
        return <input type='text'
            value={this.props.value || ''}
            placeholder={placeholder}
            onChange={this.onChange}
            onKeyPress={this.props.onEnter ? this.onKeyPress : void 0}
            disabled={this.props.isDisabled}
        />;
    }
}

export class PureSelectControl extends React.PureComponent<ParamProps<PD.Select<string | number>> & { title?: string }> {
    protected update(value: string | number) {
        this.props.onChange({ param: this.props.param, name: this.props.name, value });
    }

    onChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
        if (typeof this.props.param.defaultValue === 'number') {
            this.update(parseInt(e.target.value, 10));
        } else {
            this.update(e.target.value);
        }
    }

    render() {
        const isInvalid = this.props.value !== void 0 && !this.props.param.options.some(e => e[0] === this.props.value);
        return <select className='msp-form-control' title={this.props.title} value={this.props.value !== void 0 ? this.props.value : this.props.param.defaultValue} onChange={this.onChange} disabled={this.props.isDisabled}>
            {isInvalid && <option key={this.props.value} value={this.props.value}>{`[Invalid] ${this.props.value}`}</option>}
            {this.props.param.options.map(([value, label]) => <option key={value} value={value}>{label}</option>)}
        </select>;
    }
}

export class SelectControl extends React.PureComponent<ParamProps<PD.Select<string | number>>, { showHelp: boolean, showOptions: boolean }> {
    state = { showHelp: false, showOptions: false };

    onSelect: ActionMenu.OnSelect = item => {
        if (!item || item.value === this.props.value) {
            this.setState({ showOptions: false });
        } else {
            this.setState({ showOptions: false }, () => {
                this.props.onChange({ param: this.props.param, name: this.props.name, value: item.value });
            });
        }
    }

    toggle = () => this.setState({ showOptions: !this.state.showOptions });

    cycle = () => {
        const { options } = this.props.param;
        const current = options.findIndex(o => o[0] === this.props.value);
        const next = current === options.length - 1 ? 0 : current + 1;
        this.props.onChange({ param: this.props.param, name: this.props.name, value: options[next][0] });
    };

    items = memoizeLatest((param: PD.Select<any>) => ActionMenu.createItemsFromSelectOptions(param.options));

    renderControl() {
        const items = this.items(this.props.param);
        const current = this.props.value !== undefined ? ActionMenu.findItem(items, this.props.value) : void 0;
        const label = current
            ? current.label
            : typeof this.props.value === 'undefined'
                ? `${ActionMenu.getFirstItem(items)?.label || ''} [Default]`
                : `[Invalid] ${this.props.value}`;

        const toggle = this.props.param.cycle ? this.cycle : this.toggle;
        const textAlign = this.props.param.cycle ? 'center' : 'left';
        const icon = this.props.param.cycle
            ? (this.props.value === 'on' ? CheckSvg
                : this.props.value === 'off' ? ClearSvg : void 0)
            : void 0;

        return <ToggleButton disabled={this.props.isDisabled} style={{ textAlign, overflow: 'hidden', textOverflow: 'ellipsis' }}
            label={label} title={label as string} icon={icon} toggle={toggle} isSelected={this.state.showOptions} />;
    }

    renderAddOn() {
        if (!this.state.showOptions) return null;

        const items = this.items(this.props.param);
        const current = ActionMenu.findItem(items, this.props.value);

        return <ActionMenu items={items} current={current} onSelect={this.onSelect} />;
    }

    toggleHelp = () => this.setState({ showHelp: !this.state.showHelp });

    render() {
        return renderSimple({
            props: this.props,
            state: this.state,
            control: this.renderControl(),
            toggleHelp: this.toggleHelp,
            addOn: this.renderAddOn()
        });
    }
}

export class IntervalControl extends React.PureComponent<ParamProps<PD.Interval>, { isExpanded: boolean }> {
    state = { isExpanded: false }

    components = {
        0: PD.Numeric(0, { step: this.props.param.step }, { label: 'Min' }),
        1: PD.Numeric(0, { step: this.props.param.step }, { label: 'Max' })
    }

    change(value: PD.MultiSelect<any>['defaultValue']) {
        this.props.onChange({ name: this.props.name, param: this.props.param, value });
    }

    componentChange: ParamOnChange = ({ name, value }) => {
        const v = [...this.props.value];
        v[+name] = value;
        this.change(v);
    }

    toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isExpanded: !this.state.isExpanded });
        e.currentTarget.blur();
    }

    render() {
        const v = this.props.value;
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        const p = getPrecision(this.props.param.step || 0.01);
        const value = `[${v[0].toFixed(p)}, ${v[1].toFixed(p)}]`;
        return <>
            <ControlRow label={label} control={<button onClick={this.toggleExpanded} disabled={this.props.isDisabled}>{value}</button>} />
            {this.state.isExpanded && <div className='msp-control-offset'>
                <ParameterControls params={this.components} values={v} onChange={this.componentChange} onEnter={this.props.onEnter} />
            </div>}
        </>;
    }
}

export class BoundedIntervalControl extends SimpleParam<PD.Interval> {
    onChange = (v: [number, number]) => { this.update(v); }
    renderControl() {
        return <Slider2 value={this.props.value} min={this.props.param.min!} max={this.props.param.max!}
            step={this.props.param.step} onChange={this.onChange} disabled={this.props.isDisabled} onEnter={this.props.onEnter} />;
    }
}

export class ColorControl extends SimpleParam<PD.Color> {
    onChange = (e: React.ChangeEvent<HTMLSelectElement>) => {
        this.update(Color(parseInt(e.target.value)));
    }

    stripStyle(): React.CSSProperties {
        return {
            background: Color.toStyle(this.props.value),
            position: 'absolute',
            bottom: '0',
            height: '4px',
            right: '0',
            left: '0'
        };
    }

    renderControl() {
        return <div style={{ position: 'relative' }}>
            <select value={this.props.value} onChange={this.onChange}>
                {ColorValueOption(this.props.value)}
                {ColorOptions()}
            </select>
            <div style={this.stripStyle()} />
        </div>;
    }
}

const colorGradientInterpolated = memoize1((colors: Color[]) => {
    const styles = colors.map(c => Color.toStyle(c));
    return `linear-gradient(to right, ${styles.join(', ')})`;
});

const colorGradientBanded = memoize1((colors: Color[]) => {
    const n = colors.length;
    const styles: string[] = [`${Color.toStyle(colors[0])} ${100 * (1 / n)}%`];
    for (let i = 1, il = n - 1; i < il; ++i) {
        styles.push(
            `${Color.toStyle(colors[i])} ${100 * (i / n)}%`,
            `${Color.toStyle(colors[i])} ${100 * ((i + 1) / n)}%`
        );
    }
    styles.push(`${Color.toStyle(colors[n - 1])} ${100 * ((n - 1) / n)}%`);
    return `linear-gradient(to right, ${styles.join(', ')})`;
});

function colorStripStyle(list: PD.ColorList['defaultValue'], right = '0'): React.CSSProperties {
    return {
        background: colorGradient(list.colors, list.kind === 'set'),
        position: 'absolute',
        bottom: '0',
        height: '4px',
        right,
        left: '0'
    };
}

function colorGradient(colors: Color[], banded: boolean) {
    return banded ? colorGradientBanded(colors) : colorGradientInterpolated(colors);
}

function createColorListHelpers() {

    const addOn = (l: [ColorListName, any, any]) => {
        const preset = getColorListFromName(l[0]);
        return <div style={colorStripStyle({ kind: preset.type !== 'qualitative' ? 'interpolate' : 'set', colors: preset.list })} />;
    };

    return {
        ColorPresets: {
            all: ActionMenu.createItemsFromSelectOptions(ColorListOptions, { addOn }),
            scale: ActionMenu.createItemsFromSelectOptions(ColorListOptionsScale, { addOn }),
            set: ActionMenu.createItemsFromSelectOptions(ColorListOptionsSet, { addOn })
        },
        ColorsParam: PD.ObjectList({ color: PD.Color(0x0 as Color) }, ({ color }) => Color.toHexString(color).toUpperCase()),
        IsInterpolatedParam: PD.Boolean(false, { label: 'Interpolated' })
    };
}

let _colorListHelpers: ReturnType<typeof createColorListHelpers>;
function ColorListHelpers() {
    if (_colorListHelpers) return _colorListHelpers;
    _colorListHelpers = createColorListHelpers();
    return _colorListHelpers;
}

export class ColorListControl extends React.PureComponent<ParamProps<PD.ColorList>, { showHelp: boolean, show?: 'edit' | 'presets' }> {
    state = { showHelp: false, show: void 0 as 'edit' | 'presets' | undefined };

    protected update(value: PD.ColorList['defaultValue']) {
        this.props.onChange({ param: this.props.param, name: this.props.name, value });
    }

    toggleEdit = () => this.setState({ show: this.state.show === 'edit' ? void 0 : 'edit' });
    togglePresets = () => this.setState({ show: this.state.show === 'presets' ? void 0 : 'presets' });

    renderControl() {
        const { value } = this.props;
        // TODO: fix the button right offset
        return <>
            <button onClick={this.toggleEdit} style={{ position: 'relative', paddingRight: '33px' }}>
                {value.colors.length === 1 ? '1 color' : `${value.colors.length} colors`}
                <div style={colorStripStyle(value, '33px')} />
            </button>
            <IconButton svg={BookmarksOutlinedSvg} onClick={this.togglePresets} toggleState={this.state.show === 'presets'} title='Color Presets'
                style={{ padding: 0, position: 'absolute', right: 0, top: 0, width: '32px' }} />
        </>;
    }

    selectPreset: ActionMenu.OnSelect = item => {
        if (!item) return;
        this.setState({ show: void 0 });

        const preset = getColorListFromName(item.value as ColorListName);
        this.update({ kind: preset.type !== 'qualitative' ? 'interpolate' : 'set', colors: preset.list });
    }

    colorsChanged: ParamOnChange = ({ value }) => {
        this.update({
            kind: this.props.value.kind,
            colors: (value as (typeof _colorListHelpers)['ColorsParam']['defaultValue']).map(c => c.color)
        });
    }

    isInterpolatedChanged: ParamOnChange = ({ value }) => {
        this.update({ kind: value ? 'interpolate' : 'set', colors: this.props.value.colors });
    }

    renderColors() {
        if (!this.state.show) return null;
        const { ColorPresets, ColorsParam, IsInterpolatedParam } = ColorListHelpers();

        const preset = ColorPresets[this.props.param.presetKind];
        if (this.state.show === 'presets') return <ActionMenu items={preset} onSelect={this.selectPreset} />;

        const values = this.props.value.colors.map(color => ({ color }));
        return <div className='msp-control-offset'>
            <ObjectListControl name='colors' param={ColorsParam} value={values} onChange={this.colorsChanged} isDisabled={this.props.isDisabled} onEnter={this.props.onEnter} />
            <BoolControl name='isInterpolated' param={IsInterpolatedParam} value={this.props.value.kind === 'interpolate'} onChange={this.isInterpolatedChanged} isDisabled={this.props.isDisabled} onEnter={this.props.onEnter} />
        </div>;
    }

    toggleHelp = () => this.setState({ showHelp: !this.state.showHelp });

    render() {
        return renderSimple({
            props: this.props,
            state: this.state,
            control: this.renderControl(),
            toggleHelp: this.toggleHelp,
            addOn: this.renderColors()
        });
    }
}

export class Vec3Control extends React.PureComponent<ParamProps<PD.Vec3>, { isExpanded: boolean }> {
    state = { isExpanded: false }

    components = {
        0: PD.Numeric(0, { step: this.props.param.step }, { label: (this.props.param.fieldLabels && this.props.param.fieldLabels.x) || 'X' }),
        1: PD.Numeric(0, { step: this.props.param.step }, { label: (this.props.param.fieldLabels && this.props.param.fieldLabels.y) || 'Y' }),
        2: PD.Numeric(0, { step: this.props.param.step }, { label: (this.props.param.fieldLabels && this.props.param.fieldLabels.z) || 'Z' })
    }

    change(value: PD.MultiSelect<any>['defaultValue']) {
        this.props.onChange({ name: this.props.name, param: this.props.param, value });
    }

    componentChange: ParamOnChange = ({ name, value }) => {
        const v = Vec3.copy(Vec3.zero(), this.props.value);
        v[+name] = value;
        this.change(v);
    }

    toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isExpanded: !this.state.isExpanded });
        e.currentTarget.blur();
    }

    render() {
        const v = this.props.value;
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        const p = getPrecision(this.props.param.step || 0.01);
        const value = `[${v[0].toFixed(p)}, ${v[1].toFixed(p)}, ${v[2].toFixed(p)}]`;
        return <>
            <ControlRow label={label} control={<button onClick={this.toggleExpanded} disabled={this.props.isDisabled}>{value}</button>} />
            {this.state.isExpanded && <div className='msp-control-offset'>
                <ParameterControls params={this.components} values={v} onChange={this.componentChange} onEnter={this.props.onEnter} />
            </div>}
        </>;
    }
}

export class Mat4Control extends React.PureComponent<ParamProps<PD.Mat4>, { isExpanded: boolean }> {
    state = { isExpanded: false }

    components = {
        json: PD.Text(JSON.stringify(Mat4()), { description: 'JSON array with 4x4 matrix in a column major (j * 4 + i indexing) format' })
    }

    change(value: PD.MultiSelect<any>['defaultValue']) {
        this.props.onChange({ name: this.props.name, param: this.props.param, value });
    }

    componentChange: ParamOnChange = ({ name, value }) => {
        const v = Mat4.copy(Mat4(), this.props.value);
        if (name === 'json') {
            Mat4.copy(v, JSON.parse(value));
        } else {
            v[+name] = value;
        }

        this.change(v);
    }

    toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isExpanded: !this.state.isExpanded });
        e.currentTarget.blur();
    }

    changeValue(idx: number) {
        return (v: number) => {
            const m = Mat4.copy(Mat4(), this.props.value);
            m[idx] = v;
            this.change(m);
        };
    }

    get grid() {
        const v = this.props.value;
        const rows: React.ReactNode[] = [];
        for (let i = 0; i < 4; i++) {
            const row: React.ReactNode[] = [];
            for (let j = 0; j < 4; j++) {
                row.push(<TextInput key={j} numeric delayMs={50} value={Mat4.getValue(v, i, j)} onChange={this.changeValue(4 * j + i)} className='msp-form-control' blurOnEnter={true} isDisabled={this.props.isDisabled} />);
            }
            rows.push(<div className='msp-flex-row' key={i}>{row}</div>);
        }
        return <div className='msp-parameter-matrix'>{rows}</div>;
    }

    render() {
        const v = {
            json: JSON.stringify(this.props.value)
        };
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        return <>
            <ControlRow label={label} control={<button onClick={this.toggleExpanded} disabled={this.props.isDisabled}>{'4\u00D74 Matrix'}</button>} />
            {this.state.isExpanded && <div className='msp-control-offset'>
                {this.grid}
                <ParameterControls params={this.components} values={v} onChange={this.componentChange} onEnter={this.props.onEnter} />
            </div>}
        </>;
    }
}

export class UrlControl extends SimpleParam<PD.UrlParam> {
    onChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        const value = e.target.value;
        if (value !== Asset.getUrl(this.props.value || '')) {
            this.update(Asset.Url(value));
        }
    }

    onKeyPress = (e: React.KeyboardEvent<HTMLInputElement>) => {
        if ((e.keyCode === 13 || e.charCode === 13)) {
            if (this.props.onEnter) this.props.onEnter();
        }
        e.stopPropagation();
    }

    renderControl() {
        const placeholder = this.props.param.label || camelCaseToWords(this.props.name);
        return <input type='text'
            value={Asset.getUrl(this.props.value || '')}
            placeholder={placeholder}
            onChange={this.onChange}
            onKeyPress={this.props.onEnter ? this.onKeyPress : void 0}
            disabled={this.props.isDisabled}
        />;
    }
}

export class FileControl extends React.PureComponent<ParamProps<PD.FileParam>> {
    change(value: File) {
        this.props.onChange({ name: this.props.name, param: this.props.param, value: Asset.File(value) });
    }

    onChangeFile = (e: React.ChangeEvent<HTMLInputElement>) => {
        this.change(e.target.files![0]);
    }

    render() {
        const value = this.props.value;

        return <div className='msp-btn msp-btn-block msp-btn-action msp-loader-msp-btn-file' style={{ marginTop: '1px' }}>
            {value ? value.name : 'Select a file...'} <input disabled={this.props.isDisabled} onChange={this.onChangeFile} type='file' multiple={false} accept={this.props.param.accept} />
        </div>;
    }
}

export class FileListControl extends React.PureComponent<ParamProps<PD.FileListParam>> {
    change(value: FileList) {
        const files: Asset.File[] = [];
        if (value) {
            for (let i = 0, il = value.length; i < il; ++i) {
                files.push(Asset.File(value[i]));
            }
        }
        this.props.onChange({ name: this.props.name, param: this.props.param, value: files });
    }

    onChangeFileList = (e: React.ChangeEvent<HTMLInputElement>) => {
        this.change(e.target.files!);
    }

    render() {
        const value = this.props.value;

        const names: string[] = [];
        if (value) {
            for (const file of value) {
                names.push(file.name);
            }
        }
        const label = names.length === 0
            ? 'Select files...' : names.length === 1
                ? names[0] : `${names.length} files selected`;

        return <div className='msp-btn msp-btn-block msp-btn-action msp-loader-msp-btn-file' style={{ marginTop: '1px' }}>
            {label} <input disabled={this.props.isDisabled} onChange={this.onChangeFileList} type='file' multiple={true} accept={this.props.param.accept} />
        </div>;
    }
}

export class MultiSelectControl extends React.PureComponent<ParamProps<PD.MultiSelect<any>>, { isExpanded: boolean }> {
    state = { isExpanded: false }

    change(value: PD.MultiSelect<any>['defaultValue']) {
        this.props.onChange({ name: this.props.name, param: this.props.param, value });
    }

    toggle(key: string) {
        return (e: React.MouseEvent<HTMLButtonElement>) => {
            if (this.props.value.indexOf(key) < 0) this.change(this.props.value.concat(key));
            else this.change(this.props.value.filter(v => v !== key));
            e.currentTarget.blur();
        };
    }

    toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isExpanded: !this.state.isExpanded });
        e.currentTarget.blur();
    }

    render() {
        const current = this.props.value;
        const emptyLabel = this.props.param.emptyValue;
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        return <>
            <ControlRow label={label} control={<button onClick={this.toggleExpanded} disabled={this.props.isDisabled}>
                {current.length === 0 && emptyLabel ? emptyLabel : `${current.length} of ${this.props.param.options.length}`}
            </button>} />
            {this.state.isExpanded && <div className='msp-control-offset'>
                {this.props.param.options.map(([value, label]) => {
                    const sel = current.indexOf(value) >= 0;
                    return <Button key={value} onClick={this.toggle(value)} disabled={this.props.isDisabled} style={{ marginTop: '1px' }}>
                        <span style={{ float: sel ? 'left' : 'right' }}>{sel ? `✓ ${label}` : `${label} ✗`}</span>
                    </Button>;
                })}
            </div>}
        </>;
    }
}

export class GroupControl extends React.PureComponent<ParamProps<PD.Group<any>> & { inMapped?: boolean }, { isExpanded: boolean, showHelp: boolean }> {
    state = { isExpanded: !!this.props.param.isExpanded, showHelp: false }

    change(value: any) {
        this.props.onChange({ name: this.props.name, param: this.props.param, value });
    }

    onChangeParam: ParamOnChange = e => {
        this.change({ ...this.props.value, [e.name]: e.value });
    }

    toggleExpanded = () => this.setState({ isExpanded: !this.state.isExpanded });

    pivoted() {
        const key = this.props.param.pivot as string;
        const params = this.props.param.params;
        const pivot = params[key];
        const Control = controlFor(pivot)!;

        const ctrl = <Control name={key} param={pivot} value={this.props.value[key]} onChange={this.onChangeParam} onEnter={this.props.onEnter} isDisabled={this.props.isDisabled} />;

        if (!this.state.isExpanded) {
            return <div className='msp-mapped-parameter-group'>
                {ctrl}
                <IconButton svg={MoreHorizSvg} onClick={this.toggleExpanded} toggleState={this.state.isExpanded} title={`More Options`} />
            </div>;
        }

        const filtered = Object.create(null);
        for (const k of Object.keys(params)) {
            if (k !== key) filtered[k] = params[k];
        }

        return <div className='msp-mapped-parameter-group'>
            {ctrl}
            <IconButton svg={MoreHorizSvg} onClick={this.toggleExpanded} toggleState={this.state.isExpanded} title={`More Options`} />
            <div className='msp-control-offset'>
                <ParameterControls params={filtered} onEnter={this.props.onEnter} values={this.props.value} onChange={this.onChangeParam} isDisabled={this.props.isDisabled} />
            </div>
        </div>;
    }

    render() {
        const params = this.props.param.params;

        // Do not show if there are no params.
        if (Object.keys(params).length === 0) return null;

        if (this.props.param.pivot) return this.pivoted();

        const label = this.props.param.label || camelCaseToWords(this.props.name);

        const controls = <ParameterControls params={params} onChange={this.onChangeParam} values={this.props.value} onEnter={this.props.onEnter} isDisabled={this.props.isDisabled} />;

        if (this.props.inMapped) {
            return <div className='msp-control-offset'>{controls}</div>;
        }

        if (this.props.param.isFlat) {
            return controls;
        }

        return <div className='msp-control-group-wrapper'>
            <div className='msp-control-group-header'>
                <button className='msp-btn msp-form-control msp-btn-block' onClick={this.toggleExpanded}>
                    <Icon svg={this.state.isExpanded ? ArrowDropDownSvg : ArrowRightSvg} />
                    {label}
                </button>
            </div>
            {this.state.isExpanded && <div className='msp-control-offset'>
                {controls}
            </div>}
        </div>;
    }
}

export class MappedControl extends React.PureComponent<ParamProps<PD.Mapped<any>>, { isExpanded: boolean }> {
    state = { isExpanded: false }

    // TODO: this could lead to a rare bug where the component is reused with different mapped control.
    // I think there are currently no cases where this could happen in the UI, but still need to watch out..
    private valuesCache: { [name: string]: PD.Values<any> } = {}
    private setValues(name: string, values: PD.Values<any>) {
        this.valuesCache[name] = values;
    }
    private getValues(name: string) {
        if (name in this.valuesCache) {
            return this.valuesCache[name];
        } else {
            return this.props.param.map(name).defaultValue;
        }
    }

    change(value: PD.Mapped<any>['defaultValue']) {
        this.props.onChange({ name: this.props.name, param: this.props.param, value });
    }

    onChangeName: ParamOnChange = e => {
        this.change({ name: e.value, params: this.getValues(e.value) });
    }

    onChangeParam: ParamOnChange = e => {
        this.setValues(this.props.value.name, e.value);
        this.change({ name: this.props.value.name, params: e.value });
    }

    toggleExpanded = () => this.setState({ isExpanded: !this.state.isExpanded });

    areParamsEmpty(params: PD.Params) {
        for (const k of Object.keys(params)) {
            if (!params[k].isHidden) return false;
        }
        return true;
    }

    render() {
        const value: PD.Mapped<any>['defaultValue'] = this.props.value;
        const param = this.props.param.map(value.name);
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        const Mapped = controlFor(param);

        const help = this.props.param.help;
        const select = help
            ? {
                ...this.props.param.select,
                help: (name: any) => help({ name, params: this.getValues(name) })
            }
            : this.props.param.select;

        const Select = <SelectControl param={select}
            isDisabled={this.props.isDisabled} onChange={this.onChangeName} onEnter={this.props.onEnter}
            name={label} value={value.name} />;

        if (!Mapped) {
            return Select;
        }

        if (param.type === 'group' && !param.isFlat) {
            if (!this.areParamsEmpty(param.params)) {
                return <div className='msp-mapped-parameter-group'>
                    {Select}
                    <IconButton svg={MoreHorizSvg} onClick={this.toggleExpanded} toggleState={this.state.isExpanded} title={`${label} Properties`} />
                    {this.state.isExpanded && <GroupControl inMapped param={param} value={value.params} name={value.name} onChange={this.onChangeParam} onEnter={this.props.onEnter} isDisabled={this.props.isDisabled} />}
                </div>;
            }

            return Select;
        }

        return <>
            {Select}
            <Mapped param={param} value={value.params} name={value.name} onChange={this.onChangeParam} onEnter={this.props.onEnter} isDisabled={this.props.isDisabled} />
        </>;
    }
}

type ObjectListEditorProps = { params: PD.Params, value: object, isUpdate?: boolean, apply: (value: any) => void, isDisabled?: boolean }
class ObjectListEditor extends React.PureComponent<ObjectListEditorProps, { current: object }> {
    state = { current: this.props.value };

    onChangeParam: ParamOnChange = e => {
        this.setState({ current: { ...this.state.current, [e.name]: e.value } });
    }

    apply = () => {
        this.props.apply(this.state.current);
    }

    componentDidUpdate(prevProps: ObjectListEditorProps) {
        if (this.props.params !== prevProps.params || this.props.value !== prevProps.value) {
            this.setState({ current: this.props.value });
        }
    }

    render() {
        return <>
            <ParameterControls params={this.props.params} onChange={this.onChangeParam} values={this.state.current} onEnter={this.apply} isDisabled={this.props.isDisabled} />
            <button className={`msp-btn msp-btn-block msp-form-control msp-control-top-offset`} onClick={this.apply} disabled={this.props.isDisabled}>
                {this.props.isUpdate ? 'Update' : 'Add'}
            </button>
        </>;
    }
}

type ObjectListItemProps = { param: PD.ObjectList, value: object, index: number, actions: ObjectListControl['actions'], isDisabled?: boolean }
class ObjectListItem extends React.PureComponent<ObjectListItemProps, { isExpanded: boolean }> {
    state = { isExpanded: false };

    update = (v: object) => {
        // this.setState({ isExpanded: false }); // TODO auto update? mark changed state?
        this.props.actions.update(v, this.props.index);
    }

    moveUp = () => {
        this.props.actions.move(this.props.index, -1);
    };

    moveDown = () => {
        this.props.actions.move(this.props.index, 1);
    };

    remove = () => {
        this.setState({ isExpanded: false });
        this.props.actions.remove(this.props.index);
    };

    toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isExpanded: !this.state.isExpanded });
        e.currentTarget.blur();
    };

    render() {
        return <>
            <div className='msp-param-object-list-item'>
                <button className='msp-btn msp-btn-block msp-form-control' onClick={this.toggleExpanded}>
                    <span>{`${this.props.index + 1}: `}</span>
                    {this.props.param.getLabel(this.props.value)}
                </button>
                <div>
                    <IconButton svg={ArrowDownwardSvg} title='Move Up' onClick={this.moveUp} small={true} />
                    <IconButton svg={ArrowUpwardSvg} title='Move Down' onClick={this.moveDown} small={true} />
                    <IconButton svg={DeleteOutlinedSvg} title='Remove' onClick={this.remove} small={true} />
                </div>
            </div>
            {this.state.isExpanded && <div className='msp-control-offset'>
                <ObjectListEditor params={this.props.param.element} apply={this.update} value={this.props.value} isUpdate isDisabled={this.props.isDisabled} />
            </div>}
        </>;
    }
}

export class ObjectListControl extends React.PureComponent<ParamProps<PD.ObjectList>, { isExpanded: boolean }> {
    state = { isExpanded: false }

    change(value: any) {
        this.props.onChange({ name: this.props.name, param: this.props.param, value });
    }

    add = (v: object) => {
        this.change([...this.props.value, v]);
    };

    actions = {
        update: (v: object, i: number) => {
            const value = this.props.value.slice(0);
            value[i] = v;
            this.change(value);
        },
        move: (i: number, dir: -1 | 1) => {
            let xs = this.props.value;
            if (xs.length === 1) return;

            let j = (i + dir) % xs.length;
            if (j < 0) j += xs.length;

            xs = xs.slice(0);
            const t = xs[i];
            xs[i] = xs[j];
            xs[j] = t;
            this.change(xs);
        },
        remove: (i: number) => {
            const xs = this.props.value;
            const update: object[] = [];
            for (let j = 0; j < xs.length; j++) {
                if (i !== j) update.push(xs[j]);
            }
            this.change(update);
        }
    }

    toggleExpanded = (e: React.MouseEvent<HTMLButtonElement>) => {
        this.setState({ isExpanded: !this.state.isExpanded });
        e.currentTarget.blur();
    };

    render() {
        const v = this.props.value;
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        const value = `${v.length} item${v.length !== 1 ? 's' : ''}`;
        return <>
            <ControlRow label={label} control={<button onClick={this.toggleExpanded} disabled={this.props.isDisabled}>{value}</button>} />
            {this.state.isExpanded && <div className='msp-control-offset'>
                {this.props.value.map((v, i) => <ObjectListItem key={i} param={this.props.param} value={v} index={i} actions={this.actions} isDisabled={this.props.isDisabled} />)}
                <ControlGroup header='New Item'>
                    <ObjectListEditor params={this.props.param.element} apply={this.add} value={this.props.param.ctor()} isDisabled={this.props.isDisabled} />
                </ControlGroup>
            </div>}
        </>;
    }
}

export class ConditionedControl extends React.PureComponent<ParamProps<PD.Conditioned<any, any, any>>> {
    change(value: PD.Conditioned<any, any, any>['defaultValue']) {
        this.props.onChange({ name: this.props.name, param: this.props.param, value });
    }

    onChangeCondition: ParamOnChange = e => {
        this.change(this.props.param.conditionedValue(this.props.value, e.value));
    }

    onChangeParam: ParamOnChange = e => {
        this.change(e.value);
    }

    render() {
        const value = this.props.value;
        const condition = this.props.param.conditionForValue(value) as string;
        const param = this.props.param.conditionParams[condition];
        const label = this.props.param.label || camelCaseToWords(this.props.name);
        const Conditioned = controlFor(param);

        const select = <SelectControl param={this.props.param.select}
            isDisabled={this.props.isDisabled} onChange={this.onChangeCondition} onEnter={this.props.onEnter}
            name={`${label} Kind`} value={condition} />;

        if (!Conditioned) {
            return select;
        }

        return <>
            {select}
            <Conditioned param={param} value={value} name={label} onChange={this.onChangeParam} onEnter={this.props.onEnter} isDisabled={this.props.isDisabled} />
        </>;
    }
}

export class ConvertedControl extends React.PureComponent<ParamProps<PD.Converted<any, any>>> {
    onChange: ParamOnChange = e => {
        this.props.onChange({
            name: this.props.name,
            param: this.props.param,
            value: this.props.param.toValue(e.value)
        });
    }

    render() {
        const value = this.props.param.fromValue(this.props.value);
        const Converted = controlFor(this.props.param.converted);

        if (!Converted) return null;
        return <Converted param={this.props.param.converted} value={value} name={this.props.name} onChange={this.onChange} onEnter={this.props.onEnter} isDisabled={this.props.isDisabled} />;
    }
}

export class ScriptControl extends SimpleParam<PD.Script> {
    onChange = (e: React.ChangeEvent<HTMLInputElement>) => {
        const value = e.target.value;
        if (value !== this.props.value.expression) {
            this.update({ language: this.props.value.language, expression: value });
        }
    }

    onKeyPress = (e: React.KeyboardEvent<HTMLInputElement>) => {
        if ((e.keyCode === 13 || e.charCode === 13)) {
            if (this.props.onEnter) this.props.onEnter();
        }
        e.stopPropagation();
    }

    renderControl() {
        // TODO: improve!

        const placeholder = this.props.param.label || camelCaseToWords(this.props.name);
        return <input type='text'
            value={this.props.value.expression || ''}
            placeholder={placeholder}
            onChange={this.onChange}
            onKeyPress={this.props.onEnter ? this.onKeyPress : void 0}
            disabled={this.props.isDisabled}
        />;
    }
}