/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State, StateTransform, StateTransformer, StateAction, StateObject, StateObjectCell } from '../../mol-state';
import * as React from 'react';
import { PurePluginUIComponent } from '../base';
import { ParameterControls, ParamOnChange } from '../controls/parameters';
import { PluginContext } from '../../mol-plugin/context';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Subject } from 'rxjs';
import { Icon, RefreshSvg, CheckSvg, ArrowRightSvg, ArrowDropDownSvg, TuneSvg } from '../controls/icons';
import { ExpandGroup, ToggleButton, Button, IconButton } from '../controls/common';

export { StateTransformParameters, TransformControlBase };

class StateTransformParameters extends PurePluginUIComponent<StateTransformParameters.Props> {
    validate(params: any) {
        // TODO
        return void 0;
    }

    areInitial(params: any) {
        return PD.areEqual(this.props.info.params, params, this.props.info.initialValues);
    }

    onChange: ParamOnChange = ({ name, value }) => {
        const params = { ...this.props.params, [name]: value };
        this.props.events.onChange(params, this.areInitial(params), this.validate(params));
    };

    render() {
        return <ParameterControls params={this.props.info.params} values={this.props.params} onChange={this.onChange} onEnter={this.props.events.onEnter} isDisabled={this.props.isDisabled} />;
    }
}


namespace StateTransformParameters {
    export interface Props {
        info: {
            params: PD.Params,
            initialValues: any,
            isEmpty: boolean
        },
        events: {
            onChange: (params: any, areInitial: boolean, errors?: string[]) => void,
            onEnter: () => void,
        }
        params: any,
        isDisabled?: boolean,
        a?: StateObject,
        b?: StateObject,
        bCell?: StateObjectCell
    }

    export type Class = React.ComponentClass<Props>

    function areParamsEmpty(params: PD.Params) {
        const keys = Object.keys(params);
        for (const k of keys) {
            if (!params[k].isHidden) return false;
        }
        return true;
    }

    export function infoFromAction(plugin: PluginContext, state: State, action: StateAction, nodeRef: StateTransform.Ref): Props['info'] {
        const source = state.cells.get(nodeRef)!.obj!;
        const params = action.definition.params ? action.definition.params(source, plugin) : { };
        const initialValues = PD.getDefaultValues(params);
        return {
            initialValues,
            params,
            isEmpty: areParamsEmpty(params)
        };
    }

    export function infoFromTransform(plugin: PluginContext, state: State, transform: StateTransform): Props['info'] {
        const cell = state.cells.get(transform.ref)!;
        // const source: StateObjectCell | undefined = (cell.sourceRef && state.cells.get(cell.sourceRef)!) || void 0;
        // const create = transform.transformer.definition.params;
        // const params = create ? create((source && source.obj) as any, plugin) : { };
        const params = (cell.params && cell.params.definition) || { };
        const initialValues = (cell.params && cell.params.values) || { };
        return {
            initialValues,
            params,
            isEmpty: areParamsEmpty(params)
        };
    }
}

namespace TransformControlBase {
    export interface ComponentState {
        params: any,
        error?: string,
        busy: boolean,
        isInitial: boolean,
        simpleOnly?: boolean,
        isCollapsed?: boolean
    }

    export interface CommonProps {
        simpleApply?: { header: string, icon?: React.FC, title?: string },
        noMargin?: boolean,
        applyLabel?: string,
        onApply?: () => void,
        autoHideApply?: boolean,
        wrapInExpander?: boolean,
        expanderHeaderLeftMargin?: string
    }
}

abstract class TransformControlBase<P, S extends TransformControlBase.ComponentState> extends PurePluginUIComponent<P & TransformControlBase.CommonProps, S> {
    abstract applyAction(): Promise<void>;
    abstract getInfo(): StateTransformParameters.Props['info'];
    abstract getHeader(): StateTransformer.Definition['display'] | 'none';
    abstract canApply(): boolean;
    abstract getTransformerId(): string;
    abstract canAutoApply(newParams: any): boolean;
    abstract applyText(): string;
    abstract isUpdate(): boolean;
    abstract getSourceAndTarget(): { a?: StateObject, b?: StateObject, bCell?: StateObjectCell };
    abstract state: S;

    private busy: Subject<boolean> = new Subject();

    private onEnter = () => {
        if (this.state.error) return;
        this.apply();
    }

    private autoApplyHandle: number | undefined = void 0;
    private clearAutoApply() {
        if (this.autoApplyHandle !== void 0) {
            clearTimeout(this.autoApplyHandle);
            this.autoApplyHandle = void 0;
        }
    }

    events: StateTransformParameters.Props['events'] = {
        onEnter: this.onEnter,
        onChange: (params, isInitial, errors) => {
            this.clearAutoApply();
            this.setState({ params, isInitial, error: errors && errors[0] }, () => {
                if (!isInitial && !this.state.error && this.canAutoApply(params)) {
                    this.clearAutoApply();
                    this.autoApplyHandle = setTimeout(this.apply, 50) as any as number;
                }
            });
        }
    }

    apply = async () => {
        this.clearAutoApply();
        this.setState({ busy: true });
        try {
            await this.applyAction();
        } catch {
            // eat errors because they should be handled elsewhere
        } finally {
            this.props.onApply?.();
            this.busy.next(false);
        }
    }

    componentDidMount() {
        this.subscribe(this.plugin.behaviors.state.isBusy, b => {
            if (this.state.busy !== b) this.busy.next(b);
        });
        this.subscribe(this.busy, busy => {
            if (this.state.busy !== busy) this.setState({ busy });
        });
    }

    refresh = () => {
        this.setState({ params: this.getInfo().initialValues, isInitial: true, error: void 0 });
    }

    setDefault = () => {
        const info = this.getInfo();
        const params = PD.getDefaultValues(info.params);
        this.setState({ params, isInitial: PD.areEqual(info.params, params, info.initialValues), error: void 0 });
    }

    toggleExpanded = () => {
        this.setState({ isCollapsed: !this.state.isCollapsed });
    }

    renderApply() {
        const canApply = this.canApply();

        if (this.props.autoHideApply && (!canApply || this.canAutoApply(this.state.params))) return null;

        return <div className='msp-transform-apply-wrap'>
            <IconButton svg={RefreshSvg} className='msp-transform-default-params' onClick={this.setDefault} disabled={this.state.busy} title='Set default params' />
            <div className={`msp-transform-apply-wider`}>
                <Button icon={canApply ? CheckSvg : void 0} className={`msp-btn-commit msp-btn-commit-${canApply ? 'on' : 'off'}`} onClick={this.apply} disabled={!canApply}>
                    {this.props.applyLabel || this.applyText()}
                </Button>
            </div>
        </div>;
    }

    renderDefault() {
        const info = this.getInfo();
        const isEmpty = info.isEmpty && this.isUpdate();

        const display = this.getHeader();

        const tId = this.getTransformerId();
        const ParamEditor: StateTransformParameters.Class = this.plugin.customParamEditors.has(tId)
            ? this.plugin.customParamEditors.get(tId)!
            : StateTransformParameters;

        const wrapClass = this.state.isCollapsed
            ? 'msp-transform-wrapper msp-transform-wrapper-collapsed'
            : 'msp-transform-wrapper';

        let params = null;
        if (!isEmpty && !this.state.isCollapsed) {
            const { a, b, bCell } = this.getSourceAndTarget();
            const applyControl = this.renderApply();
            params = <>
                <ParamEditor info={info} a={a} b={b} bCell={bCell} events={this.events} params={this.state.params} isDisabled={this.state.busy} />
                {applyControl}
            </>;
        }

        const ctrl = <div className={wrapClass} style={{ marginBottom: this.props.noMargin ? 0 : void 0 }}>
            {display !== 'none' && !this.props.wrapInExpander && <div className='msp-transform-header'>
                <Button onClick={this.toggleExpanded} title={display.description}>
                    {!isEmpty && <Icon svg={this.state.isCollapsed ? ArrowRightSvg : ArrowDropDownSvg} />}
                    {display.name}
                </Button>
            </div>}
            {params}
        </div>;

        if (isEmpty || !this.props.wrapInExpander) return ctrl;

        return <ExpandGroup header={this.isUpdate() ? `Update ${display === 'none' ? '' : display.name}` : `Apply ${display === 'none' ? '' : display.name}` } headerLeftMargin={this.props.expanderHeaderLeftMargin}>
            {ctrl}
        </ExpandGroup>;
    }

    renderSimple() {
        const info = this.getInfo();
        const canApply = this.canApply();

        const apply = <div className='msp-flex-row'>
            <Button icon={this.props.simpleApply?.icon} title={this.props.simpleApply?.title} disabled={this.state.busy || !canApply} onClick={this.apply} className='msp-btn-apply-simple'>
                {this.props.simpleApply?.header}
            </Button>
            {!info.isEmpty && <ToggleButton icon={TuneSvg} label='' title='Options' toggle={this.toggleExpanded} isSelected={!this.state.isCollapsed} disabled={this.state.busy} style={{ flex: '0 0 40px', padding: 0 }} />}
        </div>;

        if (this.state.isCollapsed) return apply;

        const tId = this.getTransformerId();
        const ParamEditor: StateTransformParameters.Class = this.plugin.customParamEditors.has(tId)
            ? this.plugin.customParamEditors.get(tId)!
            : StateTransformParameters;
        const { a, b, bCell } = this.getSourceAndTarget();

        return <>
            {apply}
            <ParamEditor info={info} a={a} b={b} bCell={bCell} events={this.events} params={this.state.params} isDisabled={this.state.busy} />
        </>;
    }

    render() {
        // console.log('rendering', ((this.props as any)?.transform?.transformer || (this.props as any)?.action)?.definition.display.name, +new Date)
        return this.props.simpleApply ? this.renderSimple() : this.renderDefault();
    }
}