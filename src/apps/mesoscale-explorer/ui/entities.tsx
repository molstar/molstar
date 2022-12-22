/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { PluginUIComponent } from '../../../mol-plugin-ui/base';
import { Button, IconButton } from '../../../mol-plugin-ui/controls/common';
import { CloseSvg, VisibilityOffOutlinedSvg, VisibilityOutlinedSvg } from '../../../mol-plugin-ui/controls/icons';
import { PluginCommands } from '../../../mol-plugin/commands';
import { State, StateObject, StateObjectCell, StateSelection, StateTransform, StateTransformer } from '../../../mol-state';
import { debounceTime, filter } from 'rxjs/operators';
import { escapeRegExp } from '../../../mol-util/string';
import { StructureElement, StructureProperties } from '../../../mol-model/structure';

export class EntityControls extends PluginUIComponent<{}, { filter: string, isDisabled: boolean }> {
    state = {
        filter: '',
        isDisabled: false,
    };

    componentDidMount() {
        this.subscribe(this.plugin.state.events.object.created, e => {
            if (PluginStateObject.Molecule.Structure.is(e.obj)) this.forceUpdate();
        });

        this.subscribe(this.plugin.state.events.object.removed, e => {
            if (PluginStateObject.Molecule.Structure.is(e.obj)) this.forceUpdate();
        });

        this.subscribe(this.plugin.state.data.behaviors.isUpdating, v => {
            this.setState({ isDisabled: v });
        });
    }

    render() {
        const disabled = this.state.isDisabled;
        const reFilter = new RegExp(escapeRegExp(this.state.filter), 'gi');
        const entities = this.plugin.state.data.select(StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure).withTag('Entity')).filter(cell => {
            if (!cell.obj) return false;

            const s = cell.obj.data;
            const l = StructureElement.Location.create(s, s.units[0], s.units[0].elements[0]);
            const label = StructureProperties.entity.pdbx_description(l)[0] || '';
            return label.match(reFilter) !== null;
        });
        return <>
            <div className={`msp-flex-row msp-control-row`} style={{ margin: '5px' }}>
                <input type='text'
                    value={this.state.filter}
                    placeholder='Search'
                    onChange={e => this.setState({ filter: e.target.value.trim().replace(/\s+/gi, ' ') })}
                    disabled={disabled}
                />
                <IconButton svg={CloseSvg} toggleState={false} disabled={disabled} onClick={() => this.setState({ filter: '' })} />
            </div>
            {entities.map(s => {
                return <EntityNode cell={s} key={s.transform.ref} />;
            })}
        </>;
    }
}

type StructureCell = StateObjectCell<PluginStateObject.Molecule.Structure, StateTransform<StateTransformer<StateObject<any, StateObject.Type<any>>, StateObject<any, StateObject.Type<any>>, any>>>

export class EntityNode extends PluginUIComponent<{ cell: StructureCell }> {
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
    }

    toggleVisible = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.State.ToggleVisibility(this.plugin, { state: this.props.cell.parent!, ref: this.ref });
        e.currentTarget.blur();
    };

    highlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.Interactivity.Object.Highlight(this.plugin, { state: this.props.cell.parent!, ref: this.ref });
        e.currentTarget.blur();
    };

    clearHighlight = (e: React.MouseEvent<HTMLElement>) => {
        e.preventDefault();
        PluginCommands.Interactivity.ClearHighlights(this.plugin);
        e.currentTarget.blur();
    };

    render() {
        const cell = this.props.cell;
        const obj = cell.obj!;

        console.log(obj.data);

        const cellState = cell.state;
        const disabled = cell.status !== 'error' && cell.status !== 'ok';

        const s = obj.data;
        const l = StructureElement.Location.create(s, s.units[0], s.units[0].elements[0]);

        const label = <Button className={`msp-btn-tree-label msp-type-class-${obj.type.typeClass}`} noOverflow disabled={disabled}>
            <span>{StructureProperties.entity.pdbx_description(l)[0]}</span>
        </Button>;

        const visibility = <IconButton svg={cellState.isHidden ? VisibilityOffOutlinedSvg : VisibilityOutlinedSvg} toggleState={false} disabled={disabled} small onClick={this.toggleVisible} />;

        return <div className={`msp-flex-row`} style={{ margin: '5px' }} onMouseEnter={this.highlight} onMouseLeave={this.clearHighlight}>
            {label}
            {visibility}
        </div>;
    }
}


