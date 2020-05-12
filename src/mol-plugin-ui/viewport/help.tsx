/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { Binding } from '../../mol-util/binding';
import { PluginUIComponent } from '../base';
import { StateTransformer, StateSelection } from '../../mol-state';
import { SelectLoci } from '../../mol-plugin/behavior/dynamic/representation';
import { FocusLoci } from '../../mol-plugin/behavior/dynamic/representation';
import { Icon, ArrowDropDownSvg, ArrowRightSvg, CameraSvg } from '../controls/icons';
import { Button } from '../controls/common';

function getBindingsList(bindings: { [k: string]: Binding }) {
    return Object.keys(bindings).map(k => [k, bindings[k]] as [string, Binding]);
}

export class BindingsHelp extends React.PureComponent<{ bindings: { [k: string]: Binding } }> {
    getBindingComponents() {
        const bindingsList = getBindingsList(this.props.bindings);
        return <>
            {bindingsList.map(value => {
                const [name, binding] = value;
                return !Binding.isEmpty(binding)
                    ? <div key={name} style={{ marginBottom: '6px' }}>
                        <b>{binding.action}</b><br /><span dangerouslySetInnerHTML={{ __html: Binding.format(binding, name) }} />
                    </div>
                    : null;
            })}
        </>;
    }

    render() {
        return <HelpText>{this.getBindingComponents()}</HelpText>;
    }
}

export class HelpText extends React.PureComponent {
    render() {
        return <div className='msp-help-text'>
            <div>{this.props.children}</div>
        </div>;
    }
}

export class HelpGroup extends React.PureComponent<{ header: string, initiallyExpanded?: boolean }, { isExpanded: boolean }> {
    state = {
        header: this.props.header,
        isExpanded: !!this.props.initiallyExpanded
    }

    toggleExpanded = () => this.setState({ isExpanded: !this.state.isExpanded });

    render() {
        return <div className='msp-control-group-wrapper'>
            <div className='msp-control-group-header'>
                <Button onClick={this.toggleExpanded}>
                    <Icon svg={this.state.isExpanded ? ArrowDropDownSvg : ArrowRightSvg} />
                    {this.props.header}
                </Button>
            </div>
            {this.state.isExpanded && <div className='msp-control-offset' style={{ display: this.state.isExpanded ? 'block' : 'none' }}>
                {this.props.children}
            </div>}
        </div>;
    }
}

function HelpSection(props: { header: string }) {
    return <div className='msp-simple-help-section'>{props.header}</div>;
}

export class ViewportHelpContent extends PluginUIComponent<{ selectOnly?: boolean }> {
    componentDidMount() {
        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => this.forceUpdate());
    }

    render() {
        const interactionBindings: { [k: string]: Binding } = {};
        this.plugin.spec.behaviors.forEach(b => {
            const { bindings } = b.defaultParams;
            if (bindings) Object.assign(interactionBindings, bindings);
        });
        return <>
            {(!this.props.selectOnly && this.plugin.canvas3d) && <HelpGroup key='trackball' header='Moving in 3D'>
                <BindingsHelp bindings={this.plugin.canvas3d.props.trackball.bindings} />
            </HelpGroup>}
            <HelpGroup key='interactions' header='Mouse Controls'>
                <BindingsHelp bindings={interactionBindings} />
            </HelpGroup>
        </>;
    }
}

export class HelpContent extends PluginUIComponent {
    componentDidMount() {
        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => this.forceUpdate());
    }

    private formatTriggers(binding: Binding) {
        return binding.triggers.map(t => Binding.Trigger.format(t)).join(' or ');
    }

    private getTriggerFor(transformer: StateTransformer, name: string) {
        const state = this.plugin.state.behaviors;
        const selections = state.select(StateSelection.Generators.ofTransformer(transformer));
        const params = selections.length === 1 ? selections[0].params : undefined;
        const bindings = params ? params.values.bindings : {};
        const binding: Binding = name in bindings ? bindings[name] : Binding.Empty;
        return this.formatTriggers(binding);
    }

    render() {
        const selectToggleTriggers = this.getTriggerFor(SelectLoci, 'clickSelectToggle');
        const focusTriggers = this.getTriggerFor(FocusLoci, 'clickFocus');

        // TODO: interactive help, for example for density

        return <div>
            <HelpSection header='Interface Controls' />
            <HelpGroup header='Inline Help'>
                <HelpText>Many user interface elements show a little questionmark icon when hovered over. Clicking the icon toggles the display of an inline help text.</HelpText>
                <HelpText>Tooltips may provide additional information on a user interface element and are shown when hovering over it with the mouse.</HelpText>
            </HelpGroup>
            <HelpGroup header='Selections'>
                <HelpText>
                    The viewer allows changing colors and representations for selections of atoms, residues or chains. Selections can be created by
                    <ul style={{ paddingLeft: '20px' }}>
                        <li>picking elements on the 3D canvas or the sequence view using the mouse, e.g. toggle selection using {selectToggleTriggers} (for more see help section on <i>Mouse Controls</i>)</li>
                        <li>using the <i>Add</i>, <i>Remove</i> and <i>Only</i> dropdown buttons in the <i>Manage Selection</i> panel which allow modifing the current selection by predefined sets</li>
                    </ul>
                </HelpText>
            </HelpGroup>
            <HelpGroup header='Coloring'>
                <HelpText>
                    There are two ways to color structures. Every representation (e.g. cartoon or spacefill) has a color theme which can be changed using the dropdown for each representation in the <i>Structure Settings</i> panel. Additionally any selection atoms, residues or chains can by given a custom color. For that, first select the parts of the structure to be colored (see help section on <i>Selections</i>) and, second, choose a color from the color dropdown botton in the <i>Selection</i> row of the <i>Change Representation</i> panel. The theme color can be seen as a base color that is overpainted by the custom color. Custom colors can be removed for a selection with the 'Clear' option in the color dropdown.
                </HelpText>
            </HelpGroup>
            <HelpGroup header='Representations'>
                <HelpText>
                    Structures can be shown with many different representations (e.g. cartoon or spacefill). The <i>Change Representation</i> panel offers a collection of predefined styles which can be applied using the <i>Preset</i> dropdown button. Additionally any selection atoms, residues or chains can by shown with a custom representation. For that, first select the parts of the structure to be mofified (see help section on <i>Selections</i>) and, second, choose a representation to hide or show from the <i>Show</i> and <i>Hide</i> dropdown bottons in the <i>Selection</i> row of the <i>Change Representation</i> panel. The <i>Everything</i> row applies the action to the whole structure instead of the current selection.
                </HelpText>
            </HelpGroup>
            <HelpGroup header='Surroundings'>
                <HelpText>
                    To show the surroundings of a residue or ligand, click it in the 3D scene or in the sequence widget using {focusTriggers}.
                </HelpText>
            </HelpGroup>

            <HelpSection header='How-to Guides' />
            <HelpGroup header='Create an Image'>
                <HelpText>
                    <p>Use the <Icon svg={CameraSvg} /> icon in the viewport to bring up the screenshot controls.</p>
                    <p>To adjust the size of the image, use the <i>Resolution</i> dropdown.</p>
                </HelpText>
            </HelpGroup>

            <HelpSection header='Mouse Controls' />
            <ViewportHelpContent />
        </div>;
    }
}