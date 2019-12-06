/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as React from 'react';
import { Binding } from '../../mol-util/binding';
import { PluginUIComponent } from '../base';
import { StateTransformer, StateSelection } from '../../mol-state';
import { SelectLoci } from '../../mol-plugin/behavior/dynamic/representation';
import { StructureRepresentationInteraction } from '../../mol-plugin/behavior/dynamic/selection/structure-representation-interaction';
import { Icon } from '../controls/common';

function getBindingsList(bindings: { [k: string]: Binding }) {
    return Object.keys(bindings).map(k => [k, bindings[k]] as [string, Binding])
}

class BindingsHelp extends React.Component<{ bindings: { [k: string]: Binding } }, { isExpanded: boolean }> {
    getBindingComponents() {
        const bindingsList = getBindingsList(this.props.bindings)
        return <ul style={{ paddingLeft: '20px' }}>
            {bindingsList.map(value => {
                const [name, binding] = value
                return !Binding.isEmpty(binding)
                    ? <li key={name}>{Binding.format(binding, name)}</li>
                    : null
            })}
        </ul>
    }

    render() {
        return <HelpText>{this.getBindingComponents()}</HelpText>
    }
}

class HelpText extends React.PureComponent {
    render() {
        return <div className='msp-control-row msp-help-text'>
            <div>{this.props.children}</div>
        </div>
    }
}

class HelpGroup extends React.Component<{ header: string, initiallyExpanded?: boolean }, { isExpanded: boolean }> {
    state = {
        header: this.props.header,
        isExpanded: !!this.props.initiallyExpanded
    }

    toggleExpanded = () => this.setState({ isExpanded: !this.state.isExpanded });

    render() {
        return <div className='msp-control-group-wrapper'>
            <div className='msp-control-group-header'>
                <button className='msp-btn msp-btn-block' onClick={this.toggleExpanded}>
                    <span className={`msp-icon msp-icon-${this.state.isExpanded ? 'collapse' : 'expand'}`} />
                    {this.props.header}
                </button>
            </div>
            {this.state.isExpanded && <div className='msp-control-offset' style={{ display: this.state.isExpanded ? 'block' : 'none' }}>
                {this.props.children}
            </div>}
        </div>
    }
}

function HelpSection(props: { header: string }) {
    return <div className='msp-simple-help-section'>{props.header}</div>;
}

export class HelpContent extends PluginUIComponent {
    componentDidMount() {
        this.subscribe(this.plugin.events.canvas3d.settingsUpdated, () => this.forceUpdate());
    }

    private getMouseBindingComponents() {
        const interactionBindings: { [k: string]: Binding } = {}
        this.plugin.spec.behaviors.forEach(b => {
            const { bindings } = b.defaultParams
            if (bindings) Object.assign(interactionBindings, bindings)
        })
        return <>
            {this.plugin.canvas3d && <HelpGroup key='trackball' header='Moving in 3D'>
                <BindingsHelp bindings={this.plugin.canvas3d.props.trackball.bindings} />
            </HelpGroup>}
            <HelpGroup key='interactions' header='Select, Highlight, Focus'>
                <BindingsHelp bindings={interactionBindings} />
            </HelpGroup>
        </>
    }

    private formatTriggers(binding: Binding) {
        return binding.triggers.map(t => Binding.Trigger.format(t)).join(' or ')
    }

    private getTriggerFor(transformer: StateTransformer, name: string) {
        const state = this.plugin.state.behaviorState
        const selections = state.select(StateSelection.Generators.ofTransformer(transformer))
        const params = selections.length === 1 ? selections[0].params : undefined
        const bindings = params ? params.values.bindings : {}
        const binding: Binding = name in bindings ? bindings[name] : Binding.Empty
        return this.formatTriggers(binding)
    }

    render() {
        const selectToggleTriggers = this.getTriggerFor(SelectLoci, 'clickSelectToggle')
        const structureInteractionTriggers = this.getTriggerFor(StructureRepresentationInteraction, 'clickInteractionAroundOnly')
        // const volumeAroundTriggers = this.getTriggerFor(StructureRepresentationInteraction, 'clickInteractionAroundOnly') // TODO get from correct behavior transform

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
                    To show the surroundings of a residue or ligand, click it in the 3D scene or in the sequence widget using {structureInteractionTriggers}.
                </HelpText>
            </HelpGroup>
            {/* <HelpGroup header='Densities'>
                <HelpText>
                    Densities can be shown for both X-ray and cryo-EM structures. By default the density around an element/atom can be shown by clicking using {volumeAroundTriggers}. The <i>Density Controls</i> panel offers a variety of options to adjust the display of density maps. The absence of the <i>Density Controls</i> panel indicates that no density is available for the loaded entry which is the case for e.g. NMR structures or very old X-ray structures.
                </HelpText>
            </HelpGroup> */}

            <HelpSection header='How-to Guides' />
            {/* <HelpGroup header='RCSB Molecule of the Month Style'>
                <HelpText>
                    <ol style={{ paddingLeft: '20px' }}>
                        <li>First, hide everything, then show everything with the spacefill representation using the <i>Representation</i> panel.</li>
                        <li>Change color theme of the spacefill representation to <i>illustrative</i> using the <i>Structure Settings</i> panel.</li>
                        <li>Set render style to <i>toon</i> and activate <i>occlusion</i> in the <i>General Settings</i> panel.</li>
                    </ol>
                </HelpText>
            </HelpGroup> */}
            <HelpGroup header='Create an Image'>
                <HelpText>
                    <p>Use the <Icon name='screenshot' /> icon in the viewport or go to the <i>Create Image</i> panel and click <i>download</i> to get the same image you see on the 3D canvas.</p>
                    <p>To adjust the size of the image, select <i>Custom</i> from the <i>Size</i> dropdown in the <i>Create Image</i> panel. Adjust the <i>Width</i> and <i>Height</i> using the sliders. To see an image preview with the correct aspect ratio, activate the preview by expanding the <i>Preview</i> panel.</p>
                </HelpText>
            </HelpGroup>

            <HelpSection header='Mouse Controls' />
            {this.getMouseBindingComponents()}
        </div>
    }
}