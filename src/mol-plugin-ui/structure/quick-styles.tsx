/**
 * Copyright (c) 2022-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Andy Turner <agdturner@gmail.com>
 */

import { PostprocessingParams } from '../../mol-canvas3d/passes/postprocessing';
import { PresetStructureRepresentations } from '../../mol-plugin-state/builder/structure/representation-preset';
import { PluginConfig } from '../../mol-plugin/config';
import { PluginContext } from '../../mol-plugin/context';
import { Color } from '../../mol-util/color';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { CollapsableControls, PurePluginUIComponent } from '../base';
import { Button } from '../controls/common';
import { MagicWandSvg } from '../controls/icons';
import { StructureRepresentationPresetProviderRef } from '../../mol-plugin-state/builder/structure/representation';
//import { Representation } from '../../mol-repr/representation';
//import { StructureParams } from '../../mol-repr/structure/params';

/**
 * A React component that provides a collapsible UI section for quick styling.
 */
export class StructureQuickStylesControls extends CollapsableControls {
    /**
     * Returns the default state of the component.
     * @returns {Object} The default state object.
     * @property {boolean} isCollapsed - Indicates whether the section is collapsed (default is `false`).
     * @property {string} header - The header text of the section (default is `'Quick Styles'`).
     * @property {Object} brand - Styling information for the section header.
     * @property {string} brand.accent - The accent color (default is `'gray'`).
     * @property {JSX.Element} brand.svg - An SVG icon for the header (default is `MagicWandSvg`).
     */
    defaultState() {
        return {
            isCollapsed: false,
            header: 'Quick Styles',
            brand: { accent: 'gray' as const, svg: MagicWandSvg }
        };
    }

    /**
     * Renders the controls for the quick styles section.
     * @returns {JSX.Element} A JSX fragment containing the `QuickStyles` component.
     */
    renderControls() {
        return <>
            <QuickStyles />
        </>;
    }
}

/**
 * Type representing the names of different representations.
 * @typedef {'default' | 'cartoon' | 'spacefill' | 'surface'} RepresentationName
 */
type RepresentationName = 'default' | 'cartoon' | 'spacefill' | 'surface';

/**
 * Type representing the names of different styles.
 * @typedef {'default' | 'illustrative' | 'load' | 'save'} StyleName
 */
type StyleName = 'default' | 'illustrative' | 'load' | 'save';

/**
 * Interface representing the state of the QuickStyles component.
 * @interface
 */
interface QuickStylesState {
    /**
     * Indicates whether the component is busy.
     * @type {boolean}
     */
    busy: boolean,

    /**
     * The current style name.
     * @type {StyleName}
     */
    style: StyleName,

    /**
     * The current representation name.
     * @type {RepresentationName}
     */
    representation: RepresentationName
}

/**
 * A React component that provides quick styling options for structures.
 */
export class QuickStyles extends PurePluginUIComponent<{}, QuickStylesState> {

    /**
     * States for QuickStyles.
     */
    state: QuickStylesState = { busy: false, style: 'default', representation: 'default' };

    /**
     * Applies the specified representation to the plugin's structure component.
     * 
     * This method applies the specified representation to the plugin's structure component based on the
     * provided representation name. It sets the component state to busy while the representation is being
     * applied and resets it once the operation is complete.
     * 
     * @param representationName - The name of the representation to be applied.
     * @returns A promise that resolves when the representation has been applied.
     */
    async applyRepresentation(representationName: RepresentationName) {
        this.setState({ busy: true });
        await applyRepresentation(this.plugin, representationName);
        await applyStyle(this.plugin, this.state.style); // reapplying current style is desired because some presets come with weird params (namely spacefill comes with ignoreLight:true)
        this.setState({ busy: false });
    }

    /**
     * Applies the specified style to the plugin's structure component and canvas.
     * 
     * This method applies the specified style to the plugin's structure component and canvas based on the
     * provided style name. It sets the component state to busy while the style is being applied and resets
     * it once the operation is complete.
     * 
     * @param style - The name of the style to be applied.
     * @returns A promise that resolves when the style has been applied.
     */
    async applyStyle(style: StyleName) {
        this.setState({ busy: true });
        await applyStyle(this.plugin, style);
        this.setState({ busy: false, style });
    }

    /**
     * Handles the loading of a style from a file input event.
     * 
     * This method is triggered when a file is selected in an input element. It reads the file content,
     * parses it as JSON, and applies the loaded style. If the file is not valid JSON, it logs an error
     * and sets the component state to not busy.
     * 
     * @param event - The file input change event.
     */
    async loadRepresentationsAndStyleFromFile(event: React.ChangeEvent<HTMLInputElement>) {
        console.log('Loading representations and style from file');
        this.setState({ busy: true });

        const file = event.target.files?.[0];
        if (!file) {
            console.log('No file selected');
            this.setState({ busy: false });
            return;
        }

        const reader = new FileReader();
        reader.onload = async (e) => {
            const content = e.target?.result as string;
            try {
                const data = JSON.parse(content);
                //const { representations, style } = JSON.parse(content);
                await this.applyLoadedRepresentations(data.representations);
                await this.applyLoadedStyle(data.style);
            } catch (error) {
                console.error('Error parsing JSON:', error);
            } finally {
                this.setState({ busy: false });
            }
        };

        reader.onerror = (e) => {
            console.error('Error reading file:', e);
            this.setState({ busy: false });
        };

        reader.readAsText(file);
    };


    /**
     * Applies the loaded representation to the plugin's structure component.
     * 
     * This method updates the plugin's structure component based on the provided representation object.
     * It uses the representation object to resolve the provider and apply it to the plugin's structure
     * component.
     * 
     * @param representations - The representations object containing the properties to be applied.
     * @returns A promise that resolves when the representation has been applied.
     */
    async applyLoadedRepresentations(representations: StructureRepresentationPresetProviderRef) {

        console.log('Applying loaded representations:', representations);

        const { structures } = this.plugin.managers.structure.hierarchy.selection;

        if (structures.length === 0) {
            console.log('No structures selected');
            return;
        }

        console.log('Selected structures:', structures);

        const provider = this.plugin.builders.structure.representation.resolveProvider(representations);

        if (!provider) {
            console.log('Provider not found for representations:', representations);
            return;
        }

        console.log('Resolved provider:', provider);

        try {
            await this.plugin.managers.structure.component.applyPreset(structures, provider);
            console.log('Representation applied successfully');
        } catch (error) {
            console.error('Error applying representation:', error);
        }

    }

    /**
     * Applies the loaded style to the plugin's structure component and canvas.
     * 
     * This method updates the plugin's structure component options and canvas properties based on the
     * provided style object. It handles the `ignoreLight` option for the structure component and
     * postprocessing options (outline, occlusion, shadow) for the canvas.
     * 
     * @param style - The style object containing the properties to be applied.
     * @returns A promise that resolves when the style has been applied.
     */
    async applyLoadedStyle(style: any) {
        // Update the structure component options if ignoreLight is defined in the style
        if (style.ignoreLight !== undefined) {
            await this.plugin.managers.structure.component.setOptions({
                ...this.plugin.managers.structure.component.state.options,
                ignoreLight: style.ignoreLight
            });
        }
        // Update the canvas properties if the canvas3d instance is available
        if (this.plugin.canvas3d) {
            const p = PD.getDefaultValues(PostprocessingParams);
            this.plugin.canvas3d.setProps({
                postprocessing: {
                    outline: style.outline !== undefined ? style.outline : p.outline,
                    occlusion: style.occlusion !== undefined ? style.occlusion : p.occlusion,
                    shadow: style.shadow !== undefined ? style.shadow : p.shadow
                }
            });
        }
    }

    /**
     * Saves the current representation and style to a JSON file.
     * 
     * This method retrieves the current representation and style settings from the plugin's canvas3d and
     * structure component and saves them to a JSON file. It prompts the user to download the file.
     */
    saveRepresentationsAndStyleToFile(representationPreset: any, style: any) {
        const data = { representationPreset, style };
        const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(data, null, 2));
        const downloadAnchorNode = document.createElement('a');
        downloadAnchorNode.setAttribute("href", dataStr);
        downloadAnchorNode.setAttribute("download", "preset_style.json");
        document.body.appendChild(downloadAnchorNode); // required for firefox
        downloadAnchorNode.click();
        downloadAnchorNode.remove();
    }

    /**
     * Renders the quick styles component.
     * @returns {JSX.Element} A JSX fragment containing the quick styles controls.
     */
    render() {
        return <>
            <NoncollapsableGroup title='Apply Representation'>
                <div className='msp-flex-row'>
                    <Button title='Applies default representation (depends on structure size)'
                        onClick={() => this.applyRepresentation('default')} disabled={this.state.busy} >
                        Default
                    </Button>
                    <Button title='Applies cartoon polymer and ball-and-stick ligand representation'
                        onClick={() => this.applyRepresentation('cartoon')} disabled={this.state.busy} >
                        Cartoon
                    </Button>
                    <Button title='Applies spacefill representation'
                        onClick={() => this.applyRepresentation('spacefill')} disabled={this.state.busy} >
                        Spacefill
                    </Button>
                    <Button title='Applies molecular surface representation'
                        onClick={() => this.applyRepresentation('surface')} disabled={this.state.busy} >
                        Surface
                    </Button>
                </div>
            </NoncollapsableGroup>
            <NoncollapsableGroup title='Apply Style'>
                <div className='msp-flex-row'>
                    <Button title='Applies default appearance (no outline, no ignore-light)'
                        onClick={() => this.applyStyle('default')} disabled={this.state.busy} >
                        Default
                    </Button>
                    <Button title='Applies illustrative appearance (outline, ignore-light)'
                        onClick={() => this.applyStyle('illustrative')} disabled={this.state.busy} >
                        Illustrative
                    </Button>
                </div>
            </NoncollapsableGroup>
            <NoncollapsableGroup title='Load/Save Representations And Style'>
                <div className='msp-flex-row'>
                    <input
                        type="file"
                        accept=".json"
                        onChange={(e) => {
                            console.log('File input changed');
                            this.loadRepresentationsAndStyleFromFile(e);
                        }}
                        style={{ display: 'none' }}
                    />
                    <Button
                        title='Prompts user to load representations and style file'
                        disabled={this.state.busy}
                        onClick={() => {
                            console.log('Load button clicked');
                            const fileInput = document.querySelector('input[type="file"]');
                            if (fileInput) {
                                (fileInput as HTMLInputElement).click();
                            }
                        }}
                    >
                        Load
                    </Button>
                    <Button title='Save current representations and style to file'
                        onClick={() => {
                            console.log('Save button clicked');
                            //this.saveStyleToFile();
                            this.saveRepresentationsAndStyleToFile(
                                this.getCurrentRepresentations(),
                                this.getCurrentStyle());
                        }} disabled={this.state.busy} >
                        Save
                    </Button>
                </div>
            </NoncollapsableGroup>
        </>;
    }

    /**
     * Retrieves the current representation settings from the plugin's structure component.
     * @returns {Object} An object representing the current representation settings.
     * @returns {StructureRepresentationRef} return.representation - The representation reference.
     * @returns {StructureParams} return.reprParams - The representation parameters.
     */
    getCurrentRepresentations(): any {
        return {
            representations: this.plugin.managers.structure.component.state.representations
        };
    }

    /**
     * Retrieves the current style settings from the plugin's canvas3d and structure component.
     * @returns {Object} An object representing the current style settings.
     * @returns {boolean} return.ignoreLight - Indicates whether light is ignored for stylized rendering.
     * @returns {Object} return.outline - The outline settings from the postprocessing properties.
     * @returns {Object} return.occlusion - The occlusion settings from the postprocessing properties.
     * @returns {Object} return.shadow - The shadow settings from the postprocessing properties.
     */
    getCurrentStyle(): any {
        const pp = this.plugin.canvas3d?.props.postprocessing;
        return {
            ignoreLight: this.plugin.managers.structure.component.state.options.ignoreLight,
            outline: pp?.outline,
            occlusion: pp?.occlusion,
            shadow: pp?.shadow
        };
    }
}

/**
 * Visually imitates `ControlGroup` but is always expanded.
 */
function NoncollapsableGroup(props: { title: string, children: any }): JSX.Element {
    return <div className='msp-control-group-wrapper'>
        <div className='msp-control-group-header'><div><b>{props.title}</b></div></div>
        {props.children}
    </div>;
}

/**
 * Applies the specified representation to the plugin's structure component.
 */
async function applyRepresentation(plugin: PluginContext, representation: RepresentationName) {
    const { structures } = plugin.managers.structure.hierarchy.selection;

    switch (representation) {
        case 'default':
            const defaultPreset = plugin.config.get(PluginConfig.Structure.DefaultRepresentationPreset) || PresetStructureRepresentations.auto.id;
            const provider = plugin.builders.structure.representation.resolveProvider(defaultPreset);
            await plugin.managers.structure.component.applyPreset(structures, provider);
            break;
        case 'spacefill':
            await plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations.illustrative);
            break;
        case 'cartoon':
            await plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations['polymer-and-ligand']);
            break;
        case 'surface':
            await plugin.managers.structure.component.applyPreset(structures, PresetStructureRepresentations['molecular-surface']);
            break;
    }
}

/**
 * Applies the specified style to the plugin's structure component and canvas.
 */
async function applyStyle(plugin: PluginContext, style: StyleName) {
    if (style === 'default') {
        await plugin.managers.structure.component.setOptions({ ...plugin.managers.structure.component.state.options, ignoreLight: false });

        if (plugin.canvas3d) {
            const p = PD.getDefaultValues(PostprocessingParams);
            plugin.canvas3d.setProps({
                postprocessing: { outline: p.outline, occlusion: p.occlusion, shadow: p.shadow }
            });
        }
    }

    if (style === 'illustrative') {
        await plugin.managers.structure.component.setOptions({ ...plugin.managers.structure.component.state.options, ignoreLight: true });

        if (plugin.canvas3d) {
            const pp = plugin.canvas3d.props.postprocessing;
            plugin.canvas3d.setProps({
                postprocessing: {
                    outline: {
                        name: 'on',
                        params: pp.outline.name === 'on'
                            ? pp.outline.params
                            : {
                                scale: 1,
                                color: Color(0x000000),
                                threshold: 0.33,
                                includeTransparent: true,
                            }
                    },
                    occlusion: {
                        name: 'on',
                        params: pp.occlusion.name === 'on'
                            ? pp.occlusion.params
                            : {
                                multiScale: { name: 'off', params: {} },
                                radius: 5,
                                bias: 0.8,
                                blurKernelSize: 15,
                                blurDepthBias: 0.5,
                                samples: 32,
                                resolutionScale: 1,
                                color: Color(0x000000),
                                transparentThreshold: 0.4,
                            }
                    },
                    shadow: { name: 'off', params: {} },
                }
            });
        }
    }
}
