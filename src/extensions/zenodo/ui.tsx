/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { DownloadFile } from '../../mol-plugin-state/actions/file';
import { DownloadStructure, LoadTrajectory } from '../../mol-plugin-state/actions/structure';
import { DownloadDensity } from '../../mol-plugin-state/actions/volume';
import { CoordinatesFormatCategory } from '../../mol-plugin-state/formats/coordinates';
import { TopologyFormatCategory } from '../../mol-plugin-state/formats/topology';
import { TrajectoryFormatCategory } from '../../mol-plugin-state/formats/trajectory';
import { VolumeFormatCategory } from '../../mol-plugin-state/formats/volume';
import { CollapsableControls, CollapsableState } from '../../mol-plugin-ui/base';
import { Button } from '../../mol-plugin-ui/controls/common';
import { OpenInBrowserSvg } from '../../mol-plugin-ui/controls/icons';
import { ParameterControls } from '../../mol-plugin-ui/controls/parameters';
import { PluginContext } from '../../mol-plugin/context';
import { formatBytes } from '../../mol-util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';

type ZenodoFile = {
    bucket: string
    checksum: string
    key: string
    links: {
        [key: string]: string
        self: string
    }
    size: number
    type: string
}

type ZenodoRecord = {
    id: number
    conceptdoi: string
    conceptrecid: string
    created: string
    doi: string
    files: ZenodoFile[]
    revision: number
    updated: string
    metadata: {
        title: string
    }
}

interface State {
    busy?: boolean
    recordValues: PD.Values<typeof ZenodoImportParams>
    importValues?: PD.Values<ImportParams>
    importParams?: ImportParams
    record?: ZenodoRecord
    files?: ZenodoFile[]
}

const ZenodoImportParams = {
    record: PD.Text('', { description: 'Zenodo ID.' })
};

function createImportParams(files: ZenodoFile[], plugin: PluginContext) {
    const modelOpts: [string, string][] = [];
    const topologyOpts: [string, string][] = [];
    const coordinatesOpts: [string, string][] = [];
    const volumeOpts: [string, string][] = [];
    const compressedOpts: [string, string][] = [];

    const structureExts = new Map<string, { format: string, isBinary: boolean }>();
    const coordinatesExts = new Map<string, { format: string, isBinary: boolean }>();
    const topologyExts = new Map<string, { format: string, isBinary: boolean }>();
    const volumeExts = new Map<string, { format: string, isBinary: boolean }>();

    for (const { provider: { category, binaryExtensions, stringExtensions }, name } of plugin.dataFormats.list) {
        if (category === TrajectoryFormatCategory) {
            if (binaryExtensions) for (const e of binaryExtensions) structureExts.set(e, { format: name, isBinary: true });
            if (stringExtensions) for (const e of stringExtensions) structureExts.set(e, { format: name, isBinary: false });
        } else if (category === VolumeFormatCategory) {
            if (binaryExtensions) for (const e of binaryExtensions) volumeExts.set(e, { format: name, isBinary: true });
            if (stringExtensions) for (const e of stringExtensions) volumeExts.set(e, { format: name, isBinary: false });
        } else if (category === CoordinatesFormatCategory) {
            if (binaryExtensions) for (const e of binaryExtensions) coordinatesExts.set(e, { format: name, isBinary: true });
            if (stringExtensions) for (const e of stringExtensions) coordinatesExts.set(e, { format: name, isBinary: false });
        } else if (category === TopologyFormatCategory) {
            if (binaryExtensions) for (const e of binaryExtensions) topologyExts.set(e, { format: name, isBinary: true });
            if (stringExtensions) for (const e of stringExtensions) topologyExts.set(e, { format: name, isBinary: false });
        }
    }

    for (const file of files) {
        const label = `${file.key} (${formatBytes(file.size)})`;
        if (structureExts.has(file.type)) {
            const { format, isBinary } = structureExts.get(file.type)!;
            modelOpts.push([`${file.links.self}|${format}|${isBinary}`, label]);
            topologyOpts.push([`${file.links.self}|${format}|${isBinary}`, label]);
        } else if (volumeExts.has(file.type)) {
            const { format, isBinary } = volumeExts.get(file.type)!;
            volumeOpts.push([`${file.links.self}|${format}|${isBinary}`, label]);
        } else if (topologyExts.has(file.type)) {
            const { format, isBinary } = topologyExts.get(file.type)!;
            topologyOpts.push([`${file.links.self}|${format}|${isBinary}`, label]);
        } else if (coordinatesExts.has(file.type)) {
            const { format, isBinary } = coordinatesExts.get(file.type)!;
            coordinatesOpts.push([`${file.links.self}|${format}|${isBinary}`, label]);
        } else if (file.type === 'zip') {
            compressedOpts.push([`${file.links.self}|${file.type}|true`, label]);
        }
    }

    const params: PD.Params = {};
    let defaultType = '';

    if (modelOpts.length) {
        defaultType = 'structure';
        params.structure = PD.Select(modelOpts[0][0], modelOpts);
    }

    if (topologyOpts.length && coordinatesOpts.length) {
        if (!defaultType) defaultType = 'trajectory';
        params.trajectory = PD.Group({
            topology: PD.Select(topologyOpts[0][0], topologyOpts),
            coordinates: PD.Select(coordinatesOpts[0][0], coordinatesOpts),
        }, { isFlat: true });
    }

    if (volumeOpts.length) {
        if (!defaultType) defaultType = 'volume';
        params.volume = PD.Select(volumeOpts[0][0], volumeOpts);
    }

    if (compressedOpts.length) {
        if (!defaultType) defaultType = 'compressed';
        params.compressed = PD.Select(compressedOpts[0][0], compressedOpts);
    }

    return {
        type: PD.MappedStatic(defaultType, Object.keys(params).length ? params : { '': PD.EmptyGroup() })
    };
}
type ImportParams = ReturnType<typeof createImportParams>

export class ZenodoImportUI extends CollapsableControls<{}, State> {
    protected defaultState(): State & CollapsableState {
        return {
            header: 'Zenodo Import',
            isCollapsed: true,
            brand: { accent: 'cyan', svg: OpenInBrowserSvg },
            recordValues: PD.getDefaultValues(ZenodoImportParams),
            importValues: undefined,
            importParams: undefined,
            record: undefined,
            files: undefined,
        };
    }

    private recordParamsOnChange = (values: any) => {
        this.setState({ recordValues: values });
    };

    private importParamsOnChange = (values: any) => {
        this.setState({ importValues: values });
    };

    private loadRecord = async () => {
        try {
            this.setState({ busy: true });
            const record: ZenodoRecord = await this.plugin.runTask(this.plugin.fetch({ url: `https://zenodo.org/api/records/${this.state.recordValues.record}`, type: 'json' }));
            const importParams = createImportParams(record.files, this.plugin);
            this.setState({
                record,
                files: record.files,
                busy: false,
                importValues: PD.getDefaultValues(importParams),
                importParams
            });
        } catch (e) {
            console.error(e);
            this.plugin.log.error(`Failed to load Zenodo record '${this.state.recordValues.record}'`);
            this.setState({ busy: false });
        }
    };

    private loadFile = async (values: PD.Values<ImportParams>) => {
        try {
            this.setState({ busy: true });

            const t = values.type;
            if (t.name === 'structure') {
                const defaultParams = DownloadStructure.createDefaultParams(this.plugin.state.data.root.obj!, this.plugin);

                const [url, format, isBinary] = t.params.split('|');

                await this.plugin.runTask(this.plugin.state.data.applyAction(DownloadStructure, {
                    source: {
                        name: 'url',
                        params: {
                            url,
                            format: format as any,
                            isBinary: isBinary === 'true',
                            options: defaultParams.source.params.options,
                        }
                    }
                }));
            } else if (t.name === 'trajectory') {
                const [topologyUrl, topologyFormat, topologyIsBinary] = t.params.topology.split('|');
                const [coordinatesUrl, coordinatesFormat] = t.params.coordinates.split('|');

                await this.plugin.runTask(this.plugin.state.data.applyAction(LoadTrajectory, {
                    source: {
                        name: 'url',
                        params: {
                            model: {
                                url: topologyUrl,
                                format: topologyFormat as any,
                                isBinary: topologyIsBinary === 'true',
                            },
                            coordinates: {
                                url: coordinatesUrl,
                                format: coordinatesFormat as any,
                            },
                        }
                    }
                }));
            } else if (t.name === 'volume') {
                const [url, format, isBinary] = t.params.split('|');

                await this.plugin.runTask(this.plugin.state.data.applyAction(DownloadDensity, {
                    source: {
                        name: 'url',
                        params: {
                            url,
                            format: format as any,
                            isBinary: isBinary === 'true',
                        }
                    }
                }));
            } else if (t.name === 'compressed') {
                const [url, format, isBinary] = t.params.split('|');

                await this.plugin.runTask(this.plugin.state.data.applyAction(DownloadFile, {
                    url,
                    format: format as any,
                    isBinary: isBinary === 'true',
                    visuals: true
                }));
            }
        } catch (e) {
            console.error(e);
            this.plugin.log.error(`Failed to load Zenodo file`);
        } finally {
            this.setState({ busy: false });
        }
    };

    private clearRecord = () => {
        this.setState({
            importValues: undefined,
            importParams: undefined,
            record: undefined,
            files: undefined
        });
    };

    private renderLoadRecord() {
        return <div style={{ marginBottom: 10 }}>
            <ParameterControls params={ZenodoImportParams} values={this.state.recordValues} onChangeValues={this.recordParamsOnChange} isDisabled={this.state.busy} />
            <Button onClick={this.loadRecord} style={{ marginTop: 1 }} disabled={this.state.busy || !this.state.recordValues.record}>
                Load Record
            </Button>
        </div>;
    }

    private renderRecordInfo(record: ZenodoRecord) {
        return <div style={{ marginBottom: 10 }}>
            <div className='msp-help-text'>
                <div>Record {`${record.id}`}: <i>{`${record.metadata.title}`}</i></div>
            </div>
            <Button onClick={this.clearRecord} style={{ marginTop: 1 }} disabled={this.state.busy}>
                Clear
            </Button>
        </div>;
    }

    private renderImportFile(params: ImportParams, values: PD.Values<ImportParams>) {
        return values.type.name ? <div style={{ marginBottom: 10 }}>
            <ParameterControls params={params} values={this.state.importValues} onChangeValues={this.importParamsOnChange} isDisabled={this.state.busy} />
            <Button onClick={() => this.loadFile(values)} style={{ marginTop: 1 }} disabled={this.state.busy}>
                Import File
            </Button>
        </div> : <div className='msp-help-text' style={{ marginBottom: 10 }}>
            <div>No supported files</div>
        </div>;
    }

    protected renderControls(): JSX.Element | null {
        return <>
            {!this.state.record ? this.renderLoadRecord() : null}
            {this.state.record ? this.renderRecordInfo(this.state.record) : null}
            {this.state.importParams && this.state.importValues ? this.renderImportFile(this.state.importParams, this.state.importValues) : null}
        </>;
    }
}
