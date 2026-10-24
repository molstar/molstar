/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StateObjectSelector, StateTransformer, StateTransform } from '@molstar/core/state';
import type { PluginStateObject as SO } from '@molstar/plugin/state/objects';
import type { PluginContext } from '@molstar/plugin/context';
import { Download, ReadFile, DownloadBlob, RawData } from '@molstar/plugin/state/transforms/data/fetch';
import { type FileNameInfo, getFileNameInfo } from '@molstar/core/util/file-info';

export class DataBuilder {
  private get dataState() {
    return this.plugin.state.data;
  }

  rawData(params: StateTransformer.Params<RawData>, options?: Partial<StateTransform.Options>) {
    const data = this.dataState.build().toRoot().apply(RawData, params, options);
    return data.commit({ revertOnError: true });
  }

  // Response Data unless `asResponse: false`, which gives String or Binary Data as `isBinary` says (`.v6/data-response`)
  download(
    params: StateTransformer.Params<Download> & { asResponse: false },
    options?: Partial<StateTransform.Options>,
  ): Promise<StateObjectSelector<SO.Data.String | SO.Data.Binary>>;
  download(
    params: StateTransformer.Params<Download> & { asResponse?: true },
    options?: Partial<StateTransform.Options>,
  ): Promise<StateObjectSelector<SO.Data.Response>>;
  download(
    params: StateTransformer.Params<Download>,
    options?: Partial<StateTransform.Options>,
  ): Promise<StateObjectSelector<SO.Data.String | SO.Data.Binary | SO.Data.Response>>;
  download(params: StateTransformer.Params<Download>, options?: Partial<StateTransform.Options>) {
    const data = this.dataState.build().toRoot().apply(Download, params, options);
    return data.commit({ revertOnError: true });
  }

  downloadBlob(
    params: StateTransformer.Params<DownloadBlob> & { asResponse: false },
    options?: Partial<StateTransform.Options>,
  ): Promise<StateObjectSelector<SO.Data.Blob>>;
  downloadBlob(
    params: StateTransformer.Params<DownloadBlob> & { asResponse?: true },
    options?: Partial<StateTransform.Options>,
  ): Promise<StateObjectSelector<SO.Data.ResponseBlob>>;
  downloadBlob(
    params: StateTransformer.Params<DownloadBlob>,
    options?: Partial<StateTransform.Options>,
  ): Promise<StateObjectSelector<SO.Data.Blob | SO.Data.ResponseBlob>>;
  downloadBlob(params: StateTransformer.Params<DownloadBlob>, options?: Partial<StateTransform.Options>) {
    const data = this.dataState.build().toRoot().apply(DownloadBlob, params, options);
    return data.commit({ revertOnError: true });
  }

  readFile(
    params: StateTransformer.Params<ReadFile> & { asResponse: false },
    options?: Partial<StateTransform.Options>,
  ): Promise<{ data: StateObjectSelector<SO.Data.String | SO.Data.Binary>; fileInfo: FileNameInfo }>;
  readFile(
    params: StateTransformer.Params<ReadFile> & { asResponse?: true },
    options?: Partial<StateTransform.Options>,
  ): Promise<{ data: StateObjectSelector<SO.Data.Response>; fileInfo: FileNameInfo }>;
  readFile(
    params: StateTransformer.Params<ReadFile>,
    options?: Partial<StateTransform.Options>,
  ): Promise<{ data: StateObjectSelector<SO.Data.String | SO.Data.Binary | SO.Data.Response>; fileInfo: FileNameInfo }>;
  async readFile(params: StateTransformer.Params<ReadFile>, options?: Partial<StateTransform.Options>) {
    const data = await this.dataState.build().toRoot().apply(ReadFile, params, options).commit({ revertOnError: true });
    const fileInfo = getFileNameInfo(params.file?.file?.name ?? '');
    return { data: data, fileInfo };
  }

  constructor(public plugin: PluginContext) {}
}
