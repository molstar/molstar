/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Taken/adapted from DensityServer (https://github.com/dsehnal/DensityServer)
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import VERSION from './version'
import { LimitsConfig } from '../config';

export function getDocumentation() {
    function detail(i: number) {
       return `<span class='id'>${i}</span><small> (${Math.round(100 *      LimitsConfig.maxOutputSizeInVoxelCountByPrecisionLevel[i] / 1000 / 1000) / 100 }M voxels)</small>`;
    }
    const detailMax = LimitsConfig.maxOutputSizeInVoxelCountByPrecisionLevel.length - 1;

    // TODO get from config
    const dataSource = `Specifies the data source (determined by the experiment method). Currently, <span class='id'>x-ray</span> and <span class='id'>em</span> sources are supported.`;
    const entryId = `Id of the entry. For <span class='id'>x-ray</span>, use PDB ID (i.e. <span class='id'>1cbs</span>) and for <span class='id'>em</span> use EMDB id (i.e. <span class='id'>emd-8116</span>).`;

    return `<!DOCTYPE html>
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta charset="utf-8" />
<link rel='shortcut icon' href='data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAMAAABEpIrGAAAABGdBTUEAALGPC/xhBQAAAAFzUkdCAK7OHOkAAAAnUExURQAAAMIrHrspHr0oH7soILonHrwqH7onILsoHrsoH7soH7woILwpIKgVokoAAAAMdFJOUwAQHzNxWmBHS5XO6jdtAmoAAACZSURBVDjLxZNRCsQgDAVNXmwb9f7nXZEaLRgXloXOhwQdjMYYwpOLw55fBT46KhbOKhmRR2zLcFJQj8UR+HxFgArIF5BKJbEncC6NDEdI5SatBRSDJwGAoiFDONrEJXWYhGMIcRJGCrb1TOtDahfUuQXd10jkFYq0ViIrbUpNcVT6redeC1+b9tH2WLR93Sx2VCzkv/7NjfABxjQHksGB7lAAAAAASUVORK5CYII=' />
<title>VolumeServer (${VERSION})</title>
<style>
html { -ms-text-size-adjust: 100%; -webkit-text-size-adjust: 100%; }
body { margin: 0; font-family: "Helvetica Neue",Helvetica,Arial,sans-serif; font-weight: 300; color: #333; line-height: 1.42857143; font-size: 14px }
.container { padding: 0 15px; max-width: 970px; margin: 0 auto; }
small { font-size: 80% }
h2, h4 { font-weight: 500; line-height: 1.1; }
h2 { color: black; font-size: 24px; }
h4 { font-size: 18px; margin: 20px 0 10px 0 }
h2 small { color: #777; font-weight: 300 }
hr { box-sizing: content-box; height: 0; overflow: visible; }
a { background-color: transparent; -webkit-text-decoration-skip: objects; text-decoration: none }
a:active, a:hover { outline-width: 0; }
a:focus, a:hover { text-decoration: underline; color: #23527c }
.list-unstyled { padding: 0; list-style: none; margin: 0 0 10px 0 }
.cs-docs-query-wrap { padding: 24px 0; border-bottom: 1px solid #eee }
.cs-docs-query-wrap > h2 { margin: 0; color: black; }
.cs-docs-query-wrap > h2 > span { color: #DE4D4E; font-family: Menlo,Monaco,Consolas,"Courier New",monospace; font-size: 90% }
.cs-docs-param-name, .cs-docs-template-link { color: #DE4D4E; font-family: Menlo,Monaco,Consolas,"Courier New",monospace }
table {margin: 0; padding: 0; }
table th { font-weight: bold; border-bottom: none; text-align: left; padding: 6px 12px }
td { padding: 6px 12px }
td:not(:last-child), th:not(:last-child) { border-right: 1px dotted #ccc }
tr:nth-child(even) { background: #f9f9f9 }
span.id  { color: #DE4D4E; font-family: Menlo,Monaco,Consolas,"Courier New",monospace; }
</style>
</head>
<body>
<div class="container">
<div style='text-align: center; margin-top: 24px;'><span style='font-weight: bold; font-size: 16pt'>VolumeServer</span> <span>${VERSION}</span></div>

<div style='text-align: justify; padding: 24px 0; border-bottom: 1px solid #eee'>
  <p>
    <b>VolumeServer</b> is a service for accessing subsets of volumetric density data. It automatically downsamples the data
    depending on the volume of the requested region to reduce the bandwidth requirements and provide near-instant access to even the
    largest data sets.
  </p>
  <p>
    It uses the text based <a href='https://en.wikipedia.org/wiki/Crystallographic_Information_File'>CIF</a> and binary
    <a href='https://github.com/dsehnal/BinaryCIF' style='font-weight: bold'>BinaryCIF</a>
    formats to deliver the data to the client.
    The server support is integrated into the <a href='https://github.com/dsehnal/LiteMol' style='font-weight: bold'>LiteMol Viewer</a>.
  </p>
</div>

<div class="cs-docs-query-wrap">
  <h2>Data Header / Check Availability <span>/&lt;source&gt;/&lt;id&gt;</span><br>
  <small>Returns a JSON response specifying if data is available and the maximum region that can be queried.</small></h2>
  <div id="coordserver-documentation-ambientResidues-body" style="margin: 24px 24px 0 24px">
    <h4>Examples</h4>
    <a href="/VolumeServer/x-ray/1cbs" class="cs-docs-template-link" target="_blank" rel="nofollow">/x-ray/1cbs</a><br>
    <a href="/VolumeServer/em/emd-8116" class="cs-docs-template-link" target="_blank" rel="nofollow">/em/emd-8116</a>
    <h4>Parameters</h4>
    <table cellpadding="0" cellspacing="0" style='width: 100%'>
    <tbody><tr><th style='width: 80px'>Name</th><th>Description</th></tr>
    <tr>
    <td class="cs-docs-param-name">source</td>
    <td>${dataSource}</td>
    </tr>
    <tr>
    <td class="cs-docs-param-name">id</td>
    <td>${entryId}</td>
    </tr>
    </tbody></table>
  </div>
</div>

<div class="cs-docs-query-wrap">
  <h2>Box <span>/&lt;source&gt;/&lt;id&gt;/box/&lt;a,b,c&gt;/&lt;u,v,w&gt;?&lt;optional parameters&gt;</span><br>
  <small>Returns density data inside the specified box for the given entry. For X-ray data, returns 2Fo-Fc and Fo-Fc volumes in a single response.</small></h2>
  <div style="margin: 24px 24px 0 24px">
    <h4>Examples</h4>
    <a href="/VolumeServer/em/emd-8003/box/-2,7,10/4,10,15.5?encoding=cif&space=cartesian" class="cs-docs-template-link" target="_blank" rel="nofollow">/em/emd-8003/box/-2,7,10/4,10,15.5?excoding=cif&space=cartesian</a><br>
    <a href="/VolumeServer/x-ray/1cbs/box/0.1,0.1,0.1/0.23,0.31,0.18?space=fractional" class="cs-docs-template-link" target="_blank" rel="nofollow">/x-ray/1cbs/box/0.1,0.1,0.1/0.23,0.31,0.18?space=fractional</a>
    <h4>Parameters</h4>
    <table cellpadding="0" cellspacing="0" style='width: 100%'>
    <tbody><tr><th style='width: 80px'>Name</th><th>Description</th></tr>
    <tr>
    <td class="cs-docs-param-name">source</td>
    <td>${dataSource}</td>
    </tr>
    <tr>
    <td class="cs-docs-param-name">id</td>
    <td>${entryId}</td>
    </tr>
    <tr>
    <td class="cs-docs-param-name">a,b,c</td>
    <td>Bottom left corner of the query region in Cartesian or fractional coordinates (determined by the <span class='id'>&amp;space</span> query parameter).</td>
    </tr>
    <tr>
    <td class="cs-docs-param-name">u,v,w</td>
    <td>Top right corner of the query region in Cartesian or fractional coordinates (determined by the <span class='id'>&amp;space</span> query parameter).</td>
    </tr>
    <tr>
    <td class="cs-docs-param-name">encoding</td>
    <td>Determines if text based <span class='id'>CIF</span> or binary <span class='id'>BinaryCIF</span> encoding is used. An optional argument, default is <span class='id'>BinaryCIF</span> encoding.</td>
    </tr>
    <tr>
    <td class="cs-docs-param-name">space</td>
    <td>Determines the coordinate space the query is in. Can be <span class='id'>cartesian</span> or <span class='id'>fractional</span>. An optional argument, default values is <span class='id'>cartesian</span>.</td>
    </tr>
    <tr>
      <td class="cs-docs-param-name">detail</td>
      <td>
        Determines the maximum number of voxels the query can return. Possible values are in the range from ${detail(0)} to ${detail(detailMax)}.
        Default value is <span class='id'>0</span>. Note: different detail levels might lead to the same result.
      </td>
    </tr>
    </tbody></table>
  </div>
</div>

<div class="cs-docs-query-wrap">
  <h2>Cell <span>/&lt;source&gt;/&lt;id&gt;/cell?&lt;optional parameters&gt;</span><br>
  <small>Returns (downsampled) volume data for the entire "data cell". For X-ray data, returns unit cell of 2Fo-Fc and Fo-Fc volumes, for EM data returns everything.</small></h2>
  <div style="margin: 24px 24px 0 24px">
    <h4>Example</h4>
    <a href="/VolumeServer/em/emd-8116/cell?detail=1" class="cs-docs-template-link" target="_blank" rel="nofollow">/em/emd-8116/cell?detail=1</a><br>
    <h4>Parameters</h4>
    <table cellpadding="0" cellspacing="0" style='width: 100%'>
    <tbody><tr><th style='width: 80px'>Name</th><th>Description</th></tr>
    <tr>
    <td class="cs-docs-param-name">source</td>
    <td>${dataSource}</td>
    </tr>
    <tr>
    <td class="cs-docs-param-name">id</td>
    <td>${entryId}</td>
    </tr>
    <tr>
    <td class="cs-docs-param-name">encoding</td>
    <td>Determines if text based <span class='id'>CIF</span> or binary <span class='id'>BinaryCIF</span> encoding is used. An optional argument, default is <span class='id'>BinaryCIF</span> encoding.</td>
    </tr>
    <tr>
      <td class="cs-docs-param-name">detail</td>
      <td>
        Determines the maximum number of voxels the query can return. Possible values are in the range from ${detail(0)} to ${detail(detailMax)}.
        Default value is <span class='id'>0</span>. Note: different detail levels might lead to the same result.
      </td>
    </tr>
    </tbody></table>
  </div>
</div>


<div style="color: #999;font-size:smaller;margin: 20px 0; text-align: right">&copy; 2016 &ndash; now, David Sehnal | Node ${process.version}</div>

</body>
</html>`;
}