/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as React from 'react'
import * as ReactDOM from 'react-dom'
import * as Rx from 'rxjs'

import { QueryDefinition, QueryList } from '../../servers/model/server/api'

import './index.html'

interface State {
    query: Rx.BehaviorSubject<QueryDefinition>,
    id: Rx.BehaviorSubject<string>,
    params: Rx.BehaviorSubject<any>,
    isBinary: Rx.BehaviorSubject<boolean>,
    models: Rx.BehaviorSubject<number[]>,
    url: Rx.Subject<string>
}

class Root extends React.Component<{ state: State }, {  }> {
    render() {
        return <div>
            <div>
                Query: <QuerySelect state={this.props.state} />
            </div>
            <div>
                ID: <input type='text' onChange={t => this.props.state.id.next(t.currentTarget.value)} />
            </div>
            <div>
                Params:<br/>
                <QueryParams state={this.props.state} />
            </div>
            <div>
                Model numbers (empty for all): <ModelNums state={this.props.state} />
            </div>
            <div>
                <input type='checkbox' onChange={t => this.props.state.isBinary.next(!!t.currentTarget.checked)} /> Binary
            </div>
            <div>
                Query string:
                <QueryUrl state={this.props.state} />
            </div>
        </div>
    }
}

class QuerySelect extends React.Component<{ state: State }> {
    render() {
        return <select onChange={s => this.props.state.query.next(QueryList[+s.currentTarget.value].definition)}>
            { QueryList.map((q, i) => <option value={i} key={i} selected={i === 1}>{q.definition.niceName}</option>) }
        </select>
    }
}

class QueryParams extends React.Component<{ state: State }, { prms: string }> {
    state = { prms: '' };

    parseParams(str: string) {
        this.setState({ prms: str });
        try {
            const params = JSON.parse(str);
            this.props.state.params.next(params);
        } catch {
            this.props.state.params.next({});
        }
    }

    componentDidMount() {
        this.props.state.query.subscribe(q => this.setState({ prms: formatParams(q) }))
    }

    render() {
        return <textarea style={{height: '300px'}} value={this.state.prms} cols={80} onChange={t => this.parseParams(t.currentTarget.value)} />;
    }
}

class QueryUrl extends React.Component<{ state: State }, { queryString: string }> {
    state = { queryString: '' };

    componentDidMount() {
        this.props.state.url.subscribe(url => this.setState({ queryString: url }))
    }

    render() {
        return <input type='text' value={this.state.queryString} style={{ width: '800px' }} />
    }
}

class ModelNums extends React.Component<{ state: State }> {
    render() {
        return <input type='text' defaultValue='1' style={{ width: '300px' }} onChange={t =>
            this.props.state.models.next(t.currentTarget.value.split(',')
                .map(v => v.trim())
                .filter(v => !!v)
                .map(v => +v)
                )} />
    }
}

const state: State = {
    query: new Rx.BehaviorSubject(QueryList[1].definition),
    id: new Rx.BehaviorSubject('1cbs'),
    params: new Rx.BehaviorSubject({ }),
    isBinary: new Rx.BehaviorSubject<boolean>(false),
    models: new Rx.BehaviorSubject<number[]>([]),
    url: new Rx.Subject()
}

function formatParams(def: QueryDefinition) {
    const prms = Object.create(null);
    for (const p of def.jsonParams) {
        prms[p.name] = p.exampleValues ? p.exampleValues[0] : void 0;
    }
    return JSON.stringify(prms, void 0, 2);
}

function formatUrl() {
    const json = JSON.stringify({
        name: state.query.value.name,
        id: state.id.value,
        modelNums: state.models.value.length ? state.models.value : void 0,
        binary: state.isBinary.value,
        params: state.params.value
    });
    state.url.next(encodeURIComponent(json));
}

Rx.merge(state.query, state.id, state.params, state.isBinary, state.models).subscribe(s => formatUrl());

ReactDOM.render(<Root state={state} />, document.getElementById('app'));
