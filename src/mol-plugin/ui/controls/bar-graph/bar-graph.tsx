/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Luna <paulluna0215@gmail.com>
 */
import * as React from "react";
import { Vec2 } from "mol-math/linear-algebra";
import { scale } from "../line-graph/normilization";

import Bar from "./Bar";

interface BarGraphState {
  bars: any[];
  selected: number;
  currentX: number;
  myRef: any;
  hovering: boolean;
}

export default class BarGraph extends React.Component<any, BarGraphState> {
  height: number;
  padding: number;
  constructor(props: any) {
    super(props);
    this.state = {
      bars: [],
      selected: this.props.mean,
      currentX: -1,
      myRef: React.createRef,
      hovering: false,
    };

    this.height = this.props.height;
    this.padding = this.props.padding ? this.props.padding : 0;
  }

  componentDidMount() {
    this.setHeightOfBars();
  }

  private refCallBack = (element: any) => {
    if(element) {
      this.setState({myRef: element});
    }
  }

  private handleHover = (value: any) => {
    this.props.onHover(value);
  };

  private handleClick = (value: number) => {
    this.props.onClick(value);
    this.setState({ selected: value });
  };

  private handleLeave = () => {
    this.props.onMouseLeave();
  };

  private handleMouseMove = (event: any) => {
    const pt = this.state.myRef.createSVGPoint();
    pt.x = event.clientX;
    pt.y = event.clientY;
    let svgP = pt.matrixTransform(this.state.myRef.getScreenCTM().inverse());
    this.setState({ currentX: svgP.x});
  }

  private setHeightOfBars() {
    let bars = [];
    let max = -Infinity;

    for (let i = 0; i < this.props.counts.length; i++) {
      if (this.props.counts[i] > max) {
        max = this.props.counts[i];
      }
    }

    for (let i = 0; i < this.props.bins.length - 1; i++) {
      let position = Vec2.create(0, 0);
      position[0] = scale(
        parseFloat(this.props.bins[i]),
        this.props.bins[0],
        this.props.bins[this.props.bins.length - 1],
        this.padding / 2,
        600 + this.padding / 2
      );
      position[1] = scale(
        this.props.counts[i],
        0,
        max,
        this.padding / 2,
        this.height + this.padding / 2
      );
      position[1] = this.height + this.padding - position[1];

      bars.push(
        <Bar
          key={this.props.bins[i]}
          value={this.props.bins[i]}
          class="bar"
          x={position[0]}
          y={position[1]}
          width={6}
          height={this.height + this.padding - position[1]}
          onHover={this.handleHover}
          onMouseLeave={this.handleLeave}
          onClick={this.handleClick}
        />
      );
      this.setState({ bars });
    }
  }

  public render() {
    return (
      <div>
        <svg
          className="msp-bar-graph"
          viewBox={`0 0 670 ${this.height + this.padding}`}
          onMouseMove={this.handleMouseMove}
          ref={this.refCallBack}
        >
          <g className="cursor">
            <rect
              x={this.state.currentX}
              y={this.padding / 2}
              width={6}
              height={this.height + this.padding / 2}
              fillOpacity={0.4}
            />
          </g>
          <g className="bars">{this.state.bars}</g>
        </svg>
      </div>
    );
  }
}
