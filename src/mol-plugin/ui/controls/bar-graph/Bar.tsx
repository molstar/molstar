import * as React from "react";

export default class Bar extends React.Component<any, any> {
  private value: number;
  constructor(props: any) {
    super(props);

    this.value = this.props.value;
  }

  private onHover = (event: any) => {
    if(!this.props.onHover) { return; }
    this.props.onHover(this.value);
  };

  private offHover = (event: any) => {
    if(!this.props.onMouseLeave) { return; }
    this.props.onMouseLeave();
  }

  private onClick = (event: any) => {
    if(!this.props.onClick) { return; }
    this.props.onClick(this.value);
  }

  private onFocus = () => {
    if(!this.props.onHover) { return; }
    this.props.onHover(this.value);
    if(!this.props.onClick) { return; }
    this.props.onClick(this.value);
  }

  public render() {
    return (
      <rect
        key={this.value}
        className={this.props.class}
        tabIndex={0}
        onFocus={this.onFocus}
        x={this.props.x}
        y={this.props.y}
        width={this.props.width}
        height={this.props.height}
        onMouseOver={this.onHover}
        onClick={this.onClick}
        onMouseLeave={this.offHover}
      />
    );
  }
}
