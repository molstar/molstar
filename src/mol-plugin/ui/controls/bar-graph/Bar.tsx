import * as React from "react";

export default class Bar extends React.Component<any, any> {
  private value: number;
  constructor(props: any) {
    super(props);

    this.value = this.props.value;
  }

  private onHover = (event: any) => {
    this.props.onHover(this.value);
  };

  private onClick = (event: any) => {
    this.props.onClick(this.value);
  }

  private onFocus = () => {
    this.props.onHover(this.value);
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
        onMouseLeave={this.props.onMouseLeave}
      />
    );
  }
}
