import React from "react";
import ReactDOM from "react-dom";
import "bootstrap";

import "./bootstrap.min.css";


function Link(props) {
  return (<li className={props.active ? "nav-item active" : "nav-item"}>
    <a className="nav-link" href="#">{props.text}<span class="sr-only">(current)</span></a>
  </li>);
}

function Dropdown(props) {
  return (<ul className="navbar-nav ">
    <li className="nav-item dropdown">
      <a
        className="nav-link dropdown-toggle"
        href="#"
        id="navbarDropdown"
        role="button"
        data-toggle="dropdown"
        aria-haspopup="true"
        aria-expanded="false"
      >
        {props.title}
      </a>
      <div className="dropdown-menu" aria-labelledby="navbarDropdown">
        {props.children}
      </div>
    </li>
  </ul>);
}

function Navbar(props) {
  return (<nav className="navbar navbar-expand-lg navbar-dark bg-primary">
    <a className="navbar-brand" href="#">ACME Haplotype Reconstruction</a>
    <div class="collapse navbar-collapse" id="navbarColor01">
      <ul class="navbar-nav mr-auto">
        {props.children}
      </ul>
    </div>
  </nav>);
}

function App() {
  return (<div>
    <Navbar>
      <Link text="Error Correction" />
      <Link text="Read Graph" />
      <Link text="Haplotype Reconstruction" />
    </Navbar>
    <div style={{ maxWidth: 1140 }} className="container-fluid">
    </div>
  </div>);
}

ReactDOM.render(
  <App />,
  document.body.appendChild(document.createElement("div"))
);
