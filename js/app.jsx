import React from "react";
import ReactDOM from "react-dom";
import "bootstrap";

import "bootstrap/dist/css/bootstrap.min.css";


function Link(props) {
  return (
    <p>
      header - {props.header}
      to - {props.to}
    </p>
  );
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

function Navbar() {
  return (<nav className="navbar navbar-expand-lg navbar-light bg-light">
    <a className="navbar-brand" href="#">ACME Haplotype Reconstruction</a>
  </nav>);
}

function App() {
  return (<div>
    <Navbar />
    <div style={{ maxWidth: 1140 }} className="container-fluid">
    </div>
  </div>);
}

ReactDOM.render(
  <App />,
  document.body.appendChild(document.createElement("div"))
);
