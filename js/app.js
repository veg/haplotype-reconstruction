import React, { useState, useEffect } from "react";
import ReactDOM from "react-dom";
import { BrowserRouter, Route, Link, Switch } from "react-router-dom";
import Nav from "react-bootstrap/Nav";
import Navbar from "react-bootstrap/Navbar";
import Container from "react-bootstrap/Container";
import Row from "react-bootstrap/Row";
import Col from "react-bootstrap/Col";
import Breadcrumb from "react-bootstrap/Breadcrumb";
import axios from "axios";

import Visualization from "./visualization";

import "./styles.scss";


function BreadCrumbItem(props) {
  return (<li className="breadcrumb-item">
    <a href="#" onClick={props.onClick}>{props.show}</a>
  </li>);
}

function Home() {
  const [path, setPath] = useState("/"),
    [files, setFiles] = useState([]),
    split = Math.ceil(files.length/2),
    renderFile = file => {
      const { name, type } = file;
      if (type == 'directory') {
        return (<div key={name}>
          <a
            href="#"
            onClick={()=>setPath(path+name+"/")}
          >
            {name}
          </a>
        </div>);
      }
      return (<div key={name}>
          <a href={'/api/static'+path+name} key={name}>
          {name}
        </a>
      </div>);
    };
  useEffect(() => {
    axios.post('/api/files', { path: path })
      .then(response => {
        setFiles(response.data.files);
      }).catch(error => {
        console.log(error);
      });
  }, [path]);
  return (<Container>
    <Row>
      <Col md={12}>
        <h2>File browser</h2>
        <Breadcrumb>
          {('output/' + path).split("/")
            .filter(x=>x)
            .map((directory,i) => (<BreadCrumbItem
                onClick={() => {
                  const new_path = path.split("/").slice(0, i+1).join("/") + "/";
                  setPath(new_path)}
                }
                show={directory}
                key={directory}
              >
                <a href="#">{directory}</a>
              </BreadCrumbItem>))
          }
        </Breadcrumb>
      </Col>
      <Col md={6}>
        {files.slice(0, split).map(renderFile)}
      </Col>
      <Col md={6}>
        {files.slice(split).map(renderFile)}
      </Col>
    </Row>
  </Container>);
}


function App() {
  return (<BrowserRouter>
    <Navbar bg="dark" variant="dark">
      <Link className="navbar-brand" to="/">ACME Haplotype Reconstruction</Link>
      <Nav className="mr-auto">
        <Link className="nav-link" to="/visualization">Visualization</Link>
      </Nav>
    </Navbar>
    <Switch>
      <Route path="/visualization">
        <Visualization />
      </Route>
      <Route path="/">
        <Home />
      </Route>
    </Switch>
  </BrowserRouter>)
}

ReactDOM.render(
  <App />,
  document.body.appendChild(document.createElement('div'))
);
