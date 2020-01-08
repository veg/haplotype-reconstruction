import React from "react";
import ReactDOM from "react-dom";
import { BrowserRouter, Route, Link, Switch } from "react-router-dom";
import Nav from "react-bootstrap/Nav";
import Navbar from "react-bootstrap/Navbar";
import Container from "react-bootstrap/Container";
import Row from "react-bootstrap/Row";
import Col from "react-bootstrap/Col";
import axios from "axios";

import "./styles.scss";


function Home() {
  return (<Container>
    <Row>
      <Col md={12}>
        <h1>File browser</h1>
      </Col>
    </Row>
  </Container>);
}

function Visualization() {
  return <h1>Visualization</h1>;
}

function App() {
  axios.get('/api/test')
    .then(response => {
      console.log(response);
    }).catch(error => {
      console.log(error);
    });
  return (<BrowserRouter>
    <Navbar bg="dark" variant="dark">
      <Link className="navbar-brand" to="/">ACME Haplotype Reconstruction</Link>
      <Nav className="mr-auto">
        <Link className="nav-link" to="/visualiztion">Visualization</Link>
      </Nav>
    </Navbar>
    <Switch>
      <Route path="/visualiztion">
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
