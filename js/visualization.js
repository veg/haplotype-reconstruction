import React, { useState, useEffect } from "react";
import Alignment from "alignment.js";
import { Route, Link, Switch, useRouteMatch, useParams } from "react-router-dom";
import Container from "react-bootstrap/Container";
import Row from "react-bootstrap/Row";
import Col from "react-bootstrap/Col";
import axios from "axios";


function Dataset() {
  const { dataset } = useParams();
  const [fasta, setFasta] = useState();
  useEffect(() => {
    axios.get(`/api/static/${dataset}/fastp/bowtie2/pol/sc/superreads-cvs.fasta`)
      .then(response => setFasta(response.data))
      .catch(error => console.log(error));
  });
  return (<div>
    <h4>{dataset}</h4>
    <Alignment fasta={fasta} />
  </div>);
}

function DatasetLink(props) {
  return (<div>
    <Link to={"/visualization/" + props.dataset}>{props.dataset}</Link>
  </div>);
}

function AllDatasets() {
  const [datasets, setDatasets] = useState([]),
    split = Math.ceil(datasets.length/2);
  useEffect(() => {
    axios.post('/api/files', { path: "/" })
      .then(response => {
        const discards = ['lanl', 'references', 'truth'];
        setDatasets(response.data.files.map(dataset => dataset.name)
          .filter(dataset => {
            return discards.indexOf(dataset) == -1;
          }));
      }).catch(error => {
        console.log(error);
      });
    }, [0]);
  return (<Container>
    <Row>
      <Col md={12}>
        <h4>Select a dataset.</h4>
      </Col>
      <Col md={6}>
        {datasets.slice(0, split).map(dataset => (
          <DatasetLink dataset={dataset} key={dataset} />
        ))}
      </Col>
      <Col md={6}>
        {datasets.slice(split).map(dataset => (
          <DatasetLink dataset={dataset} key={dataset} />
        ))}
      </Col>
    </Row>
  </Container>)
}

export default function() {
  const match = useRouteMatch();
  return (<Container>
    <Row>
      <Col>
        <h2>Visualization</h2>
        <Switch>
          <Route path={match.path + "/:dataset/"}>
            <Dataset />
          </Route>
          <Route path={match.url}>
            <AllDatasets />
          </Route>
        </Switch>
      </Col>
    </Row>
  </Container>);
}
