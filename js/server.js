const express = require('express'),
  fs = require('fs'),
  bodyParser = require('body-parser'),
  config = require('../config.json');

app = express();
app.use(bodyParser());

app.get('/api/test', (req, res) => {
  res.status(200).json({status: "okay"});
});

app.post('/api/files', (req, res) => {
  const files = fs.readdirSync('./output' + req.body.path)
    .map(file => {
      const path = './output' + req.body.path + file;
      return {
        name: file,
        type: fs.lstatSync(path).isDirectory() ? 'directory' : 'file'
      };
    });
  res.status(200).json({files: files});
});

app.use('/api/static', express.static('output'));

app.listen(config.api_port, () => {
  console.log(`Listening on port ${config.api_port}...`);
});
