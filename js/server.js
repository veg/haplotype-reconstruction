const express = require('express');

app = express();

app.get('/api/test', (req, res) => {
  res.status(200).json({status: "okay"});
});

app.use('/api/static', express.static('output'));

app.listen(8000, () => {
  console.log('Listening on port 8000...');
});
