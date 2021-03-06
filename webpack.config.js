const HtmlWebpackPlugin = require('html-webpack-plugin'),
  config = require('./config.json');

module.exports = {
  entry: './js/app.js',
  output: {
    publicPath: '/'
  },
  plugins: [
    new HtmlWebpackPlugin({
      title: "ACME Haplotype Reconstruction"
    })
  ],
  module: {
    rules: [
      {
        test: /\.js$/,
        loader: 'babel-loader'
      },
      {
        test: /\.scss$/,
        use: [
          'style-loader',
          'css-loader',
          'sass-loader'
        ] 
      }
    ],
  },
  devServer: {
    host: '0.0.0.0',
    historyApiFallback: true,
    disableHostCheck: true,
    proxy: {
      '/api': `http://localhost:${config.api_port}`,
      changeOrigin: true
    }
  }
}
